/**
 * The 3D globe view (`src/gpu/globe.ts`) against a real, headless device —
 * the one thing `tsc`/`vite build` cannot check, since a WGSL template
 * string is just a string to them. Not a parity suite like `gpu.test.ts`:
 * this view has no CPU reference to compare against, so what these assert is
 * narrower — the scene and particle pipelines compile and build, a frame
 * mid-transition (ortho/perspective and reveal both mid-blend) draws without
 * a validation error, and the tween's own state machine lands where it says
 * it will.
 */

import { describe, it, expect, afterEach } from "vitest";
import { gpuDevice, gpuErrors } from "./gpu";
import { GpuSimulation } from "../src/gpu/sim";
import { GpuParticles } from "../src/gpu/particles";
import { Globe3D } from "../src/gpu/globe";
import { ANNULUS } from "../src/geometry";

const OPT = {
  nr: 16, na: 32, gnr: 33, gna: 64,
  ri: ANNULUS.lo, ro: ANNULUS.hi, Ra: 1e4, dt: 1e-3,
};

const device = await gpuDevice();

afterEach(() => {
  expect(gpuErrors.splice(0).join(" | ")).toBe("");
});

(device ? describe : describe.skip)("globe view smoke test", () => {
  it("scene pass compiles, builds, and draws something", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno");
    globe.setViewport(N);
    const now = performance.now();
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    // One realistic frame in, not `tick(performance.now())` called back to
    // back with `toggle` — real elapsed time between those two calls can be
    // a handful of *microseconds*, well under the scene's own `t > 1e-4`
    // ray self-intersection guard for the cut face the orthographic ray
    // starts essentially on, which the real frame loop never hits (frames
    // are ~16 ms apart, not microseconds) but a synchronous test call can.
    globe.tick(now + 16);

    const tex = device!.createTexture({
      size: [N, N], format: "rgba8unorm",
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
    });
    globe.draw(tex.createView(), null);
    await device!.queue.onSubmittedWorkDone();

    const px = device!.createBuffer({
      size: N * N * 4, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
    });
    const enc = device!.createCommandEncoder();
    enc.copyTextureToBuffer({ texture: tex }, { buffer: px, bytesPerRow: N * 4 }, [N, N]);
    device!.queue.submit([enc.finish()]);
    await px.mapAsync(GPUMapMode.READ);
    const img = new Uint8Array(px.getMappedRange().slice(0));
    px.unmap();
    px.destroy();
    tex.destroy();

    // Not every pixel is the clear colour — something was actually rasterised.
    let distinct = 0;
    for (let i = 4; i < img.length; i += 4) {
      if (img[i] !== img[0] || img[i + 1] !== img[1] || img[i + 2] !== img[2]) { distinct++; }
    }
    expect(distinct).toBeGreaterThan(0);

    globe.destroy();
    sim.destroy();
  });

  it("mid-transition frame draws without validation errors", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "viridis");
    globe.setViewport(N);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 2, panX: 0.1, panY: -0.2 });
    // A frame partway through the tween — the ortho/perspective blend and
    // the reveal fade both sit strictly between their endpoints here.
    globe.tick(performance.now() + 500);

    const tex = device!.createTexture({
      size: [N, N], format: "rgba8unorm",
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
    });
    globe.draw(tex.createView(), null);
    await device!.queue.onSubmittedWorkDone();
    tex.destroy();
    globe.destroy();
    sim.destroy();
  });

  it("particles-in-3D pass compiles, builds, and draws with an attached cloud", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    sim.particles = new GpuParticles(sim, { count: 256 });
    sim.particles.setViewport(N);
    const globe = new Globe3D(sim, "inferno");
    globe.setViewport(N);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    globe.tick(performance.now() + 16);   // see the first test's own note on why not `tick(performance.now())`

    const tex = device!.createTexture({
      size: [N, N], format: "rgba8unorm",
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
    });
    globe.draw(tex.createView(), sim.particles);
    await device!.queue.onSubmittedWorkDone();
    tex.destroy();

    sim.particles.destroy();
    globe.destroy();
    sim.destroy();
  });

  it("full 2D -> 3D -> 2D tween runs cleanly frame by frame", async () => {
    const N = 32;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno");
    globe.setViewport(N);
    const tex = device!.createTexture({
      size: [N, N], format: "rgba8unorm",
      usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
    const view = tex.createView();

    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    const t0 = performance.now();
    for (let i = 0; i <= 10; i++) {
      globe.tick(t0 + i * 150);
      globe.draw(view, null);
    }
    expect(globe.viewMode).toBe("3d");
    expect(globe.inTransition).toBe(false);

    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    for (let i = 0; i <= 10; i++) {
      globe.tick(t0 + 2000 + i * 150);
      globe.draw(view, null);
    }
    expect(globe.viewMode).toBe("2d");
    expect(globe.inTransition).toBe(false);

    await device!.queue.onSubmittedWorkDone();
    tex.destroy();
    globe.destroy();
    sim.destroy();
  });
});
