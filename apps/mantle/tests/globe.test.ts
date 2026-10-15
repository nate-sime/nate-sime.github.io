// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The 3D globe view (`src/gpu/globe.ts`) against a real, headless device —
 * the one thing `tsc`/`vite build` cannot check, since a WGSL template
 * string is just a string to them. Not a parity suite like `gpu.test.ts`:
 * this view has no CPU reference to compare against, so what these assert is
 * narrower — the scene and particle pipelines compile and build, a frame
 * mid-transition (ortho/perspective and reveal both mid-blend) draws without
 * a validation error, the tween's own state machine lands where it says it
 * will, and a handful of pixels read back are what the shading claims they
 * are ("globe view pixels", below).
 */

import { describe, it, expect, afterEach } from "vitest";
import { gpuDevice, gpuErrors } from "./gpu";
import { GpuSimulation } from "../src/gpu/sim";
import { GpuParticles } from "../src/gpu/particles";
import { Globe3D, sunDirection, timeOfDay } from "../src/gpu/globe";
import { DEFAULT_SUN, MILKY_WAY, NOISE3 } from "../src/gpu/wgsl";
import { milkyWayCube } from "../src/gpu/sky";
import { COLORMAPS } from "../src/colormaps";
import { SURFACE_MATERIALS, toSurfaceTexture } from "../src/gpu/surfaceAssets";
import { EARTH, MARS, VENUS } from "../src/planets";
import { ANNULUS } from "../src/geometry";
import type { ProceduralSurfaceSettings } from "../src/ui/presets";

const OPT = {
  nr: 16, na: 32, gnr: 33, gna: 64,
  ri: ANNULUS.lo, ro: ANNULUS.hi, Ra: 1e4, dt: 1e-3,
};

const device = await gpuDevice();
// No network fetch in a headless test — the 1x1 fallback exercises exactly
// the same bind group shape a real image would, just through the procedural
// shading branch (`gp.hasSurface === 0`).
const surface = device ? toSurfaceTexture(device, null) : null;
const material = SURFACE_MATERIALS[EARTH.visual.surface];

afterEach(() => {
  expect(gpuErrors.splice(0).join(" | ")).toBe("");
});

// Pure geometry, so it runs with or without an adapter.
describe("time of day", () => {
  const tilts = [EARTH, VENUS, MARS].map((p) => p.visual.axialTiltDeg);

  it("hour zero is the default sun, and a full day comes back to it", () => {
    for (const tilt of tilts) {
      for (const hours of [0, 24, -24]) {
        const sun = sunDirection(tilt, hours);
        for (let i = 0; i < 3; i++) expect(sun[i]).toBeCloseTo(DEFAULT_SUN[i], 9);
      }
    }
  });

  it("turns the sun without changing its length or its latitude", () => {
    for (const tilt of tilts) {
      const t = tilt * Math.PI / 180;
      const axis = [0, Math.cos(t), -Math.sin(t)];
      const lat = (v: readonly number[]) => v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
      for (const hours of [3, 7.5, -11]) {
        const sun = sunDirection(tilt, hours);
        expect(Math.hypot(...sun)).toBeCloseTo(1, 9);
        expect(lat(sun)).toBeCloseTo(lat(DEFAULT_SUN), 9);
      }
    }
  });

  it("advances local time by exactly the hours asked for", () => {
    for (const tilt of tilts) {
      const t0 = timeOfDay(tilt, 0);
      for (const hours of [1, 5.25, -3, 11.75]) {
        const dt = ((timeOfDay(tilt, hours) - t0) % 24 + 24) % 24;
        expect(dt).toBeCloseTo(((hours % 24) + 24) % 24, 6);
      }
    }
  });
});

(device ? describe : describe.skip)("globe view smoke test", () => {
  it("scene pass compiles, builds, and draws something", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
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
    const globe = new Globe3D(sim, "viridis", surface!, material, EARTH);
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

  it("binds Venus's procedural exterior fallback", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!,
      SURFACE_MATERIALS[VENUS.visual.surface], VENUS);
    globe.setViewport(32);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    globe.tick(performance.now() + 16);
    const tex = device!.createTexture({
      size: [32, 32], format: "rgba8unorm", usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
    globe.draw(tex.createView(), null);
    await device!.queue.onSubmittedWorkDone();
    tex.destroy();
    globe.destroy();
    sim.destroy();
  });

  it("binds Mars's procedural exterior fallback", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!,
      SURFACE_MATERIALS[MARS.visual.surface], MARS);
    globe.setViewport(32);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    globe.tick(performance.now() + 16);
    const tex = device!.createTexture({
      size: [32, 32], format: "rgba8unorm", usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
    globe.draw(tex.createView(), null);
    await device!.queue.onSubmittedWorkDone();
    tex.destroy(); globe.destroy(); sim.destroy();
  });

  it("draws a reader-authored procedural exterior", async () => {
    const generated: ProceduralSurfaceSettings = {
      kind: "procedural", seed: 12345, rockiness: 0.8, terrainScale: 0.6,
      oceanCoverage: 0.45, plantLife: 0.7, iceCaps: 0.2, cloudCover: 0.35,
      atmosphereHue: "#73b8ff", atmosphereDensity: 0.4,
    };
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH, generated);
    globe.setViewport(32);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 });
    globe.tick(performance.now() + 16);
    const tex = device!.createTexture({
      size: [32, 32], format: "rgba8unorm", usage: GPUTextureUsage.RENDER_ATTACHMENT,
    });
    globe.draw(tex.createView(), null);
    await device!.queue.onSubmittedWorkDone();
    tex.destroy(); globe.destroy(); sim.destroy();
  });

  it("particles-in-3D pass compiles, builds, and draws with an attached cloud", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    sim.particles = new GpuParticles(sim, { count: 256 });
    sim.particles.setViewport(N);
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
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
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
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

    // Planet travel closes the thermal wedge without changing this settled
    // 3-D camera, then opens it again after the destination is built.
    globe.setCutaway(false);
    globe.tick(t0 + 1800);
    globe.draw(view, null);
    expect(globe.inTransition).toBe(false);
    globe.setCutaway(true, 300);
    const cutawayStart = performance.now();
    globe.tick(cutawayStart + 150);
    globe.draw(view, null);
    expect(globe.inTransition).toBe(true);
    globe.tick(cutawayStart + 400);
    globe.draw(view, null);
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

/**
 * What a frame actually shows, read back pixel by pixel. Each of these pins
 * one claim the smoke tests above cannot see: the sky is written at all (it
 * once was not — a background depth of 1.0 failed the "less" test against a
 * buffer cleared to 1), a cut face with its shading off is exactly the
 * colour map, closing the wedge leaves the atmosphere whole whatever "cut
 * atmosphere" says, the baked Milky Way is the band it was baked from, and
 * the edges are anti-aliased.
 */
(device ? describe : describe.skip)("globe view pixels", () => {
  const settle = (globe: Globe3D, sim: GpuSimulation, N: number): number => {
    globe.setViewport(N);
    globe.toggle({ halfExtent: sim.halfExtent, zoom: 1, panX: 0, panY: 0 }, true);
    const now = performance.now() + 16;
    globe.tick(now);
    return now;
  };

  const render = async (globe: Globe3D, N: number): Promise<Uint8Array> => {
    globe.setViewport(N);
    const tex = device!.createTexture({
      size: [N, N], format: "rgba8unorm",
      usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
    });
    globe.draw(tex.createView(), null);
    const px = device!.createBuffer({
      size: N * N * 4, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ,
    });
    const enc = device!.createCommandEncoder();
    enc.copyTextureToBuffer({ texture: tex }, { buffer: px, bytesPerRow: N * 4 }, [N, N]);
    device!.queue.submit([enc.finish()]);
    await px.mapAsync(GPUMapMode.READ);
    const img = new Uint8Array(px.getMappedRange().slice(0));
    px.unmap(); px.destroy(); tex.destroy();
    return img;
  };

  /**
   * Globe3D's settled perspective camera (`CAMERA` in `wgsl.ts`): the pixel a
   * world point — or, with `atInfinity`, a world direction — lands on.
   */
  const toPixel = (globe: Globe3D, P: readonly number[], N: number, atInfinity = false): [number, number] => {
    const [az, el, dist] = globe.cameraOrientation;
    const dir = [Math.cos(el) * Math.sin(az), Math.sin(el), Math.cos(el) * Math.cos(az)];
    const fwd = dir.map((v) => -v);
    const right = [-fwd[2], 0, fwd[0]].map((v, _, r) => v / Math.hypot(...r));
    const up = [
      right[1] * fwd[2] - right[2] * fwd[1],
      right[2] * fwd[0] - right[0] * fwd[2],
      right[0] * fwd[1] - right[1] * fwd[0],
    ];
    const rel = atInfinity ? P : P.map((v, i) => v - dir[i] * dist);
    const dot = (a: readonly number[], b: readonly number[]) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const t = Math.tan((46 * Math.PI) / 360);
    const cz = dot(rel, fwd);
    const ndc = [dot(rel, right) / (cz * t), dot(rel, up) / (cz * t)];
    return [Math.floor((ndc[0] + 1) / 2 * N), Math.floor((1 - ndc[1]) / 2 * N)];
  };

  const at = (img: Uint8Array, N: number, x: number, y: number): number[] =>
    [...img.subarray((y * N + x) * 4, (y * N + x) * 4 + 3)];

  it("draws the sky: looking sunward, the sun's disc is in the frame", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
    settle(globe, sim, N);
    // Turn round to face the sun, a little off-axis so the planet is not in
    // front of it; the background would otherwise be the bare clear colour.
    globe.zoomBy(Math.log(14 / 6.5) / 0.0015);
    const az = -0.65 + Math.PI + 0.25, el = -0.67 + 0.12;
    globe.orbit(-(az - 0.55) / 2.4, (el - 0.32) / 2.4);
    globe.tick(performance.now() + 32);
    const img = await render(globe, N);
    const [x, y] = toPixel(globe, DEFAULT_SUN, N, true);
    for (const c of at(img, N, x, y)) expect(c).toBeGreaterThan(230);
    globe.destroy(); sim.destroy();
  });

  it("with face shading off, a cut face is exactly the colour map", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    // A uniform field, so every point of both faces maps to the same colour:
    // inferno's middle stop, at T = 0.5.
    const T = sim.buffer("T");
    device!.queue.writeBuffer(T, 0, new Float32Array(T.size / 4).fill(0.5));
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
    settle(globe, sim, N);
    globe.setLighting({ faceShading: false, cutAtmosphere: true, sunHours: 0 });
    // Mid-mantle on the θ = 0 face: clear of the core's glow and the crust line.
    const rho = OPT.ri + 0.6 * (OPT.ro - OPT.ri);
    const [x, y] = toPixel(globe, [rho * Math.cos(0.6), rho * Math.sin(0.6), 0], N);
    const img = await render(globe, N);
    const expected = COLORMAPS.inferno[2].map((v) => v * 255);
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        at(img, N, x + dx, y + dy).forEach((c, i) => expect(Math.abs(c - expected[i])).toBeLessThanOrEqual(1));
      }
    }
    globe.destroy(); sim.destroy();
  });

  it("with the wedge closed, cutting the atmosphere changes nothing", async () => {
    const N = 64;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const globe = new Globe3D(sim, "inferno", surface!, material, EARTH);
    const now = settle(globe, sim, N);
    globe.setCutaway(false);
    globe.tick(now);
    globe.setLighting({ faceShading: true, cutAtmosphere: true, sunHours: 0 });
    const cut = await render(globe, N);
    globe.setLighting({ faceShading: true, cutAtmosphere: false, sunHours: 0 });
    const whole = await render(globe, N);
    expect(cut).toEqual(whole);
    globe.destroy(); sim.destroy();
  });

  it("the baked Milky Way is the band evaluated directly", async () => {
    // A Fibonacci sphere of directions, each sampled from the cube map and
    // evaluated in closed form. A face baked in the wrong orientation would
    // be off by the band's whole brightness, not by filtering error.
    const M = 8192;
    const sky = milkyWayCube(device!);
    const module = device!.createShaderModule({ code: NOISE3 + MILKY_WAY + /* wgsl */ `
@group(0) @binding(0) var skyTex: texture_cube<f32>;
@group(0) @binding(1) var skySamp: sampler;
@group(0) @binding(2) var<storage, read_write> out: array<vec2f>;
@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) id: vec3u) {
  if (id.x >= ${M}u) { return; }
  let y = 1.0 - 2.0 * (f32(id.x) + 0.5) / ${M}.0;
  let r = sqrt(1.0 - y * y);
  let a = f32(id.x) * 2.39996323;
  let D = vec3f(r * cos(a), y, r * sin(a));
  let direct = milkyWay(D);
  out[id.x] = vec2f(length(textureSampleLevel(skyTex, skySamp, D, 0.0).rgb - direct), length(direct));
}` });
    const pipeline = device!.createComputePipeline({ layout: "auto", compute: { module, entryPoint: "main" } });
    const out = device!.createBuffer({ size: M * 8, usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC });
    const read = device!.createBuffer({ size: M * 8, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
    const enc = device!.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, device!.createBindGroup({ layout: pipeline.getBindGroupLayout(0), entries: [
      { binding: 0, resource: sky.view }, { binding: 1, resource: sky.sampler },
      { binding: 2, resource: { buffer: out } },
    ] }));
    pass.dispatchWorkgroups(M / 64);
    pass.end();
    enc.copyBufferToBuffer(out, 0, read, 0, M * 8);
    device!.queue.submit([enc.finish()]);
    await read.mapAsync(GPUMapMode.READ);
    const v = new Float32Array(read.getMappedRange().slice(0));
    read.unmap(); read.destroy(); out.destroy();
    let err = 0, peak = 0;
    for (let i = 0; i < M; i++) { err = Math.max(err, v[2 * i]); peak = Math.max(peak, v[2 * i + 1]); }
    expect(peak).toBeGreaterThan(0.005);
    expect(err).toBeLessThan(0.08 * peak);
  });

  /**
   * Anti-aliasing, against a 3×-supersampled render of the same frame: each
   * pixel's 3×3 block mean is the reference, and the block's middle sample —
   * the very ray the pixel itself would trace — is what the frame would show
   * without anti-aliasing. Over the pixels the silhouette crosses, the real
   * frame must sit much closer to the reference than that. The wedge is
   * closed and the planet airless, so the silhouette is the only edge, and
   * only a ring of pixels about it is compared: stars are sized in pixels,
   * so the two resolutions draw them differently.
   */
  it("anti-aliases the silhouette", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const airless = { label: "airless", visual: { ...EARTH.visual, atmosphere: undefined } };
    const globe = new Globe3D(sim, "inferno", surface!, material, airless);
    const now = settle(globe, sim, N);
    globe.setCutaway(false);
    globe.tick(now);
    const img = await render(globe, N);
    const fine = await render(globe, 3 * N);
    const lum = (p: number[]) => 0.3 * p[0] + 0.59 * p[1] + 0.11 * p[2];
    // The camera looks at the planet's centre, so the silhouette is a circle
    // about the frame's centre.
    const dist = globe.cameraOrientation[2];
    const limb = Math.tan(Math.asin(OPT.ro / dist)) / Math.tan((46 * Math.PI) / 360) * N / 2;
    let edges = 0, smoothed = 0, aliased = 0;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        if (Math.abs(Math.hypot(x + 0.5 - N / 2, y + 0.5 - N / 2) - limb) > 1.5) continue;
        const block: number[][] = [];
        for (let j = 0; j < 3; j++) for (let i = 0; i < 3; i++) block.push(at(fine, 3 * N, 3 * x + i, 3 * y + j));
        const ls = block.map(lum);
        if (Math.max(...ls) - Math.min(...ls) < 48) continue;
        edges++;
        const ref = [0, 1, 2].map((c) => block.reduce((s, p) => s + p[c], 0) / 9);
        const own = at(img, N, x, y);
        for (let c = 0; c < 3; c++) {
          smoothed += Math.abs(own[c] - ref[c]);
          aliased += Math.abs(block[4][c] - ref[c]);
        }
      }
    }
    expect(edges).toBeGreaterThan(50);
    expect(smoothed).toBeLessThan(0.6 * aliased);
    globe.destroy(); sim.destroy();
  });
});
