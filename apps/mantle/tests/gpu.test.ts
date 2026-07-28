/**
 * GPU parity with the CPU reference, and the runtime controls layered on it.
 *
 * The ladder is built the same way as the rest: each test isolates one kernel
 * against the f64 routine it was ported from, and only the last runs the whole
 * coupled loop. A failure therefore names the stage rather than the symptom.
 *
 * Tolerances are stated as *relative* to the field's own scale, because f32
 * carries ~7 decimal digits and every quantity here spans several orders of
 * magnitude. Absolute thresholds would be meaningless for ψ, whose magnitude
 * scales with Ra.
 */

import { describe, it, expect, afterEach } from "vitest";
import { gpuDevice, gpuErrors, maxDiff } from "./gpu";
import { GpuSimulation } from "../src/gpu/sim";
import { Simulation, buoyancyLoad } from "../src/solver/step";
import { Temperature } from "../src/solver/temperature";
import { ANNULUS, box } from "../src/geometry";
import { clampedAxis, periodicAxis, Field } from "../src/spline";
import {
  quadTable, tabulatedLoad, operatorTables, strainRate, W0, W1,
} from "../src/solver/assembly";
import { VariableStokes } from "../src/solver/stokes";
import { viscosity, gammaFor, meanViscosity, strainScale } from "../src/solver/rheology";
import * as dft from "../src/dft";
import { mat } from "../src/linalg";
import { MESH, PRESETS } from "../src/ui/presets";

const OPT = {
  nr: 16, na: 32, gnr: 33, gna: 64,
  ri: ANNULUS.lo, ro: ANNULUS.hi, Ra: 1e4, dt: 1e-3,
};

/** The CPU twin: same spaces, same fixed dt, no CFL sizing. */
const reference = () => new Simulation({
  nr: OPT.nr, na: OPT.na, gnr: OPT.gnr, gna: OPT.gna,
  geom: ANNULUS, Ra: OPT.Ra,
});

const device = await gpuDevice();

/** Any uncaptured validation error is a failure, whatever the numbers say. */
afterEach(() => {
  expect(gpuErrors.splice(0).join(" | ")).toBe("");
});

/** Draw one frame into an offscreen N×N target and bring the pixels back. */
async function snapshot(sim: GpuSimulation, N: number): Promise<Uint8Array> {
  // copyTextureToBuffer requires bytesPerRow to be a multiple of 256. Violating
  // it is a validation error, not a throw: the copy is dropped and the buffer
  // stays zeroed, so the test reads a blank image and blames the renderer.
  expect((N * 4) % 256).toBe(0);
  const tex = device!.createTexture({
    size: [N, N], format: "rgba8unorm",
    usage: GPUTextureUsage.RENDER_ATTACHMENT | GPUTextureUsage.COPY_SRC,
  });
  sim.draw(tex.createView());
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
  return img;
}

describe("gather-form assembly", () => {
  // Pure f64, no GPU: the tables the GPU is handed must reproduce the reference
  // scatter assembly exactly, so a parity failure downstream cannot be blamed on
  // the refactoring of the loop nest.
  it("reproduces the reference buoyancy load", () => {
    const rAx = clampedAxis(16, OPT.ri, OPT.ro), aAx = periodicAxis(32);
    const temp = new Temperature(ANNULUS, 33, 64);
    temp.reset(0.3, 3);
    const T = (r: number, phi: number) => temp.sample(temp.T, r, phi);

    const rt = quadTable(rAx, W0), at = quadTable(aAx, W1);
    const Tq = new Float64Array(rt.x.length * at.x.length);
    for (let q = 0; q < rt.x.length; q++)
      for (let a = 0; a < at.x.length; a++) Tq[q * at.x.length + a] = T(rt.x[q], at.x[a]);

    const ref = buoyancyLoad(rAx, aAx, T, OPT.Ra);
    const got = tabulatedLoad(rt, at, rAx.n, aAx.n, Tq, OPT.Ra);
    let scale = 0, err = 0;
    for (let i = 0; i < rAx.n; i++) {
      scale = Math.max(scale, ...Array.from(ref[i], Math.abs));
      err = Math.max(err, maxDiff(ref[i], got[i]));
    }
    expect(err / scale).toBeLessThan(1e-13);
  });
});

describe.skipIf(!device)("WebGPU pipeline", () => {
  it("has a Stockham FFT matching the reference DFT", async () => {
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    // b is a real field of nr rows; bRe/bIm are its azimuthal modes.
    const b = await sim.read("b");
    const re = await sim.read("bRe"), im = await sim.read("bIm");
    const rows = mat(OPT.nr, OPT.na);
    for (let i = 0; i < OPT.nr; i++) rows[i].set(b.subarray(i * OPT.na, (i + 1) * OPT.na));
    const ref = dft.forward(rows);

    let scale = 0, err = 0;
    for (let i = 0; i < OPT.nr; i++)
      for (let k = 0; k < OPT.na; k++) {
        scale = Math.max(scale, Math.abs(ref.re[i][k]), Math.abs(ref.im[i][k]));
        err = Math.max(err, Math.abs(ref.re[i][k] - re[i * OPT.na + k]),
          Math.abs(ref.im[i][k] - im[i * OPT.na + k]));
      }
    expect(scale).toBeGreaterThan(0);
    expect(err / scale).toBeLessThan(1e-6);
  });

  // The whole point of the stream-function formulation: exact incompressibility
  // is structural, so it must survive the port to f32 untouched — unlike every
  // other quantity here, this one is not allowed to degrade at all.
  it("keeps the GPU velocity divergence-free", async () => {
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    for (let n = 0; n < 5; n++) sim.step();
    const psi = await sim.read("psi");

    const f = new Field(clampedAxis(OPT.nr, OPT.ri, OPT.ro), periodicAxis(OPT.na));
    for (let i = 0; i < OPT.nr; i++) f.c[i].set(psi.subarray(i * OPT.na, (i + 1) * OPT.na));

    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++) {
      const r = OPT.ri + ((OPT.ro - OPT.ri) * i) / 30;
      for (let j = 0; j < 40; j++) {
        const phi = (2 * Math.PI * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(r, phi)));
        const v = f.velocity(r, phi);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    }
    expect(speed).toBeGreaterThan(1); // convecting, not a null field
    expect(div / speed).toBeLessThan(1e-5);
  });

  it("matches the CPU reference over a fixed short run", async () => {
    const cpu = reference();
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    sim.writeTemperature(cpu.temp.T); // identical initial state, exactly
    const T0 = cpu.temp.T.map((r) => Float64Array.from(r));

    for (let n = 0; n < 25; n++) { cpu.step(OPT.dt); sim.step(); }

    const T = await sim.read("T");
    let err = 0, moved = 0;
    for (let i = 0; i < OPT.gnr; i++) {
      const row = T.subarray(i * OPT.gna, (i + 1) * OPT.gna);
      err = Math.max(err, maxDiff(cpu.temp.T[i], row));
      moved = Math.max(moved, maxDiff(T0[i], cpu.temp.T[i]));
    }

    // Guard against a vacuously easy comparison: the field must have evolved by
    // far more than the two paths differ. T ∈ [0, 1], so these are already
    // relative figures.
    expect(moved).toBeGreaterThan(1e-2);
    // f32 noise enters through the biharmonic solve and is then carried by the
    // flow, so it grows with the run; 25 steps is the "fixed short run" of §13.
    // Measured 1.1e-6 here, and still only 4.2e-6 after 200 steps — the bound is
    // loose enough to absorb a different GPU's rounding, not the scheme drifting.
    expect(err).toBeLessThan(1e-5);
  });

  // The render pass reads the temperature buffer directly, so it is the last
  // link in the "no readback" chain and worth exercising: a frame that never
  // leaves the GPU can only be checked by drawing one and looking at it.
  it("renders the annulus from the temperature buffer", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", OPT);
    const N = 64;
    const img = await snapshot(sim, N);

    const at = (x: number, y: number) => img.subarray((y * N + x) * 4, (y * N + x) * 4 + 3);
    const lum = (p: Uint8Array) => p[0] * 0.30 + p[1] * 0.59 + p[2] * 0.11;
    // The hole and the corners are outside [r_i, r_o] and must be background.
    expect([...at(N / 2, N / 2)]).toEqual([5, 5, 12]);
    expect([...at(0, 0)]).toEqual([5, 5, 12]);
    // Inferno is monotone in lightness, so hot (inner, T = 1) reads brighter
    // than cold (outer, T = 0) — the map is not applied upside down.
    const inner = lum(at(N / 2, 13)); // r ≈ 0.63, just outside the hot radius
    const outer = lum(at(N / 2, 4));  // r ≈ 0.93, just inside the cold one
    expect(inner).toBeGreaterThan(outer + 50);
  });

  it("computes the same Nusselt number as the CPU reduction", async () => {
    const cpu = reference();
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    sim.writeTemperature(cpu.temp.T);
    for (let n = 0; n < 10; n++) { cpu.step(OPT.dt); sim.step(); }

    const stat = await sim.read("stat");
    const ref = cpu.temp.nusselt();
    expect(stat[0]).toBeCloseTo(ref.inner, 3);
    expect(stat[1]).toBeCloseTo(ref.outer, 3);
  });
});

/**
 * The runtime controls. Each of these asserts something the UI would
 * otherwise only *appear* to do: a slider that writes a uniform nothing reads,
 * or a dt change that leaves the diffusion factors stale, both look fine on
 * screen and are wrong.
 */
describe.skipIf(!device)("runtime controls", () => {
  // Stokes is linear in the buoyancy load and Ra is a scalar multiplier on it,
  // so doubling Ra must double ψ *exactly*. That makes this both a check that
  // the uniform write reaches the solve and a check that the solve is linear.
  it("scales ψ linearly with Ra", async () => {
    const t = new Temperature(ANNULUS, OPT.gnr, OPT.gna);
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    sim.writeTemperature(t.T);
    const a = await sim.read("psi");

    sim.Ra = 2 * OPT.Ra;
    sim.writeTemperature(t.T);   // same T, re-solved
    const b = await sim.read("psi");

    let scale = 0, err = 0;
    for (let i = 0; i < a.length; i++) {
      scale = Math.max(scale, Math.abs(a[i]));
      err = Math.max(err, Math.abs(2 * a[i] - b[i]));
    }
    expect(scale).toBeGreaterThan(0);
    expect(err / scale).toBeLessThan(1e-5);
  });

  // The contour scale comes from a GPU reduction over the ψ *coefficients*
  // (convex hull, see psiMaxSource) — so it must equal max|c| exactly, not
  // approximately, and must track the flow rather than sitting at its seed value.
  it("reduces max|ψ| over the coefficients", async () => {
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    for (let n = 0; n < 20; n++) sim.step();
    const psi = await sim.read("psi");
    const stat = await sim.read("stat");
    const ref = psi.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    expect(ref).toBeGreaterThan(0);
    expect(stat[2]).toBeCloseTo(ref, 6);
  });

  // ψ isocontours are the streamlines. They must appear when asked for, and
  // must not leak outside the annulus — the overlay is computed for every pixel
  // (fwidth needs uniform control flow), so "does it respect the mask" is a real
  // question, not a formality.
  it("draws streamlines only inside the annulus", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", { ...OPT, levels: 0 });
    for (let n = 0; n < 40; n++) sim.step();
    const plain = await snapshot(sim, N);

    sim.setStreamlines(12, 1.2);
    const lined = await snapshot(sim, N);

    let changed = 0, outside = 0;
    for (let y = 0; y < N; y++)
      for (let x = 0; x < N; x++) {
        const o = (y * N + x) * 4;
        if (plain[o] === lined[o] && plain[o + 1] === lined[o + 1]) continue;
        changed++;
        // Same screen → domain map as the vertex shader (fill = 0.92).
        const px = (2 * (x + 0.5) / N - 1) * (OPT.ro / 0.92);
        const py = (1 - 2 * (y + 0.5) / N) * (OPT.ro / 0.92);
        const r = Math.hypot(px, py);
        if (r < OPT.ri || r > OPT.ro) outside++;
      }
    expect(changed).toBeGreaterThan(200);   // the overlay is visible
    expect(outside).toBe(0);                // and confined to the domain
  });

  // The mesh overlay draws the discretisation, so the two modes must draw two
  // *different* meshes — ψ has 13×32 elements at this size and T 32×64 cells,
  // and a mode that reached the shader but selected the wrong element count
  // would still produce a plausible-looking lattice. Same mask question as the
  // contours, for the same reason: it is computed for every pixel.
  it("draws each mesh only inside the annulus, and the two differ", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", { ...OPT, levels: 0 });
    for (let n = 0; n < 40; n++) sim.step();
    const plain = await snapshot(sim, N);

    sim.mesh = MESH["ψ elements"];
    const spline = await snapshot(sim, N);
    sim.mesh = MESH["T grid"];
    const grid = await snapshot(sim, N);

    // Same screen → domain map as the vertex shader (fill = 0.92).
    const radius = (x: number, y: number) => Math.hypot(
      (2 * (x + 0.5) / N - 1) * (OPT.ro / 0.92),
      (1 - 2 * (y + 0.5) / N) * (OPT.ro / 0.92));

    let onSpline = 0, onGrid = 0, between = 0, outside = 0;
    for (let y = 0; y < N; y++)
      for (let x = 0; x < N; x++) {
        const o = (y * N + x) * 4;
        const a = plain[o] !== spline[o] || plain[o + 1] !== spline[o + 1];
        const b = plain[o] !== grid[o] || plain[o + 1] !== grid[o + 1];
        if (a) onSpline++;
        if (b) onGrid++;
        if (spline[o] !== grid[o] || spline[o + 1] !== grid[o + 1]) between++;
        if (!a && !b) continue;
        const r = radius(x, y);
        if (r < OPT.ri || r > OPT.ro) outside++;
      }
    expect(onSpline).toBeGreaterThan(200);
    expect(onGrid).toBeGreaterThan(200);
    expect(between).toBeGreaterThan(200);   // the mode picks the mesh, not just "a" mesh
    expect(outside).toBe(0);
  });

  // The overlay is a distance field over the element widths, so the *count* of
  // elements is the whole of what it knows about the discretisation. It comes
  // from the axes rather than from a `nr − 3` written into the shader; this is
  // the check that the two still agree.
  it("carries the element counts the spline axes report", async () => {
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    const p = new Int32Array((await sim.read("params")).buffer);
    expect(p[13]).toBe(clampedAxis(OPT.nr, OPT.ri, OPT.ro).elements().length);
    expect(p[14]).toBe(periodicAxis(OPT.na).elements().length);
  });

  // dt has a table behind it (the Thomas factors) as well as a uniform, so the
  // only convincing test is parity against a CPU reference stepping at the *new*
  // dt: stale factors would diffuse by the old amount and drift immediately.
  it("re-factorises the diffusion operator when dt changes", async () => {
    const dt2 = OPT.dt / 2;
    const cpu = reference();
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    sim.setDt(dt2);
    // dt is read back through the f32 uniform — that *is* the value the GPU
    // stepped with, so the round-trip error belongs here rather than being hidden.
    expect(sim.dt / dt2 - 1).toBeLessThan(1e-6);
    sim.writeTemperature(cpu.temp.T);

    for (let n = 0; n < 20; n++) { cpu.step(dt2); sim.step(); }

    const T = await sim.read("T");
    let err = 0;
    for (let i = 0; i < OPT.gnr; i++)
      err = Math.max(err, maxDiff(cpu.temp.T[i], T.subarray(i * OPT.gna, (i + 1) * OPT.gna)));
    expect(err).toBeLessThan(1e-5);
  });

  it("restores the initial condition on reseed", async () => {
    const sim = GpuSimulation.create(device!, "bgra8unorm", OPT);
    for (let n = 0; n < 30; n++) sim.step();
    sim.reseed(0.05, 4);
    expect(sim.steps).toBe(0);
    expect(sim.time).toBe(0);

    const t = new Temperature(ANNULUS, OPT.gnr, OPT.gna);
    t.reset(0.05, 4);
    const T = await sim.read("T");
    let err = 0;
    for (let i = 0; i < OPT.gnr; i++)
      err = Math.max(err, maxDiff(t.T[i], T.subarray(i * OPT.gna, (i + 1) * OPT.gna)));
    expect(err).toBeLessThan(1e-6);
  });
});

/**
 * The μ(T) tier on the GPU.
 *
 * The load-bearing test is the residual one. Everything else here compares the
 * GPU against something that could, in principle, share a mistake with it; that
 * one takes the GPU's *own* ψ, load and viscosity samples and evaluates the
 * variational form against them in f64 on the CPU. It therefore checks the
 * matrix-free kernels, the preconditioner and the CG scalars at once, and it
 * checks them against the definition of the problem rather than against another
 * implementation of it.
 *
 * Its tolerance is set by f32, not by CG, and cannot be tightened by iterating
 * harder. ψ is *stored* in f32, so it is already ε_f32 from whatever it
 * approximates, and the operator maps that to a residual κ ~ h⁻⁴ times larger:
 * measured 1.5e-4 at the 16×32 used here, rising to 1.4e-2 at 48×128 and flat in
 * the budget at 4, 12 and 30 iterations. Taking an f64-converged ψ and rounding
 * only its storage to f32 reproduces the same floor, which is what identifies
 * the cause. So this test asserts that the GPU solves the right system;
 * *convergence* is asserted in f64 in `rheology.test.ts`, where it means
 * something.
 */
describe.skipIf(!device)("variable-μ tier", () => {
  const GAMMA = gammaFor(1e3);
  const variable = (o: Partial<typeof OPT> & { gamma?: number; iters?: number } = {}) =>
    GpuSimulation.create(device!, "bgra8unorm",
      { ...OPT, variable: true, gamma: GAMMA, iters: 24, ...o });

  const axes = () =>
    [clampedAxis(OPT.nr, OPT.ri, OPT.ro), periodicAxis(OPT.na)] as const;
  const rows = (flat: Float32Array, n: number, m: number) =>
    mat(n, m).map((r, i) => { r.set(flat.subarray(i * m, (i + 1) * m)); return r; });

  it("solves the variable-μ system it claims to", async () => {
    const sim = variable({ iters: 30 });
    for (let n = 0; n < 5; n++) sim.step();

    const [rAx, aAx] = axes();
    const vs = new VariableStokes(rAx, aAx, meanViscosity(ANNULUS, GAMMA));
    // μ from the GPU's own quadrature samples, so the comparison isolates the
    // operator rather than re-testing the bicubic sampler.
    const mu = Float64Array.from(await sim.read("Tq"), (t) => viscosity(t, GAMMA));
    const psi = rows(await sim.read("psi"), OPT.nr, OPT.na);
    const b = rows(await sim.read("b"), OPT.nr, OPT.na);

    const Apsi = vs.apply(psi, mu);
    let res = 0, scale = 0;
    for (let i = 1; i < OPT.nr - 1; i++)
      for (let j = 0; j < OPT.na; j++) {
        res = Math.max(res, Math.abs(Apsi[i][j] - b[i][j]));
        scale = Math.max(scale, Math.abs(b[i][j]));
      }
    expect(scale).toBeGreaterThan(0);   // there is a load to balance at all
    expect(res / scale).toBeLessThan(1e-3);
  });

  // The isoviscous limiting case, on the hardware: at zero activation the
  // Krylov tier must land on the direct tier's answer. The two share no solver
  // code beyond the radial matvec, so this is a real check of the CG wrapper —
  // and of the k = 0 mode, which the two tiers treat differently.
  it("reproduces the constant-μ tier at zero activation", async () => {
    const t = new Temperature(ANNULUS, OPT.gnr, OPT.gna);
    const direct = GpuSimulation.create(device!, "bgra8unorm", OPT);
    direct.writeTemperature(t.T);
    const a = await direct.read("psi");

    const krylov = variable({ gamma: 0, iters: 4 });
    krylov.writeTemperature(t.T);
    const c = await krylov.read("psi");

    let scale = 0, err = 0;
    for (let i = 0; i < a.length; i++) {
      scale = Math.max(scale, Math.abs(a[i]));
      err = Math.max(err, Math.abs(a[i] - c[i]));
    }
    expect(scale).toBeGreaterThan(0);
    expect(err / scale).toBeLessThan(1e-4);
  });

  it("matches the CPU reference over a fixed short run", async () => {
    const cpu = new Simulation({
      nr: OPT.nr, na: OPT.na, gnr: OPT.gnr, gna: OPT.gna,
      geom: ANNULUS, Ra: OPT.Ra,
      variable: true, gamma: GAMMA, iters: 24,
    });
    const sim = variable();
    sim.writeTemperature(cpu.temp.T);   // identical state, and ψ reset on both
    const T0 = cpu.temp.T.map((r) => Float64Array.from(r));

    for (let n = 0; n < 15; n++) { cpu.step(OPT.dt); sim.step(); }

    const T = await sim.read("T");
    let err = 0, moved = 0;
    for (let i = 0; i < OPT.gnr; i++) {
      const row = T.subarray(i * OPT.gna, (i + 1) * OPT.gna);
      err = Math.max(err, maxDiff(cpu.temp.T[i], row));
      moved = Math.max(moved, maxDiff(T0[i], cpu.temp.T[i]));
    }
    expect(moved).toBeGreaterThan(1e-2);   // not a vacuous comparison
    // Looser than the tier-1 bound of 1e-5: an iterative solve stopped at a
    // fixed budget leaves a residual, and the f32 and f64 iterations do not
    // leave the *same* one, so the two ψ differ by more than round-off before
    // the flow carries it. Still two orders below the field's own evolution.
    expect(err).toBeLessThan(1e-4);
  });

  // Structural, and the reason for the whole formulation: exact incompressibility
  // comes from u = ∇×ψ, so neither f32 nor a variable rheology may erode it.
  it("keeps the variable-μ GPU velocity divergence-free", async () => {
    const sim = variable();
    for (let n = 0; n < 5; n++) sim.step();
    const psi = await sim.read("psi");

    const [rAx, aAx] = axes();
    const f = new Field(rAx, aAx);
    for (let i = 0; i < OPT.nr; i++)
      f.c[i].set(psi.subarray(i * OPT.na, (i + 1) * OPT.na));

    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++) {
      const r = OPT.ri + ((OPT.ro - OPT.ri) * i) / 30;
      for (let j = 0; j < 40; j++) {
        const phi = (2 * Math.PI * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(r, phi)));
        const v = f.velocity(r, phi);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    }
    expect(speed).toBeGreaterThan(1);
    expect(div / speed).toBeLessThan(1e-5);
  });

  // The contrast slider re-inverts the μ̄(r) blocks. A stale preconditioner is
  // invisible on screen — CG still converges, just more slowly — so the check is
  // that the *solution* moves, and keeps satisfying its own (new) system.
  it("re-inverts the preconditioner when the contrast changes", async () => {
    const sim = variable({ iters: 30 });
    const before = await sim.read("psi");

    const g2 = gammaFor(1e2);
    sim.setGamma(g2);
    expect(sim.gamma / g2 - 1).toBeLessThan(1e-6);
    sim.reseed(0.05, 4);
    for (let n = 0; n < 5; n++) sim.step();

    const [rAx, aAx] = axes();
    const vs = new VariableStokes(rAx, aAx, meanViscosity(ANNULUS, g2));
    const mu = Float64Array.from(await sim.read("Tq"), (t) => viscosity(t, g2));
    const psi = rows(await sim.read("psi"), OPT.nr, OPT.na);
    const b = rows(await sim.read("b"), OPT.nr, OPT.na);

    const Apsi = vs.apply(psi, mu);
    let res = 0, scale = 0, moved = 0;
    for (let i = 1; i < OPT.nr - 1; i++)
      for (let j = 0; j < OPT.na; j++) {
        res = Math.max(res, Math.abs(Apsi[i][j] - b[i][j]));
        scale = Math.max(scale, Math.abs(b[i][j]));
        moved = Math.max(moved, Math.abs(psi[i][j] - before[i * OPT.na + j]));
      }
    expect(moved).toBeGreaterThan(0);
    expect(res / scale).toBeLessThan(1e-3);
  });
});

/**
 * The power law on the GPU.
 *
 * The two rungs that matter are the two ends of the reduction. **At n = 1** the
 * whole strain-rate apparatus — two workgroup reductions, a `pow`, a clamp —
 * must collapse to exactly the μ(T) law, because that is what lets both laws share
 * one set of pipelines with `n` as a uniform; anything approximate there would
 * make the μ(T) tier a second numerical path. **At n = 3** the μ field itself is
 * compared against the f64 law evaluated on the GPU's own ψ and Tq samples,
 * which checks the deformation kernel, both reductions and the law at once, and
 * checks them against the definition rather than against another implementation.
 */
describe.skipIf(!device)("power-law tier", () => {
  const GAMMA = gammaFor(1e3), N = 3;
  const power = (o: { n?: number; iters?: number; picard?: number } = {}) =>
    GpuSimulation.create(device!, "bgra8unorm",
      { ...OPT, variable: true, gamma: GAMMA, n: N, iters: 24, ...o });

  const axes = () =>
    [clampedAxis(OPT.nr, OPT.ri, OPT.ro), periodicAxis(OPT.na)] as const;
  const rows = (flat: Float32Array, n: number, m: number) =>
    mat(n, m).map((r, i) => { r.set(flat.subarray(i * m, (i + 1) * m)); return r; });

  it("evaluates the same μ field as the f64 law", async () => {
    const sim = power();
    for (let k = 0; k < 5; k++) sim.step();

    // The rheology is time-lagged: the μ a step solves with comes from the ψ the
    // *previous* step left. So ψ is read first and the μ compared against it is
    // the one the following step builds — reading both after the same step would
    // be comparing μ with the ψ it produced rather than the ψ it came from.
    const psi = rows(await sim.read("psi"), OPT.nr, OPT.na);
    sim.step();
    const gpu = await sim.read("mu");
    const Tq = await sim.read("Tq");

    const [rAx, aAx] = axes();
    const e = strainRate(operatorTables(rAx, aAx), psi);
    const { d, g } = strainScale(e);

    let err = 0;
    for (let i = 0; i < gpu.length; i++) {
      const ref = viscosity(Tq[i], GAMMA, (e[i] + d) / g, N);
      err = Math.max(err, Math.abs(gpu[i] - ref) / ref);
    }
    // f32 throughout, including two reductions over ~10⁵ terms and a `pow`; the
    // exponent −2/3 damps the strain rate's own error rather than amplifying it.
    expect(err).toBeLessThan(1e-3);
  });

  it("leaves the strain-rate machinery inert at n = 1", async () => {
    const sim = power({ n: 1, iters: 8 });
    for (let k = 0; k < 3; k++) sim.step();
    const mu = await sim.read("mu");
    const Tq = await sim.read("Tq");

    let err = 0;
    for (let i = 0; i < mu.length; i++)
      err = Math.max(err, Math.abs(mu[i] - viscosity(Tq[i], GAMMA)) / mu[i]);
    // Not "close to" μ(T): the exponent is exactly zero, so `pow` returns exactly
    // one and the clamp is inactive. What is left is the f32 `exp` alone.
    expect(err).toBeLessThan(1e-6);
  });

  // n is a uniform, so it is exactly the kind of control that can be written to
  // a buffer nothing reads. This is the ψ ∝ Ra test's counterpart for n.
  it("takes n as a live uniform", async () => {
    const sim = power({ n: 1 });
    for (let k = 0; k < 3; k++) sim.step();
    const before = await sim.read("mu");
    sim.n = N;
    expect(sim.n).toBe(N);
    sim.step();
    const after = await sim.read("mu");

    let moved = 0;
    for (let i = 0; i < after.length; i++)
      moved = Math.max(moved, Math.abs(after[i] - before[i]) / before[i]);
    expect(moved).toBeGreaterThan(0.1);
  });

  it("matches the CPU reference over a fixed short run", async () => {
    const cpu = new Simulation({
      nr: OPT.nr, na: OPT.na, gnr: OPT.gnr, gna: OPT.gna,
      geom: ANNULUS, Ra: OPT.Ra,
      variable: true, gamma: GAMMA, iters: 24, n: N,
    });
    const sim = power();
    sim.writeTemperature(cpu.temp.T);   // identical state, and ψ reset on both
    const T0 = cpu.temp.T.map((r) => Float64Array.from(r));

    for (let k = 0; k < 15; k++) { cpu.step(OPT.dt); sim.step(); }

    const T = await sim.read("T");
    let err = 0, moved = 0;
    for (let i = 0; i < OPT.gnr; i++) {
      const row = T.subarray(i * OPT.gna, (i + 1) * OPT.gna);
      err = Math.max(err, maxDiff(cpu.temp.T[i], row));
      moved = Math.max(moved, maxDiff(T0[i], cpu.temp.T[i]));
    }
    expect(moved).toBeGreaterThan(1e-2);   // not a vacuous comparison
    // Looser than the μ(T) tier's 1e-4, and necessarily so: μ now depends on ψ,
    // so the f32 and f64 runs solve slightly *different* systems each step and
    // the difference feeds back through the coefficient rather than only through
    // the flow. Still two orders below the field's own evolution.
    expect(err).toBeLessThan(1e-3);
  });

  // The whole formulation's reason for existing, under the most demanding
  // rheology it supports: incompressibility comes from u = ∇×ψ, so neither f32
  // nor a solution-dependent μ may erode it.
  it("keeps the nonlinear GPU velocity divergence-free", async () => {
    const sim = power();
    for (let k = 0; k < 5; k++) sim.step();
    const psi = await sim.read("psi");

    const [rAx, aAx] = axes();
    const f = new Field(rAx, aAx);
    for (let i = 0; i < OPT.nr; i++)
      f.c[i].set(psi.subarray(i * OPT.na, (i + 1) * OPT.na));

    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++) {
      const r = OPT.ri + ((OPT.ro - OPT.ri) * i) / 30;
      for (let j = 0; j < 40; j++) {
        const phi = (2 * Math.PI * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(r, phi)));
        const v = f.velocity(r, phi);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    }
    expect(speed).toBeGreaterThan(1);
    expect(div / speed).toBeLessThan(1e-5);
  });
});

/**
 * Every entry in the resolution list must actually build and run. The table
 * invariants are checked without a GPU in `presets.test.ts`; what that cannot
 * see is a device limit — the inverse table reaches 37 MB at the top of the
 * ladder, and the azimuthal transform's shared memory reaches the 16 KB
 * workgroup ceiling at N_φ = 1024. A preset that trips either fails only when a
 * user selects it.
 */
describe.skipIf(!device)("resolution ladder", () => {
  it.each(Object.entries(PRESETS))("builds and steps %s", async (_name, p) => {
    const sim = GpuSimulation.create(device!, "rgba8unorm",
      { nr: p.nr, na: p.na, gnr: p.gnr, gna: p.gna, geom: ANNULUS, Ra: 2e4, dt: p.dt });
    for (let n = 0; n < 10; n++) sim.step();

    const T = await sim.read("T");
    const psi = await sim.read("psi");
    expect(T.length).toBe(p.gnr * p.gna);
    expect(psi.length).toBe(p.nr * p.na);
    // The monotone limiter and the isothermal boundaries bound T; a broken FFT
    // length or an overflowed table shows up here as NaN or as excursions.
    expect(T.every(Number.isFinite)).toBe(true);
    expect(psi.every(Number.isFinite)).toBe(true);
    // Reduce rather than spread: these arrays reach 10⁵ entries, well past the
    // argument limit of Math.min/max.
    const span = (a: Float32Array, f: (x: number) => number) =>
      a.reduce((m, v) => Math.max(m, f(v)), -Infinity);
    expect(-span(T, (v) => -v)).toBeGreaterThan(-1e-3);
    expect(span(T, (v) => v)).toBeLessThan(1 + 1e-3);
    expect(span(psi, Math.abs)).toBeGreaterThan(0); // flow, not a null field
    sim.destroy();
  });

  /**
   * The Krylov tier at the two ends of the ladder, with the power law on so the
   * heaviest path is the one exercised. It allocates the largest buffers in the
   * simulation after the inverse table — the quadrature-point stresses and the
   * viscosity field reach ~10 MB each at ψ 192×512, on top of the 37 MB of
   * inverses — and it is the only path that binds seven storage buffers to one
   * kernel, against a guaranteed limit of eight. Both are failures a user would
   * find by choosing a preset, so both extremes are built rather than assumed.
   */
  it.each([
    ["coarse · ψ 48×128", PRESETS["coarse · ψ 48×128"]],
    ["finest · ψ 192×512", PRESETS["finest · ψ 192×512"]],
  ] as const)("builds and steps %s with μ(T, ε̇)", async (_name, p) => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", {
      nr: p.nr, na: p.na, gnr: p.gnr, gna: p.gna, geom: ANNULUS,
      Ra: 2e4, dt: p.dt, variable: true, gamma: gammaFor(1e3), iters: 3, n: 3,
    });
    for (let n = 0; n < 3; n++) sim.step();

    const T = await sim.read("T");
    const psi = await sim.read("psi");
    expect(T.every(Number.isFinite)).toBe(true);
    expect(psi.every(Number.isFinite)).toBe(true);
    expect(psi.reduce((m, v) => Math.max(m, Math.abs(v)), 0)).toBeGreaterThan(0);
    sim.destroy();
  });
});

/**
 * The Cartesian box, on the GPU.
 *
 * The metric is *compiled into* the WGSL rather than branched on a uniform (see
 * `wgsl.ts`), so the box's kernels are a second set of shader modules that the
 * annulus suite above never touches. Everything here therefore has to be
 * exercised on the box in its own right: parity against the f64 reference for
 * the coupled loop, the rendering's screen → (x, z) map, and the Krylov tier's
 * seven-buffer operator kernel — each of which contains a metric term that would
 * compile perfectly well while being wrong.
 */
describe.skipIf(!device)("Cartesian box", () => {
  const L = 4;
  const BOX = { nr: 16, na: 32, gnr: 33, gna: 64, geom: box(L), Ra: 1e4, dt: 1e-3 };

  /**
   * Parity with the f64 reference. This is the test that actually pins the box's
   * kernels: the CPU path reaches the metric through `Geometry` and the GPU path
   * through emitted WGSL, so the two agreeing is a statement that the *port* of
   * every metric term is right, not merely that the box runs.
   *
   * **Two checkpoints, answering two different questions.** A mis-ported metric
   * term is *systematic*: it makes the two paths solve different equations, so it
   * shows on the very first step at a size f32 round-off never reaches. Five
   * steps with a tight bound is what catches that.
   *
   * Twenty-five steps is not the same claim, and cannot be. This box convects at
   * max|u| ≈ 56 — nine times the annulus at the same Ra, whose layer is only 0.45
   * deep against this one's 1 — so it turns over more than once in that time, and
   * two trajectories seeded ε apart in f32 separate at the rate the flow mixes.
   * Measured: 1.0e-6 at one step, 8.7e-6 at five, 1.9e-3 at twenty-five, against
   * a field that has changed by O(1). The late bound says they have not
   * *diverged*; it does not say they still agree digit for digit, and a tolerance
   * that claimed otherwise would be measuring the seed, not the port.
   */
  it("matches the CPU reference over a fixed short run", async () => {
    const cpu = new Simulation({
      nr: BOX.nr, na: BOX.na, gnr: BOX.gnr, gna: BOX.gna, geom: box(L), Ra: BOX.Ra,
    });
    const sim = GpuSimulation.create(device!, "rgba8unorm", BOX);
    sim.writeTemperature(cpu.temp.T);
    const T0 = cpu.temp.T.map((r) => Float64Array.from(r));

    const drift = async () => {
      const T = await sim.read("T");
      let e = 0;
      for (let i = 0; i < BOX.gnr; i++)
        e = Math.max(e, maxDiff(cpu.temp.T[i], T.subarray(i * BOX.gna, (i + 1) * BOX.gna)));
      return e;
    };

    for (let k = 0; k < 5; k++) { cpu.step(BOX.dt); sim.step(); }
    expect(await drift()).toBeLessThan(1e-4);   // still round-off, not a scheme

    // Both Nusselt numbers, taken at the *early* checkpoint and against the f64
    // reduction. The GPU's normalisation is `1/L` read out of the uniform block;
    // a box still carrying the annulus' `ln(r_o/r_i)/2π` would give a number of
    // the right order and the wrong value, which the T comparison cannot see.
    //
    // Read here rather than at the end because Nu is a boundary *gradient*: it
    // divides by dr, so it amplifies the trajectory separation below by a factor
    // of ~30 and would be comparing the seed rather than the reduction. `stat`
    // is written inside `step`, after diffusion, so it describes exactly the
    // field just compared.
    const stat = await sim.read("stat");
    const ref = cpu.temp.nusselt();
    expect(ref.outer).not.toBeCloseTo(1, 2);   // not the trivial conductive value
    expect(stat[0] / ref.inner).toBeCloseTo(1, 3);
    expect(stat[1] / ref.outer).toBeCloseTo(1, 3);

    for (let k = 5; k < 25; k++) { cpu.step(BOX.dt); sim.step(); }
    const err = await drift();
    let moved = 0;
    for (let i = 0; i < BOX.gnr; i++)
      moved = Math.max(moved, maxDiff(T0[i], cpu.temp.T[i]));
    expect(moved).toBeGreaterThan(1e-2);   // not a vacuous comparison
    expect(err).toBeLessThan(1e-2);
    expect(err).toBeLessThan(0.05 * moved); // and small against the run itself
    sim.destroy();
  });

  /**
   * Pure conduction reads Nu = 1 at both ends. This is what `nuScale` is *for*,
   * and the one value a wrong normalising constant cannot happen to land on: any
   * factor other than `1/L` scales both numbers off 1 together.
   */
  it("normalises the Nusselt numbers by the box length", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", BOX);
    const t = new Temperature(box(L), BOX.gnr, BOX.gna);
    for (let i = 0; i < BOX.gnr; i++) t.T[i].fill(t.conduction(i));
    sim.writeTemperature(t.T);
    // Conduction is a steady state with no buoyancy forcing (∂_x T ≡ 0, so
    // ψ ≡ 0 and nothing advects), and `stat` is only written by `step` — so one
    // step is how the reduction is reached, and it reaches it on the same field.
    sim.step();
    const stat = await sim.read("stat");
    expect(stat[0]).toBeCloseTo(1, 3);
    expect(stat[1]).toBeCloseTo(1, 3);
    sim.destroy();
  });

  /**
   * The render pass's screen → (x, z) map. The box occupies a rectangle of the
   * square viewport rather than a disk, so the overlay mask is a different
   * predicate from the annulus', and it is computed for every pixel before the
   * in-domain test — the same `fwidth` uniformity constraint, and the same
   * chance of drawing lines across the background.
   */
  it("draws its overlays only inside the rectangle", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", { ...BOX, levels: 0 });
    for (let n = 0; n < 40; n++) sim.step();
    const plain = await snapshot(sim, N);

    sim.setStreamlines(12, 1.2);
    const lined = await snapshot(sim, N);
    sim.mesh = MESH["ψ elements"];
    const meshed = await snapshot(sim, N);

    // Same map as `domainFn` for a box (fill = 0.92): the half-extent is half
    // the longer side over `fill`, and the origin is the rectangle's centre.
    const half = (0.5 * Math.max(L, 1)) / 0.92;
    let changed = 0, outside = 0;
    for (let y = 0; y < N; y++)
      for (let x = 0; x < N; x++) {
        const o = (y * N + x) * 4;
        const moved = plain[o] !== lined[o] || plain[o + 1] !== lined[o + 1]
          || plain[o] !== meshed[o] || plain[o + 1] !== meshed[o + 1];
        if (!moved) continue;
        changed++;
        const bx = (2 * (x + 0.5) / N - 1) * half + 0.5 * L;
        const bz = (1 - 2 * (y + 0.5) / N) * half + 0.5;
        if (bx < 0 || bx > L || bz < 0 || bz > 1) outside++;
      }
    expect(changed).toBeGreaterThan(200);
    expect(outside).toBe(0);
    sim.destroy();
  });

  /**
   * The box is *wider* than it is deep, so a viewport that sized it on the depth
   * would run it off both sides of the frame. Checking that the field occupies a
   * horizontal band and not the whole square is what distinguishes a fitted
   * rectangle from one drawn at the annulus' scale.
   */
  it("fits the rectangle to its longer side", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", BOX);
    for (let n = 0; n < 20; n++) sim.step();
    const px = await snapshot(sim, N);

    // Out-of-domain pixels are the shader's background, (0.02, 0.02, 0.047) — so
    // "not that colour" is the in-domain mask. A brightness threshold would not
    // do: the cold boundary maps to near-black, which is *darker* than the
    // background, and the top of the box would read as outside it.
    const inside = (x: number, y: number) => {
      const o = (y * N + x) * 4;
      return px[o] !== 5 || px[o + 1] !== 5 || px[o + 2] !== 12;
    };
    const mid = N >> 1;
    let row = 0, col = 0;
    for (let i = 0; i < N; i++) {
      if (inside(i, mid)) row++;
      if (inside(mid, i)) col++;
    }
    // The width is the fitted dimension, so it spans `fill` of the frame; the
    // height is 1/L of that. Fitting on the depth instead would fill both, and
    // fitting as though it were the annulus would run it off the sides.
    expect(row).toBeGreaterThan(0.85 * N);
    expect(col).toBeGreaterThan((0.8 * 0.92 * N) / L);
    expect(col).toBeLessThan((1.2 * 0.92 * N) / L);
    sim.destroy();
  });

  /** The Krylov tier's own kernels, on the box's metric. */
  it("runs the variable-μ tier", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", {
      ...BOX, variable: true, gamma: gammaFor(1e3), iters: 8, n: 3,
    });
    for (let n = 0; n < 5; n++) sim.step();
    const psi = await sim.read("psi");
    const T = await sim.read("T");
    expect(psi.every(Number.isFinite)).toBe(true);
    expect(T.every(Number.isFinite)).toBe(true);
    expect(psi.reduce((m, v) => Math.max(m, Math.abs(v)), 0)).toBeGreaterThan(0);

    // And the property the whole formulation exists for, on the second metric:
    // in a box `u = (ψ_x, −ψ_z)`, whose divergence is structurally zero.
    const g = box(L);
    const f = new Field(clampedAxis(BOX.nr, 0, 1), periodicAxis(BOX.na, L), g);
    for (let i = 0; i < BOX.nr; i++)
      f.c[i].set(psi.subarray(i * BOX.na, (i + 1) * BOX.na));
    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++)
      for (let j = 0; j < 40; j++) {
        const z = i / 30, x = (L * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(z, x)));
        const v = f.velocity(z, x);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    expect(speed).toBeGreaterThan(1);
    expect(div / speed).toBeLessThan(1e-5);
    sim.destroy();
  });
});

/**
 * Free-slip side walls, on the GPU.
 *
 * The wall is a projection made in *mode space*, inside the tridiagonal kernel —
 * so it is one emitted branch in one shader, and the only way to see it is to
 * look at what comes out. Two things are worth checking and they are different
 * claims: that the projection is the same one the f64 reference makes (parity),
 * and that it survives the frame loop rather than merely surviving the initial
 * condition.
 */
describe.skipIf(!device)("free-slip side walls", () => {
  const WALLED = {
    nr: 16, na: 32, gnr: 33, gna: 64,
    geom: box(1, "free-slip"), Ra: 1e4, dt: 5e-4,
  };

  /**
   * Parity, on the same two-checkpoint reasoning as the periodic box above: the
   * early bound is what would catch a projection emitted into the wrong branch,
   * the late one only that the two runs have not diverged. Seeded on the mode
   * the walls make critical, so the flow is vigorous and the separation is fast.
   */
  it("matches the CPU reference over a fixed short run", async () => {
    const cpu = new Simulation({
      nr: WALLED.nr, na: WALLED.na, gnr: WALLED.gnr, gna: WALLED.gna,
      geom: box(1, "free-slip"), Ra: WALLED.Ra, seed: { amp: 0.1, mode: 1 },
    });
    const sim = GpuSimulation.create(device!, "rgba8unorm", WALLED);
    sim.writeTemperature(cpu.temp.T);
    const T0 = cpu.temp.T.map((r) => Float64Array.from(r));

    const drift = async () => {
      const T = await sim.read("T");
      let e = 0;
      for (let i = 0; i < WALLED.gnr; i++)
        e = Math.max(e,
          maxDiff(cpu.temp.T[i], T.subarray(i * WALLED.gna, (i + 1) * WALLED.gna)));
      return e;
    };

    for (let k = 0; k < 5; k++) { cpu.step(WALLED.dt); sim.step(); }
    expect(await drift()).toBeLessThan(1e-4);

    for (let k = 5; k < 25; k++) { cpu.step(WALLED.dt); sim.step(); }
    const err = await drift();
    let moved = 0;
    for (let i = 0; i < WALLED.gnr; i++)
      moved = Math.max(moved, maxDiff(T0[i], cpu.temp.T[i]));
    expect(moved).toBeGreaterThan(1e-2);
    expect(err).toBeLessThan(1e-2);
    expect(err).toBeLessThan(0.05 * moved);
    sim.destroy();
  });

  /**
   * The wall, as the frame loop sees it. An antisymmetric field is written in —
   * a state the walls forbid outright — and one step must remove it. A
   * projection applied only where the initial condition is built, or emitted
   * into the wrong branch of the kernel, would let it through; f32 sets the
   * floor at round-off rather than zero, which is what the tolerance is.
   */
  it("annihilates the antisymmetric component every step", async () => {
    const sim = GpuSimulation.create(device!, "rgba8unorm", WALLED);
    const { gnr, gna } = WALLED;
    const t = new Temperature(box(1, "free-slip"), gnr, gna);
    for (let i = 1; i < gnr - 1; i++)
      for (let j = 0; j < gna; j++)
        t.T[i][j] = t.conduction(i)
          + 0.2 * Math.sin((Math.PI * i) / (gnr - 1)) * Math.sin((2 * Math.PI * j) / gna);
    t.applyBC();
    sim.writeTemperature(t.T);

    const odd = async () => {
      const T = await sim.read("T");
      let m = 0;
      for (let i = 0; i < gnr; i++)
        for (let j = 0; j < gna; j++)
          m = Math.max(m, Math.abs(T[i * gna + j] - T[i * gna + ((gna - j) % gna)]));
      return m;
    };
    expect(await odd()).toBeGreaterThan(0.1);   // the seed really is antisymmetric
    sim.step();
    expect(await odd()).toBeLessThan(1e-5);     // and one step is all it survives
    for (let k = 0; k < 20; k++) sim.step();
    expect(await odd()).toBeLessThan(1e-5);     // and it does not creep back
    sim.destroy();
  });

  /**
   * The render draws half the period it solves, so a walled box of width 1 fills
   * the frame where the length-2 period it is computed on would overflow it. The
   * mask is the same predicate the periodic box uses, against a different width —
   * getting the factor of two backwards here is invisible in every other test.
   */
  it("draws the width, not the period", async () => {
    const N = 128;
    const sim = GpuSimulation.create(device!, "rgba8unorm", WALLED);
    for (let n = 0; n < 20; n++) sim.step();
    const px = await snapshot(sim, N);

    const inside = (x: number, y: number) => {
      const o = (y * N + x) * 4;
      return px[o] !== 5 || px[o + 1] !== 5 || px[o + 2] !== 12;
    };
    const mid = N >> 1;
    let row = 0, col = 0;
    for (let i = 0; i < N; i++) {
      if (inside(i, mid)) row++;
      if (inside(mid, i)) col++;
    }
    // Width 1 against depth 1: square, so `fill` of the frame in both directions.
    // Drawing the period instead would halve the height.
    expect(row).toBeGreaterThan(0.85 * N);
    expect(col).toBeGreaterThan(0.85 * N);
    sim.destroy();
  });
});
