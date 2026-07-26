/**
 * The validated CPU pipeline, resident on the GPU.
 *
 * A step is fourteen compute dispatches in one command buffer, and nothing is
 * read back. The buoyancy → Stokes → ψ steps come at the *end*, matching the CPU
 * reference: the constructor solves Stokes once from the initial T, and each
 * step advects, diffuses, then re-solves — so the velocity a step advects with
 * is always the one balancing the T it starts from.
 *
 *   advect  A: T  → TA   (semi-Lagrangian, +dt)
 *           B: TA → TB   (reverse pass, −dt; boundary rows zeroed as on the CPU)
 *           C: T,TB → TA (BFECC correction)
 *           D: TA → TB   (final pass, limited against T at the departure point)
 *   diffuse    TB → modes → Thomas per mode → T
 *   stokes     T  → quadrature samples → gather to the load → DFT →
 *              per-mode A_k⁻¹ matvec → inverse DFT → ψ coefficients
 *
 * Buffers never swap roles, so every bind group is built once at init.
 *
 * **Fixed dt.** The CPU reference sizes dt from the advective CFL, which needs a
 * max-reduction of |u| on the host — a readback in the hot loop. Semi-Lagrangian
 * advection and implicit diffusion are both unconditionally stable, so dt is an
 * accuracy knob, not a stability limit: fixing it costs nothing
 * physical and lets the Thomas factors be built once, in f64, at init.
 */

import { Axis, P, clampedAxis, periodicAxis } from "../spline";
import { modeInverses } from "../solver/operators";
import { operatorTables, quadTable, SLOTS, W0, W1 } from "../solver/assembly";
import { meanViscosity } from "../solver/rheology";
import { Temperature, diffusionFactors } from "../solver/temperature";
import * as w from "./wgsl";

export interface GpuSimOptions {
  nr?: number; na?: number;      // spline space for ψ
  gnr?: number; gna?: number;    // grid for T
  ri?: number; ro?: number;
  Ra?: number; dt?: number;
  fill?: number;                 // fraction of the half-viewport r_o spans
  levels?: number;               // streamline contours across [−ψmax, ψmax]; 0 = off
  lineW?: number;                // contour half-width, in pixels
  mesh?: number;                 // mesh overlay: 0 = off, 1 = ψ elements, 2 = T grid
  /**
   * Tier 2: μ(T) by matrix-free PCG instead of the direct DFT solve. A
   * construction-time choice, not a uniform — the Krylov path allocates the
   * quadrature-point stress buffers, which are the largest in the simulation
   * after the inverse table, so a constant-μ run must not pay for them.
   */
  variable?: boolean;
  gamma?: number;                // ln of the viscosity contrast
  iters?: number;                // Krylov budget per solve
  /**
   * Power-law index. `n = 1` is μ(T) exactly, so the two variable laws are
   * one tier and one set of pipelines — this is a uniform, not a rebuild.
   */
  n?: number;
  picard?: number;               // rheology updates per solve; 1 = pure time-lagging
}

const S = Float32Array.BYTES_PER_ELEMENT;

/**
 * Indices into the uniform block (see `PARAMS` in wgsl.ts). Only these change at
 * runtime; everything else is fixed by the resolution.
 */
const F = { Ra: 6, dt: 7, levels: 11, lineW: 12, gamma: 13, n: 14, mesh: 15 } as const;

export class GpuSimulation {
  readonly nr: number; readonly na: number;
  readonly gnr: number; readonly gna: number;
  time = 0;
  steps = 0;
  /** Krylov budget per solve; free to change (it only alters the dispatch count). */
  iters: number;
  /** Rheology updates per solve; likewise only a dispatch count. */
  picard: number;
  /** Latest diagnostics; refreshed by `pollStats`, stale between polls. */
  stats = { nuInner: NaN, nuOuter: NaN, psiMax: NaN };

  private readonly rAx: Axis;
  private readonly aAx: Axis;
  private readonly buf: Record<string, GPUBuffer> = {};
  private readonly pipe: Record<string, GPUComputePipeline> = {};
  private readonly bind: Record<string, GPUBindGroup> = {};
  private render!: GPURenderPipeline;
  private renderBind!: GPUBindGroup;
  private statPending = false;
  private readonly params = new ArrayBuffer(128);
  private readonly pf: Float32Array;

  private constructor(readonly device: GPUDevice, readonly o: Required<GpuSimOptions>) {
    this.nr = o.nr; this.na = o.na; this.gnr = o.gnr; this.gna = o.gna;
    this.iters = o.iters;
    this.picard = o.picard;
    for (const n of [o.na, o.gna])
      if (!Number.isInteger(Math.log2(n))) throw new Error(`${n} is not a power of two`);

    const rAx = clampedAxis(o.nr, o.ri, o.ro), aAx = periodicAxis(o.na);
    this.rAx = rAx; this.aAx = aAx;
    // The buoyancy load's two tables, `w·B` and `w·N′`. `ot` adds the six the
    // variable-μ operator gathers against, and is built only in that tier — the
    // constant path would otherwise construct twelve tables to upload none of
    // them. Its eval half is unused here either way: the GPU recomputes the
    // basis per quadrature point rather than tabulating it (see `qevalSource`).
    const rt = quadTable(rAx, W0), at = quadTable(aAx, W1);
    const ot = o.variable ? operatorTables(rAx, aAx) : null;
    const temp = new Temperature(o.gnr, o.gna, o.ri, o.ro);
    const dr = temp.dr, dphi = temp.dphi;

    // ---- static tables, all computed in f64 and uploaded as f32 -------------
    const params = this.params;
    // The two element counts are what the mesh overlay divides the annulus into.
    // Taken from the axes rather than derived in the shader — a clamped cubic
    // axis has `n − 3` spans and a periodic one has `n`, and neither relation is
    // something a fragment shader should be asserting about `spline.ts`.
    new Int32Array(params, 0, 16).set([
      o.nr, o.na, o.gnr, o.gna, o.nr - 2, o.gnr - 2, rt.x.length, at.x.length,
      0, rAx.nLast, rAx.U.length, aAx.nLast,
      o.variable ? 1 : 0, rAx.elements().length, aAx.elements().length, 0,
    ]);
    this.pf = new Float32Array(params, 64, 16);
    this.pf.set([
      o.ri, o.ro, dr, dphi, temp.tIn, temp.tOut, o.Ra, o.dt,
      aAx.U[P], aAx.U[aAx.nLast + 1] - aAx.U[P], o.fill, o.levels, o.lineW,
      o.variable ? o.gamma : 0, o.variable ? o.n : 1, o.mesh,
    ]);

    const tri = diffusionFactors(o.gnr, o.gna, o.ri, dr, o.dt);
    const gni = o.gnr - 2, triFlat = new Float32Array(tri.length * 3 * gni);
    tri.forEach((f, k) => {
      triFlat.set(f.a, (k * 3 + 0) * gni);
      triFlat.set(f.cp, (k * 3 + 1) * gni);
      triFlat.set(f.m, (k * 3 + 2) * gni);
    });

    const T0 = new Float32Array(o.gnr * o.gna);
    for (let i = 0; i < o.gnr; i++) T0.set(temp.T[i], i * o.gna);

    const knots = new Float32Array(rAx.U.length + aAx.U.length);
    knots.set(rAx.U); knots.set(aAx.U, rAx.U.length);

    // ---- buffers ------------------------------------------------------------
    // COPY_SRC throughout so `read` can lift any buffer for verification.
    const ST = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC;
    const CD = GPUBufferUsage.COPY_DST;
    const upload = (name: string, data: ArrayBufferView, usage = ST) => {
      const b = device.createBuffer({ size: data.byteLength, usage: usage | CD });
      device.queue.writeBuffer(b, 0, data);
      this.buf[name] = b;
    };
    const scratch = (name: string, floats: number, usage = ST) => {
      this.buf[name] = device.createBuffer({ size: floats * S, usage: usage | CD });
    };

    /** Three per-DOF tables sharing one index array, interleaved for one binding. */
    const inter3 = (a: Float64Array, b: Float64Array, c: Float64Array) => {
      const out = new Float32Array(a.length * 3);
      for (let i = 0; i < a.length; i++) {
        out[3 * i] = a[i]; out[3 * i + 1] = b[i]; out[3 * i + 2] = c[i];
      }
      return out;
    };

    // COPY_SRC for the same reason as the storage buffers: the uniform block is
    // as much a thing the tests must be able to check as the fields are.
    upload("params", new Uint8Array(params),
      GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_SRC);
    upload("knots", knots);
    upload("rq", new Float32Array(rt.x));
    upload("phiq", new Float32Array(at.x));
    upload("rIdx", rt.idx);
    upload("rVal", new Float32Array(rt.val));
    upload("aIdx", at.idx);
    upload("aVal", new Float32Array(at.val));
    upload("inv", o.variable
      ? modeInverses(rAx, aAx, meanViscosity(o.ri, o.ro, o.gamma), true)
      : modeInverses(rAx, aAx));
    upload("tri", triFlat);
    upload("T", T0);

    const nq = rt.x.length * at.x.length;
    scratch("Tq", nq);
    scratch("G", rt.x.length * o.na);
    for (const n of ["b", "bRe", "bIm", "pRe", "pIm", "psi"]) scratch(n, o.nr * o.na);
    for (const n of ["TA", "TB", "tRe", "tIm", "dRe", "dIm"]) scratch(n, o.gnr * o.gna);
    scratch("stat", 4);   // [Nu inner, Nu outer, max|ψ|, —]
    this.buf.statRead = device.createBuffer({
      size: 4 * S, usage: GPUBufferUsage.MAP_READ | CD,
    });

    if (ot) {
      upload("rOp", inter3(ot.rA, ot.rG, ot.rB));
      upload("aOp", inter3(ot.aN0, ot.aN1, ot.aN2));
      // `mu` holds ε̇_II first and μ after, transformed in place: the law is
      // pointwise apart from ε̇_ref, and a second array of this size would be the
      // simulation's third largest for no reason.
      for (const n of ["Pq", "Qq", "mu"]) scratch(n, nq);
      scratch("Gop", 3 * rt.x.length * o.na);
      for (const n of ["res", "zv", "dir", "Aq"]) scratch(n, o.nr * o.na);
      scratch("sc", 8);   // [⟨r,z⟩, ⟨p,Ap⟩, ⟨r,z⟩′, —, δ, G, —, —]
    }

    // The three semi-Lagrangian passes differ only in (dt, limiter, bc).
    const pass = (name: string, dt: number, limiter: number, bc: number) => {
      this.buf[name] = device.createBuffer({ size: 16, usage: GPUBufferUsage.UNIFORM | CD });
      this.writePass(name, dt, limiter, bc);
    };
    pass("passA", o.dt, 0, 1);
    pass("passB", -o.dt, 0, 0);
    pass("passD", o.dt, 1, 1);

    // ---- pipelines and bind groups ------------------------------------------
    const group = (name: string, pipeline: GPUComputePipeline, res: string[]) => {
      this.pipe[name] = pipeline;
      this.bind[name] = device.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: res.map((r, binding) => ({ binding, resource: { buffer: this.buf[r] } })),
      });
    };
    const kernel = (name: string, code: string, ...res: string[]) => {
      group(name, device.createComputePipeline({
        layout: "auto",
        compute: { module: device.createShaderModule({ code }), entryPoint: "main" },
      }), res);
    };
    /** A second bind group over an existing pipeline — same kernel, other buffers. */
    const alias = (name: string, of: string, ...res: string[]) =>
      group(name, this.pipe[of], res);

    kernel("tq", w.tqSource(), "params", "T", "rq", "phiq", "Tq");
    kernel("g", w.gSource(SLOTS), "params", "Tq", "aIdx", "aVal", "G");
    kernel("b", w.bSource(SLOTS), "params", "G", "rIdx", "rVal", "b");
    kernel("fftA", w.fftForwardSource(o.na), "b", "bRe", "bIm");
    kernel("radial", w.radialSource(), "params", "bRe", "bIm", "inv", "pRe", "pIm");
    kernel("ifftA", w.fftInverseSource(o.na), "pRe", "pIm", "psi");

    const adv = w.advectSource();
    kernel("advA", adv, "params", "passA", "knots", "psi", "T", "T", "TA");
    kernel("advB", adv, "params", "passB", "knots", "psi", "TA", "TA", "TB");
    kernel("advD", adv, "params", "passD", "knots", "psi", "TA", "T", "TB");
    kernel("bfecc", w.bfeccSource(), "params", "T", "TB", "TA");

    kernel("fftG", w.fftForwardSource(o.gna), "TB", "tRe", "tIm");
    kernel("tridiag", w.tridiagSource(), "params", "tri", "tRe", "tIm", "dRe", "dIm");
    kernel("ifftG", w.fftInverseSource(o.gna), "dRe", "dIm", "T");
    kernel("nusselt", w.nusseltSource(), "params", "T", "stat");
    kernel("psiMax", w.psiMaxSource(), "params", "psi", "stat");

    if (!o.variable) return;

    // ---- the rheology, the matrix-free operator and the CG iteration --------
    //
    // The preconditioner is the tier-1 kernel with other buffers bound: the FFT
    // and the per-mode matvec are literally the same pipelines — "one kernel,
    // two jobs" made concrete rather than described.
    kernel("strain", w.strainSource(), "params", "knots", "psi", "rq", "phiq", "mu");
    kernel("sref", w.srefSource(), "params", "mu", "sc");
    kernel("muEval", w.muSource(), "params", "Tq", "sc", "mu");

    kernel("qevalPsi", w.qevalSource(),
      "params", "knots", "psi", "rq", "phiq", "mu", "Pq", "Qq");
    alias("qevalDir", "qevalPsi",
      "params", "knots", "dir", "rq", "phiq", "mu", "Pq", "Qq");
    kernel("gphiOp", w.gphiOpSource(SLOTS), "params", "Pq", "Qq", "aIdx", "aOp", "Gop");
    kernel("grOp", w.grOpSource(SLOTS), "params", "Gop", "rIdx", "rOp", "Aq");

    alias("preF", "fftA", "res", "bRe", "bIm");
    alias("preI", "ifftA", "pRe", "pIm", "zv");

    kernel("cgInit", w.cgInitSource(), "params", "b", "Aq", "res");
    kernel("cgCopy", w.cgCopySource(), "params", "zv", "dir");
    kernel("dotRZ0", w.dotSource(0), "params", "res", "zv", "sc");
    kernel("dotPAP", w.dotSource(1), "params", "dir", "Aq", "sc");
    kernel("dotRZ1", w.dotSource(2), "params", "res", "zv", "sc");
    kernel("cgUpdateX", w.cgUpdateXSource(),
      "params", "sc", "dir", "Aq", "psi", "res");
    kernel("cgUpdateP", w.cgUpdatePSource(), "params", "sc", "zv", "dir");
    kernel("cgRoll", w.cgRollSource(), "sc");
  }

  /**
   * Builds the tables (an O(n³) f64 inverse per mode) and compiles every
   * pipeline, so that no shader is ever created on a button press.
   */
  static create(device: GPUDevice, format: GPUTextureFormat, o: GpuSimOptions = {}) {
    const sim = new GpuSimulation(device, {
      nr: 32, na: 64, gnr: 65, gna: 128, ri: 0.55, ro: 1.0,
      Ra: 1e4, dt: 1e-3, fill: 0.92, levels: 0, lineW: 1.1, mesh: 0,
      variable: false, gamma: 0, iters: 12, n: 1, picard: 1, ...o,
    });
    sim.buildRender(format);
    sim.encode((enc) => sim.stokes(enc)); // ψ balancing the initial T
    return sim;
  }

  private buildRender(format: GPUTextureFormat): void {
    const module = this.device.createShaderModule({ code: w.renderSource() });
    this.render = this.device.createRenderPipeline({
      layout: "auto",
      vertex: { module, entryPoint: "vs" },
      fragment: { module, entryPoint: "fs", targets: [{ format }] },
    });
    this.renderBind = this.device.createBindGroup({
      layout: this.render.getBindGroupLayout(0),
      entries: ["params", "T", "knots", "psi", "stat"].map((r, binding) =>
        ({ binding, resource: { buffer: this.buf[r] } })),
    });
  }

  // ---- runtime controls ------------------------------------------------------
  //
  // Ra, the contour density, the line width, the mesh overlay and the power-law
  // index are pure uniform writes: nothing downstream of them is precomputed, so
  // they cost one 128-byte upload and take effect on the
  // next frame. `iters` and `picard` are
  // free in a different way — they are host-side loop bounds, so they change the
  // dispatch count and nothing else. dt and γ are the two knobs with an f64
  // factorisation behind them.

  private writePass(name: string, dt: number, limiter: number, bc: number): void {
    const d = new ArrayBuffer(16);
    new Float32Array(d, 0, 1)[0] = dt;
    new Int32Array(d, 4, 2).set([limiter, bc]);
    this.device.queue.writeBuffer(this.buf[name], 0, d);
  }

  private syncParams(): void {
    this.device.queue.writeBuffer(this.buf.params, 0, this.params);
  }

  get Ra(): number { return this.pf[F.Ra]; }
  set Ra(v: number) { this.pf[F.Ra] = v; this.syncParams(); }

  get dt(): number { return this.pf[F.dt]; }

  /**
   * Changing dt re-factorises (I − dt∇²) — 60k f64 flops on the CPU, so it is a
   * knob rather than a rebuild, but it is not free the way Ra is.
   */
  setDt(dt: number): void {
    const { gnr, gna, ri } = this.o, dr = (this.o.ro - ri) / (gnr - 1);
    const tri = diffusionFactors(gnr, gna, ri, dr, dt);
    const gni = gnr - 2, flat = new Float32Array(tri.length * 3 * gni);
    tri.forEach((f, k) => {
      flat.set(f.a, (k * 3 + 0) * gni);
      flat.set(f.cp, (k * 3 + 1) * gni);
      flat.set(f.m, (k * 3 + 2) * gni);
    });
    this.device.queue.writeBuffer(this.buf.tri, 0, flat);
    this.writePass("passA", dt, 0, 1);
    this.writePass("passB", -dt, 0, 0);
    this.writePass("passD", dt, 1, 1);
    this.pf[F.dt] = dt;
    this.syncParams();
  }

  get gamma(): number { return this.pf[F.gamma]; }

  /**
   * Set the viscosity contrast (γ = ln contrast) in the variable-μ tier.
   *
   * γ enters in two places and they cost very different things. In the operator
   * it is a uniform, read per quadrature point — free. In the **preconditioner**
   * it changes μ̄(r), and that means re-inverting one dense radial block per
   * azimuthal mode in f64: the same O(n³) job that dominates start-up. So this is
   * an announced rebuild, not a slider that tracks the pointer, and it is the
   * reason `rheology.ts` keeps μ̄ on a fixed profile rather than the running mean.
   */
  setGamma(gamma: number): void {
    if (!this.o.variable) throw new Error("γ has no effect in the constant-μ tier");
    this.device.queue.writeBuffer(this.buf.inv, 0,
      modeInverses(this.rAx, this.aAx,
        meanViscosity(this.o.ri, this.o.ro, gamma), true));
    this.pf[F.gamma] = gamma;
    this.syncParams();
  }

  get n(): number { return this.pf[F.n]; }

  /**
   * The power-law index. Unlike γ this is a *pure* uniform: n appears only
   * inside the pointwise law, never in the preconditioner — μ̄(r) has no
   * strain-rate profile to average, and the clamp keeps μ inside the same
   * interval the thermal term spans either way (see `rheology.ts`). So switching
   * between μ(T) and μ(T, ε̇) costs one 128-byte write, and only the *tier*
   * (which buffers exist at all) is a rebuild.
   */
  set n(v: number) { this.pf[F.n] = v; this.syncParams(); }

  /** `levels` contours across [−ψmax, ψmax]; 0 turns the streamlines off. */
  setStreamlines(levels: number, lineW = this.pf[F.lineW]): void {
    this.pf[F.levels] = levels;
    this.pf[F.lineW] = lineW;
    this.syncParams();
  }

  get mesh(): number { return this.pf[F.mesh]; }

  /**
   * Mesh overlay: 0 off, 1 the ψ spline elements, 2 the T grid cells. Both are
   * drawn from the element counts already in the uniform, so this is one 128-byte
   * write and no geometry — the mesh is a distance field in the fragment shader,
   * the same as the contours.
   */
  set mesh(mode: number) { this.pf[F.mesh] = mode; this.syncParams(); }

  /** Restart from the settled initial condition, clock included. */
  reseed(amp = 0.05, wavenumber = 4): void {
    const t = new Temperature(this.o.gnr, this.o.gna, this.o.ri, this.o.ro);
    t.reset(amp, wavenumber);
    this.writeTemperature(t.T);
    this.time = 0;
    this.steps = 0;
  }

  /** Release every GPU allocation. Required before building a replacement. */
  destroy(): void {
    for (const b of Object.values(this.buf)) b.destroy();
  }

  private dispatch(pass: GPUComputePassEncoder, name: string, n: number): void {
    pass.setPipeline(this.pipe[name]);
    pass.setBindGroup(0, this.bind[name]);
    pass.dispatchWorkgroups(Math.ceil(n / w.WG));
  }

  /** Rows-of-length-N transforms run one workgroup per row. */
  private rows(pass: GPUComputePassEncoder, name: string, rows: number): void {
    pass.setPipeline(this.pipe[name]);
    pass.setBindGroup(0, this.bind[name]);
    pass.dispatchWorkgroups(rows);
  }

  private encode(body: (enc: GPUCommandEncoder) => void): void {
    const enc = this.device.createCommandEncoder();
    body(enc);
    this.device.queue.submit([enc.finish()]);
  }

  /** Buoyancy load → Stokes solve → ψ. */
  private stokes(enc: GPUCommandEncoder): void {
    const p = enc.beginComputePass();
    this.dispatch(p, "tq", this.nrq * this.naq);
    this.dispatch(p, "g", this.nrq * this.na);
    this.dispatch(p, "b", this.nr * this.na);
    if (this.o.variable) {
      for (let s = 0; s < this.picard; s++) {
        this.rheology(p);
        this.krylov(p);
      }
    } else {
      this.rows(p, "fftA", this.nr);
      this.dispatch(p, "radial", this.nr * this.na);
      this.rows(p, "ifftA", this.nr);
    }
    this.dispatch(p, "psiMax", 1);   // contour scale for the streamline overlay
    p.end();
  }

  /**
   * μ at every quadrature point, from the *current* ψ.
   *
   * Three dispatches — strain rate, its RMS, the law — and then it is a fixed
   * coefficient field for the whole Krylov solve that follows. That lag is what
   * keeps the operator linear and symmetric under a law that is neither; at
   * n = 1 it is not a lag at all, since μ does not depend on ψ.
   *
   * Running it again (`picard` > 1) re-lags against the ψ the solve just wrote.
   * Dispatches within a compute pass are ordered and their writes visible to the
   * next, so the second sweep genuinely sees the updated field.
   */
  private rheology(p: GPUComputePassEncoder): void {
    const nq = this.nrq * this.naq;
    this.dispatch(p, "strain", nq);
    this.dispatch(p, "sref", 1);
    this.dispatch(p, "muEval", nq);
  }

  /** A·(psi | dir) → Aq, the three-pass matrix-free apply. */
  private applyOp(p: GPUComputePassEncoder, src: "qevalPsi" | "qevalDir"): void {
    this.dispatch(p, src, this.nrq * this.naq);
    this.dispatch(p, "gphiOp", this.nrq * this.na);
    this.dispatch(p, "grOp", this.nr * this.na);
  }

  /** z = M⁻¹ res — the tier-1 direct solve, other buffers bound. */
  private precondition(p: GPUComputePassEncoder): void {
    this.rows(p, "preF", this.nr);
    this.dispatch(p, "radial", this.nr * this.na);
    this.rows(p, "preI", this.nr);
  }

  /**
   * Preconditioned conjugate gradients, `iters` iterations, entirely on the GPU
   * (tier 2). The twin of `VariableStokes.solve`.
   *
   * α and β live in a four-float storage buffer and are consumed by the update
   * kernels directly, so no scalar ever crosses to the host — a convergence test
   * would need exactly that readback, which is why the budget is fixed instead.
   * Dispatches within one compute pass are ordered and their storage writes
   * visible to the next, so the data dependencies below need no explicit barrier.
   *
   * ψ is *not* cleared first: it is the warm start, and after the first frame it
   * is already the answer to within the O(dt) change in T.
   *
   * Nothing here reports a residual, and `pollStats` explains why none is
   * reported anywhere: the recursive `r` this loop carries drifts to a
   * meaningless 1e-17 in f32, and even a properly recomputed `b − Aψ` measures
   * ψ's f32 storage rather than the solve.
   */
  private krylov(p: GPUComputePassEncoder): void {
    const n = this.nr * this.na;
    this.applyOp(p, "qevalPsi");        // A x₀
    this.dispatch(p, "cgInit", n);      // res = b − A x₀
    this.precondition(p);               // z = M⁻¹ res
    this.dispatch(p, "cgCopy", n);      // dir = z
    this.dispatch(p, "dotRZ0", 1);      // ⟨r,z⟩
    for (let i = 0; i < this.iters; i++) {
      this.applyOp(p, "qevalDir");      // A p
      this.dispatch(p, "dotPAP", 1);
      this.dispatch(p, "cgUpdateX", n); // x += α p, res −= α Ap
      this.precondition(p);
      this.dispatch(p, "dotRZ1", 1);
      this.dispatch(p, "cgUpdateP", n); // p = z + β p
      this.dispatch(p, "cgRoll", 1);    // ⟨r,z⟩ ← ⟨r,z⟩′
    }
  }

  private get nrq(): number { return this.buf.rq.size / S; }
  private get naq(): number { return this.buf.phiq.size / S; }

  /** Encode and submit one full time step. */
  step(): void {
    const grid = this.gnr * this.gna;
    this.encode((enc) => {
      const p = enc.beginComputePass();
      this.dispatch(p, "advA", grid);
      this.dispatch(p, "advB", grid);
      this.dispatch(p, "bfecc", grid);
      this.dispatch(p, "advD", grid);
      this.rows(p, "fftG", this.gnr);
      this.dispatch(p, "tridiag", this.gna * 2);
      this.rows(p, "ifftG", this.gnr);
      this.dispatch(p, "nusselt", 1);
      p.end();
      // psiMax runs inside `stokes`, right after ψ is written.
      this.stokes(enc);
    });
    this.time += this.dt;
    this.steps++;
  }

  /** Draw the current temperature field straight from its storage buffer. */
  draw(view: GPUTextureView): void {
    this.encode((enc) => {
      const p = enc.beginRenderPass({
        colorAttachments: [{
          view,
          clearValue: { r: 0.02, g: 0.02, b: 0.05, a: 1 },
          loadOp: "clear", storeOp: "store",
        }],
      });
      p.setPipeline(this.render);
      p.setBindGroup(0, this.renderBind);
      p.draw(3);
      p.end();
    });
  }

  /**
   * Refresh `stats` from the GPU-side reductions. Off the frame's dependency
   * chain: the copy is queued behind work already submitted and the map resolves
   * whenever it resolves, so the hot loop never blocks on it. `max|ψ|` is read
   * only for the HUD — the renderer takes it from the buffer directly.
   *
   * **There is deliberately no residual here.** It was tried, both as CG's
   * recursive `r` and as a properly recomputed `b − Aψ`, and neither is a
   * convergence diagnostic in this solver. The recursive one drifts to a
   * meaningless 1e-17 in f32. The recomputed one is worse than useless: ψ is
   * *stored* in f32, so perturbing it by ε_f32 already moves the residual by
   * κ‖b‖, and with κ ~ h⁻⁴ that floor is ~3e-2 at ψ 96×256 — measured flat in
   * the iteration count from 4 to 20, because it is reporting the width of an
   * f32 mantissa rather than anything the solver did. This is the f64-at-init
   * argument seen from the other side: the *solution* is accurate because the
   * matvec is backward stable, and precisely that decoupling makes the residual a bad
   * proxy for it. Convergence is asserted in f64, in `tests/rheology.test.ts`,
   * where the quantity means what it says.
   */
  pollStats(): void {
    if (this.statPending) return;
    this.statPending = true;
    this.encode((enc) =>
      enc.copyBufferToBuffer(this.buf.stat, 0, this.buf.statRead, 0, 4 * S));
    void this.buf.statRead.mapAsync(GPUMapMode.READ).then(() => {
      const a = new Float32Array(this.buf.statRead.getMappedRange().slice(0));
      this.buf.statRead.unmap();
      this.stats = { nuInner: a[0], nuOuter: a[1], psiMax: a[2] };
    }).catch(() => {
      // A failed map is a diagnostic, not a simulation fault — drop it and let
      // the next poll try again rather than wedging the flag.
    }).finally(() => { this.statPending = false; });
  }

  /** Copy a solver buffer to the host. Verification only — never in the loop. */
  async read(name: string): Promise<Float32Array> {
    const src = this.buf[name];
    const dst = this.device.createBuffer({
      size: src.size, usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST,
    });
    this.encode((enc) => enc.copyBufferToBuffer(src, 0, dst, 0, src.size));
    await dst.mapAsync(GPUMapMode.READ);
    const out = new Float32Array(dst.getMappedRange().slice(0));
    dst.unmap();
    dst.destroy();
    return out;
  }

  /**
   * Overwrite the temperature grid and re-solve — used to align initial states
   * in tests, and by `reseed`.
   *
   * ψ is cleared first. The Krylov tier warm-starts from it, and a ψ balancing
   * the *old* field is not a guess about the new one, it is a wrong answer the
   * fixed budget would only partly walk back. Zeroing makes the solve depend on
   * nothing but the T just written, which is also what makes it comparable to a
   * freshly constructed CPU reference. Tier 1 is unaffected — the direct solve
   * never reads ψ.
   */
  writeTemperature(T: Float64Array[]): void {
    const f = new Float32Array(this.gnr * this.gna);
    for (let i = 0; i < this.gnr; i++) f.set(T[i], i * this.gna);
    this.device.queue.writeBuffer(this.buf.T, 0, f);
    if (this.o.variable)
      this.device.queue.writeBuffer(this.buf.psi, 0,
        new Float32Array(this.nr * this.na));
    this.encode((enc) => this.stokes(enc));
  }
}
