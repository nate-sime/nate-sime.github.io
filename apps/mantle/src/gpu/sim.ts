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
 * **Adaptive dt, host-lagged.** The CPU reference sizes dt from the advective
 * CFL by reducing max|u| on the host every step — a readback the GPU frame loop
 * cannot afford. Semi-Lagrangian advection and implicit diffusion are both
 * unconditionally stable, so dt is an accuracy knob, not a stability limit: the
 * same max-speed measure is reduced on the GPU instead (`cflSource`, into
 * `stat[4]`) and carried back through the existing asynchronous `pollStats`.
 * `main.ts` turns it into a step through `adaptiveDt` and calls `setDt` —
 * which re-factorises the Thomas tables in f64 — only when that step has
 * drifted past a hysteresis band, so the factorisation stays occasional rather
 * than becoming a per-step cost.
 */

import { type ColormapName } from "../colormaps";
import { ANNULUS, radialMargin, type Geometry } from "../geometry";
import { Axis, P, clampedAxis, periodicAxis } from "../spline";
import { modeInverses } from "../solver/operators";
import { operatorTables, quadTable, SLOTS, W0, W1 } from "../solver/assembly";
import {
  meanViscosity, meanTackleyViscosity, meanBlankenbachViscosity,
  meanVanKekenViscosity,
} from "../solver/rheology";
import { Temperature, diffusionFactors } from "../solver/temperature";
import { DEFAULT_LAYER_DEPTH } from "../particles";
import { GpuParticles } from "./particles";
import * as w from "./wgsl";

export interface GpuSimOptions {
  nr?: number; na?: number;      // spline space for ψ
  gnr?: number; gna?: number;    // grid for T
  /**
   * Domain and metric — the annulus, or a Cartesian box of some length. A
   * construction-time choice for the same reason `variable` is: the metric is
   * emitted into the WGSL rather than branched on a uniform (see `wgsl.ts`), and
   * the box length reaches the knot vector, so the quadrature tables and the
   * per-mode inverses are both built against it.
   */
  geom?: Geometry;
  Ra?: number; dt?: number;
  fill?: number;                 // fraction of the half-viewport the domain spans
  levels?: number;               // streamline contours across [−ψmax, ψmax]; 0 = off
  lineW?: number;                // contour half-width, in pixels
  mesh?: number;                 // mesh overlay: 0 = off, 1 = ψ elements, 2 = T grid
  /**
   * Temperature colour map. A render-pipeline-only rebuild (`setColormap`),
   * not a uniform: the control points are compiled into the fragment shader
   * rather than read from a buffer, so a run of five `mix`es never touches an
   * indirection the GPU would otherwise pay for every pixel.
   */
  colormap?: ColormapName;
  /**
   * Tier 2: μ(T, d) by matrix-free PCG instead of the direct DFT solve. A
   * construction-time choice, not a uniform — the Krylov path allocates the
   * quadrature-point stress buffers, which are the largest in the simulation
   * after the inverse table, so a constant-μ run must not pay for them.
   */
  variable?: boolean;
  gamma?: number;                // ln of the viscosity contrast across T
  /**
   * ln of the viscosity contrast across the depth of the layer. Like γ it is a
   * uniform *and* a preconditioner rebuild — μ̄(r) carries the depth term
   * exactly — so it is set through `setContrast`, not written on its own.
   */
  cz?: number;
  iters?: number;                // Krylov budget per solve
  /**
   * Power-law index. `n = 1` is μ(T) exactly, so the two variable laws are
   * one tier and one set of pipelines — this is a uniform, not a rebuild.
   */
  n?: number;
  picard?: number;               // rheology updates per solve; 1 = pure time-lagging
  /**
   * Tackley pseudo-plastic law instead of the power law. A construction-time
   * choice like `variable`: it selects a different pointwise-law kernel and
   * skips the `sref` reduction the power law needs (Tackley reads ε̇ raw).
   */
  tackley?: boolean;
  /**
   * Tosi et al. (2015) viscoplastic law instead of the power law. Also a
   * construction-time choice, for the same reason `tackley` is — its own
   * pointwise kernel, and no `sref` pass either. Unlike Tackley, `gamma`/`cz`
   * are load-bearing here (they are the paper's own γ_T/γ_z), so it keeps
   * the contrast sliders and their preconditioner path — see `setContrast`.
   */
  tosi?: boolean;
  /**
   * Blankenbach et al. (1989)'s own μ(T, d) = exp(−bT + cd) instead of the
   * power law. A construction-time choice like `tackley` and `tosi`: its own
   * pointwise kernel, no `sref` pass (no strain-rate dependence at all, so
   * there is nothing to normalise). Reads `gamma`/`cz` as the paper's own b,
   * c — the same two uniforms μ(T, d) and Tosi read — so the app's contrast
   * sliders and Ra enter this case exactly as the paper states them, with no
   * rescaling for the reader to work out first (see `rheology.ts`).
   */
  blankenbach?: boolean;
  /** Constant ductile yield stress. Pure uniform: μ̄(r) does not depend on it. Shared by Tackley and Tosi. */
  sigmaY?: number;
  /** Gradient of brittle yield stress with depth. Pure uniform. Shared by Tackley and Tosi. */
  sigmaB?: number;
  /** Minimum plastic viscosity. Pure uniform. Shared by Tackley and Tosi. */
  etaStar?: number;
  /**
   * van Keken et al. (1997)'s composition-linear law, μ(φ) = η_light +
   * φ(η_dense − η_light), instead of the power law. A construction-time
   * choice like `tackley`/`tosi`/`blankenbach`: its own pointwise kernel, no
   * `sref` pass. Unlike every other law here it has no T dependence at all —
   * it reads the tracer cloud's composition, not `Tq` — which is what makes
   * it the one law that combines safely with a nonzero `Rb` (see `Rb`'s own
   * header in `wgsl.ts`).
   */
  vanKeken?: boolean;
  /** Viscosity of the light material, van Keken's own law. Pure uniform: μ̄(r) uses it only as a reference profile at construction. */
  etaLight?: number;
  /** Viscosity of the dense material, van Keken's own law. Pure uniform, same caveat as `etaLight`. */
  etaDense?: number;
  /**
   * Reference interface height μ̄(r) is built from for the van Keken law — the
   * *unperturbed* version of whatever the tracer cloud's own `layerDepth` is,
   * since the preconditioner wants a φ-independent radial profile and the
   * actual composition is not one. Duplicated from the particle system's own
   * option rather than read from it, because a preconditioner profile is
   * needed at construction, before any `GpuParticles` may exist to read it
   * from — see `main.ts`'s `attachParticles`, which is what keeps the two
   * numbers in sync in practice (both come from `state.layerDepth`).
   */
  layerDepth?: number;
}

const S = Float32Array.BYTES_PER_ELEMENT;

/**
 * Indices into the uniform block (see `PARAMS` in wgsl.ts). Only these change at
 * runtime; everything else is fixed by the resolution.
 */
const F = {
  Ra: 6, dt: 7, levels: 11, lineW: 12, gamma: 13, n: 14, mesh: 15, cz: 16,
  zoom: 17, panX: 18, panY: 19, sigmaY: 20, sigmaB: 21, etaStar: 22, Rb: 24,
  etaLight: 25, etaDense: 26,
} as const;

export class GpuSimulation {
  readonly nr: number; readonly na: number;
  readonly gnr: number; readonly gna: number;
  time = 0;
  steps = 0;
  /** Krylov budget per solve; free to change (it only alters the dispatch count). */
  iters: number;
  /** Rheology updates per solve; likewise only a dispatch count. */
  picard: number;
  /**
   * Latest diagnostics; refreshed by `pollStats`, stale between polls.
   *
   * `at` and `atStep` are the simulation time and step count they describe, not
   * the ones current when they arrived. The readout does not need either — a
   * number on screen is understood to be the current one — but the Nusselt and
   * v_rms plots are time series: `at` is a sample's abscissa and `atStep` is
   * what its display window is measured in. Recording them when the copy is
   * *encoded* rather than when the map resolves is what makes them right: the reduction in
   * `stat` was written by the last step submitted before that point, and `step`
   * advances both counters after submitting, so they are already that step's end.
   */
  stats = {
    nuInner: NaN, nuOuter: NaN, psiMax: NaN, vrms: NaN, maxSpeed: NaN,
    at: NaN, atStep: NaN,
  };

  /**
   * The tracer overlay, when one has been attached (`sim.particles = new
   * GpuParticles(sim, opts)`) — null for a plain thermal run. A public,
   * ordinary field rather than something this class manages the lifetime of:
   * `GpuParticles` owns its own buffers and pipelines and hooks itself into
   * the buoyancy load on construction and back out on `destroy` (see that
   * file's header), so this class only ever needs to know whether one is
   * currently present at the four points a step, a draw or a re-seed touches it.
   */
  particles: GpuParticles | null = null;

  private readonly rAx: Axis;
  private readonly aAx: Axis;
  private readonly buf: Record<string, GPUBuffer> = {};
  private readonly pipe: Record<string, GPUComputePipeline> = {};
  private readonly bind: Record<string, GPUBindGroup> = {};
  private render!: GPURenderPipeline;
  private renderBind!: GPUBindGroup;
  private canvasFormat!: GPUTextureFormat;
  private statPending = false;
  // 176, not 160: `PARAMS` (wgsl.ts) carries one more f32 row than this class
  // otherwise uses, reserved for `Rb` — the particle feature's coupling
  // switch — and `etaLight`/`etaDense`, the van Keken tier's own two
  // parameters. Sized to match now so the uniform buffer is never smaller
  // than the struct the shaders declare; `Rb` itself is left at its
  // zero-initialised default (uncoupled) until `GpuParticles` writes it.
  private readonly params = new ArrayBuffer(176);
  private readonly pf: Float32Array;

  private constructor(readonly device: GPUDevice, readonly o: Required<GpuSimOptions>) {
    this.nr = o.nr; this.na = o.na; this.gnr = o.gnr; this.gna = o.gna;
    this.iters = o.iters;
    this.picard = o.picard;
    for (const n of [o.na, o.gna])
      if (!Number.isInteger(Math.log2(n))) throw new Error(`${n} is not a power of two`);

    const g = o.geom;
    const rAx = clampedAxis(o.nr, g.lo, g.hi), aAx = periodicAxis(o.na, g.span);
    this.rAx = rAx; this.aAx = aAx;
    // The buoyancy load's two tables, `w·B` and `w·N′`. `ot` adds the six the
    // variable-μ operator gathers against, and is built only in that tier — the
    // constant path would otherwise construct twelve tables to upload none of
    // them. Its eval half is unused here either way: the GPU recomputes the
    // basis per quadrature point rather than tabulating it (see `qevalSource`).
    const rt = quadTable(rAx, W0), at = quadTable(aAx, W1);
    const ot = o.variable ? operatorTables(rAx, aAx, g) : null;
    const temp = new Temperature(g, o.gnr, o.gna);
    const dr = temp.dr, dphi = temp.dphi;

    // ---- static tables, all computed in f64 and uploaded as f32 -------------
    const params = this.params;
    // The two element counts are what the mesh overlay divides the domain into.
    // Taken from the axes rather than derived in the shader — a clamped cubic
    // axis has `n − 3` spans and a periodic one has `n`, and neither relation is
    // something a fragment shader should be asserting about `spline.ts`.
    new Int32Array(params, 0, 16).set([
      o.nr, o.na, o.gnr, o.gna, o.nr - 2 * radialMargin(g), o.gnr - 2, rt.x.length, at.x.length,
      0, rAx.nLast, rAx.U.length, aAx.nLast,
      o.variable ? 1 : 0, rAx.elements().length, aAx.elements().length, 0,
    ]);
    this.pf = new Float32Array(params, 64, 28);
    this.pf.set([
      g.lo, g.hi, dr, dphi, temp.tIn, temp.tOut, o.Ra, o.dt,
      aAx.U[P], aAx.U[aAx.nLast + 1] - aAx.U[P], o.fill, o.levels, o.lineW,
      o.variable ? o.gamma : 0, o.variable ? o.n : 1, o.mesh,
      o.variable ? o.cz : 0,
      1, 0, 0,   // zoom, panX, panY — the identity view; see `setView`
      (o.tackley || o.tosi) ? o.sigmaY : 0,
      (o.tackley || o.tosi) ? o.sigmaB : 0,
      (o.tackley || o.tosi) ? o.etaStar : 0,
      0,   // fpad0
      // Rb left at its zero-initialised default (uncoupled) until
      // `GpuParticles` writes it; etaLight/etaDense are read only by the van
      // Keken law, but a viscosity of 0 is not an inert default the way 0 is
      // for the others, so they get real values regardless of `vanKeken`.
      0, o.etaLight, o.etaDense, 0,
    ]);

    const tri = diffusionFactors(g, o.gnr, o.gna, dr, o.dt);
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
      ? modeInverses(rAx, aAx,
        o.tackley ? meanTackleyViscosity(g)
          : o.blankenbach ? meanBlankenbachViscosity(g, o.gamma, o.cz)
          : o.vanKeken ? meanVanKekenViscosity(g, o.etaLight, o.etaDense, o.layerDepth)
          : meanViscosity(g, o.gamma, o.cz), true, g)
      : modeInverses(rAx, aAx, () => 1, false, g));
    upload("tri", triFlat);
    upload("T", T0);

    const nq = rt.x.length * at.x.length;
    scratch("Tq", nq);
    scratch("Cq", nq);
    // Placeholder `Part` uniform and composition buffer, so the composition
    // gather (`cqSource`) and, in the van Keken tier, `muEval` always have
    // something satisfiable to bind even though nothing here turns coupling
    // on: `Rb` stays at its zero-initialised default, so this placeholder
    // `C` is read but its (all-zero) value never reaches the buoyancy load —
    // and the van Keken tier's `mu[g]` reduces to the isoviscous `etaLight`
    // wherever no cloud has projected anything. A future `GpuParticles` owns
    // and overwrites both of these; this is only what keeps a thermal-only
    // run buildable without it.
    const partDefault = new ArrayBuffer(48);
    new Int32Array(partDefault, 0, 4).set([0, 2, 2, 0]);
    upload("part", new Uint8Array(partDefault),
      GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_SRC);
    scratch("C", 4); // cnr · cna = 2 · 2, zero-initialised
    scratch("G", rt.x.length * o.na);
    scratch("GC", rt.x.length * o.na);
    for (const n of ["b", "bRe", "bIm", "pRe", "pIm", "psi"]) scratch(n, o.nr * o.na);
    for (const n of ["TA", "TB", "tRe", "tIm", "dRe", "dIm"]) scratch(n, o.gnr * o.gna);
    scratch("stat", 5);   // [Nu inner, Nu outer, max|ψ|, v_rms, max CFL speed]
    this.buf.statRead = device.createBuffer({
      size: 5 * S, usage: GPUBufferUsage.MAP_READ | CD,
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
    // The compositional half of the load: its own gather chain (see
    // `bcSource`'s own header on why this is not folded into `tq`/`b`
    // above), bound to the placeholder `part`/`C` at construction and
    // rebound to a live cloud's buffers by `bindBuoyancy` — the same
    // treatment `tq` itself used to need before this split.
    kernel("cq", w.cqSource(), "params", "rq", "phiq", "part", "C", "Cq");
    alias("gC", "g", "params", "Cq", "aIdx", "aVal", "GC");
    kernel("bC", w.bcSource(SLOTS), "params", "GC", "rIdx", "rVal", "b");
    kernel("fftA", w.fftForwardSource(o.na), "b", "bRe", "bIm");
    kernel("radial", w.radialSource(), "params", "bRe", "bIm", "inv", "pRe", "pIm");
    kernel("ifftA", w.fftInverseSource(o.na), "pRe", "pIm", "psi");

    const adv = w.advectSource(g);
    kernel("advA", adv, "params", "passA", "knots", "psi", "T", "T", "TA");
    kernel("advB", adv, "params", "passB", "knots", "psi", "TA", "TA", "TB");
    kernel("advD", adv, "params", "passD", "knots", "psi", "TA", "T", "TB");
    kernel("bfecc", w.bfeccSource(), "params", "T", "TB", "TA");

    kernel("fftG", w.fftForwardSource(o.gna), "TB", "tRe", "tIm");
    kernel("tridiag", w.tridiagSource(g), "params", "tri", "tRe", "tIm", "dRe", "dIm");
    kernel("ifftG", w.fftInverseSource(o.gna), "dRe", "dIm", "T");
    kernel("nusselt", w.nusseltSource(g), "params", "T", "stat");
    kernel("psiMax", w.psiMaxSource(), "params", "psi", "stat");
    kernel("rms", w.rmsSource(g), "params", "knots", "psi", "stat");
    kernel("cfl", w.cflSource(g), "params", "knots", "psi", "stat");

    if (!o.variable) return;

    // ---- the rheology, the matrix-free operator and the CG iteration --------
    //
    // The preconditioner is the tier-1 kernel with other buffers bound: the FFT
    // and the per-mode matvec are literally the same pipelines — "one kernel,
    // two jobs" made concrete rather than described.
    kernel("strain", w.strainSource(g), "params", "knots", "psi", "rq", "phiq", "mu");
    if (o.tackley) {
      // Tackley reads ε̇ raw — no normalising reduction, so no `sref` pass.
      kernel("muEval", w.tackleyMuSource(), "params", "Tq", "rq", "mu");
    } else if (o.tosi) {
      // Tosi reads ε̇ raw too, for the same reason — no `sref` pass.
      kernel("muEval", w.tosiMuSource(), "params", "Tq", "rq", "mu");
    } else if (o.blankenbach) {
      // No ε̇ at all, so no `sref` pass either — see `rheology()` below.
      kernel("muEval", w.blankenbachMuSource(), "params", "Tq", "rq", "mu");
    } else if (o.vanKeken) {
      // No ε̇ and no `Tq` either — reads the composition directly, bound to
      // the same placeholder `part`/`C` as `cq` until a cloud attaches (see
      // `bindBuoyancy`).
      kernel("muEval", w.vanKekenMuSource(), "params", "rq", "phiq", "part", "C", "mu");
    } else {
      kernel("sref", w.srefSource(), "params", "mu", "sc");
      kernel("muEval", w.muSource(), "params", "Tq", "sc", "rq", "mu");
    }

    kernel("qevalPsi", w.qevalSource(g),
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
      nr: 32, na: 64, gnr: 65, gna: 128, geom: ANNULUS,
      Ra: 1e4, dt: 1e-3, fill: 0.92, levels: 0, lineW: 1.1, mesh: 0,
      variable: false, gamma: 0, cz: 0, iters: 12, n: 1, picard: 1,
      tackley: false, sigmaY: 1, sigmaB: 1, etaStar: 1e-3,
      tosi: false, blankenbach: false,
      vanKeken: false, etaLight: 1, etaDense: 1, layerDepth: DEFAULT_LAYER_DEPTH,
      colormap: "inferno", ...o,
    });
    sim.buildRender(format);
    sim.encode((enc) => sim.stokes(enc)); // ψ balancing the initial T
    return sim;
  }

  /** The canvas format this simulation renders into — what a `GpuParticles` attached to it must target too. */
  get format(): GPUTextureFormat { return this.canvasFormat; }

  /**
   * Buffer lookup by name — the small, explicit surface `GpuParticles`
   * borrows through (`params`, `passA`, `knots`, `psi`, `T`, `stat`) rather
   * than a broader handle onto this class's internals.
   */
  buffer(name: string): GPUBuffer { return this.buf[name]; }

  /**
   * Repoint the composition-reading kernels at a live tracer cloud's
   * `part`/`C` buffers — called once when a `GpuParticles` attaches. Bind
   * groups are immutable once built, so this rebuilds them rather than
   * writing into them. Two kernels ever read composition: `cq`, the
   * compositional half of the buoyancy load, always; and `muEval`, only in
   * the van Keken tier, where the pointwise law reads composition instead of
   * temperature (`vanKekenMuSource`) — every other tier's `muEval` never
   * touches `part`/`C` at all, so rebinding it here would be rebuilding a
   * bind group nothing uses.
   */
  bindBuoyancy(part: GPUBuffer, C: GPUBuffer): void {
    this.bind.cq = this.device.createBindGroup({
      layout: this.pipe.cq.getBindGroupLayout(0),
      entries: [this.buf.params, this.buf.rq, this.buf.phiq, part, C, this.buf.Cq]
        .map((buffer, binding) => ({ binding, resource: { buffer } })),
    });
    if (this.o.vanKeken) {
      this.bind.muEval = this.device.createBindGroup({
        layout: this.pipe.muEval.getBindGroupLayout(0),
        entries: [this.buf.params, this.buf.rq, this.buf.phiq, part, C, this.buf.mu]
          .map((buffer, binding) => ({ binding, resource: { buffer } })),
      });
    }
  }

  /** Hand the composition reads back to the zeroed placeholders — called when a `GpuParticles` detaches. */
  unbindBuoyancy(): void { this.bindBuoyancy(this.buf.part, this.buf.C); }

  get Rb(): number { return this.pf[F.Rb]; }

  /**
   * The compositional Rayleigh number. A pure uniform write, like Ra — the
   * host-side check `step`/`stokes` makes before dispatching the
   * compositional gather chain at all (`this.particles && this.Rb !== 0`) is
   * what lets the thermochemical coupling be switched on and off without
   * touching a pipeline, the same property `pp.buoy != 0.0` used to give
   * `tqSource` before the compositional term moved to its own dispatches
   * (see `bcSource`). Left at its zero-initialised default until this is
   * called, which is also what keeps a plain thermal run — with or without a
   * tracer overlay attached — unperturbed by particles existing at all.
   *
   * Independent of `Ra`, on purpose: the buoyancy load reads `Ra·T − Rb·C`,
   * not `Ra·(T − B·C)`, so `Rb` can be nonzero while `Ra = 0` — the purely
   * compositional (isothermal) limit the Rayleigh–Taylor case needs, which a
   * single ratio `B = Rb/Ra` could never reach.
   */
  set Rb(v: number) { this.pf[F.Rb] = v; this.syncParams(); }

  get etaLight(): number { return this.pf[F.etaLight]; }
  get etaDense(): number { return this.pf[F.etaDense]; }

  /**
   * Set the van Keken law's two viscosities together — **not** a pure
   * uniform write, unlike `sigmaY`/`sigmaB`/`etaStar`: `etaLight`/`etaDense`
   * also set μ̄(r) (`meanVanKekenViscosity`), so changing either re-inverts
   * the preconditioner in f64, the same announced cost `setContrast` pays
   * for `gamma`/`cz`. Set together for the same reason those two are: one
   * f64 job over a profile that is a function of both.
   */
  setViscosity(etaLight = this.etaLight, etaDense = this.etaDense): void {
    if (!this.o.vanKeken) throw new Error("η_light/η_dense have no effect outside the van Keken law");
    this.device.queue.writeBuffer(this.buf.inv, 0,
      modeInverses(this.rAx, this.aAx,
        meanVanKekenViscosity(this.o.geom, etaLight, etaDense, this.o.layerDepth),
        true, this.o.geom));
    this.pf[F.etaLight] = etaLight;
    this.pf[F.etaDense] = etaDense;
    this.syncParams();
  }

  private buildRender(format: GPUTextureFormat): void {
    this.canvasFormat = format;
    const module = this.device.createShaderModule({
      code: w.renderSource(this.o.geom, this.o.colormap),
    });
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
  // dispatch count and nothing else. dt and the two contrasts are the knobs with
  // an f64 factorisation behind them.

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
   * knob rather than a rebuild, but it is not free the way Ra is. Called from
   * `main.ts`'s adaptive loop through `adaptiveDt`'s hysteresis band, not on
   * every step.
   */
  setDt(dt: number): void {
    const { gnr, gna, geom } = this.o, dr = (geom.hi - geom.lo) / (gnr - 1);
    const tri = diffusionFactors(geom, gnr, gna, dr, dt);
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

  get cz(): number { return this.pf[F.cz]; }

  /**
   * Set the two viscosity contrasts — γ = ln of the contrast across T, c = ln of
   * the contrast across depth — in the variable-μ tier.
   *
   * Both enter in two places and they cost very different things. In the
   * operator they are uniforms, read per quadrature point — free. In the
   * **preconditioner** they change μ̄(r), and that means re-inverting one dense
   * radial block per azimuthal mode in f64: the same O(n³) job that dominates
   * start-up. So this is an announced rebuild, not a slider that tracks the
   * pointer, and it is the reason `rheology.ts` keeps μ̄ on a fixed profile
   * rather than the running mean.
   *
   * They are set together because they share that one f64 job: writing them
   * separately would pay for it twice to reach a state either could have
   * described. The defaults are the current values, so moving one slider leaves
   * the other where it was.
   *
   * In the Tosi tier these two sliders are the paper's own γ_T and γ_z, not
   * this app's power-law γ and c — the same contrast, read by a different
   * formula, but the same `meanViscosity` profile: Tosi's linear branch *is*
   * `viscosity` at `strain = 1, n = 1` (see `rheology.ts`), so there is
   * nothing for this method to special-case.
   *
   * The Blankenbach tier *does* need a different profile — `meanViscosity`
   * is the centred law's ε̇ → 0 limit, and Blankenbach's own law is not
   * centred — so this is the one place `blankenbach` is read, rather than
   * being handled implicitly the way Tosi's reuse of `meanViscosity` is.
   */
  setContrast(gamma = this.gamma, cz = this.cz): void {
    if (!this.o.variable) throw new Error("γ has no effect in the constant-μ tier");
    const profile = this.o.blankenbach
      ? meanBlankenbachViscosity(this.o.geom, gamma, cz)
      : meanViscosity(this.o.geom, gamma, cz);
    this.device.queue.writeBuffer(this.buf.inv, 0,
      modeInverses(this.rAx, this.aAx, profile, true, this.o.geom));
    this.pf[F.gamma] = gamma;
    this.pf[F.cz] = cz;
    this.syncParams();
  }

  get n(): number { return this.pf[F.n]; }

  /**
   * The power-law index. Unlike the two contrasts this is a *pure* uniform: n
   * appears only inside the pointwise law, never in the preconditioner — μ̄(r)
   * has no strain-rate profile to average, and the clamp keeps μ inside the same
   * interval the thermal and depth terms span either way (see `rheology.ts`). So
   * switching between μ(T, d) and μ(T, d, ε̇) costs one uniform write, and only
   * the *tier* (which buffers exist at all) is a rebuild.
   */
  set n(v: number) { this.pf[F.n] = v; this.syncParams(); }

  get sigmaY(): number { return this.pf[F.sigmaY]; }

  /** Constant ductile yield stress, Tackley and Tosi laws. A pure uniform — see `GpuSimOptions.sigmaY`. */
  set sigmaY(v: number) { this.pf[F.sigmaY] = v; this.syncParams(); }

  get sigmaB(): number { return this.pf[F.sigmaB]; }

  /** Gradient of brittle yield stress with depth, Tackley and Tosi laws. A pure uniform. */
  set sigmaB(v: number) { this.pf[F.sigmaB] = v; this.syncParams(); }

  get etaStar(): number { return this.pf[F.etaStar]; }

  /** Minimum plastic viscosity, Tackley and Tosi laws. A pure uniform. */
  set etaStar(v: number) { this.pf[F.etaStar] = v; this.syncParams(); }

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

  /**
   * Half the world-space width of the fixed `[-halfExtent, halfExtent]²` window
   * that a zoom of 1 shows — the same quantity `halfExtent()` in the vertex
   * shader computes, kept here so the host can turn a screen-pixel offset into
   * the world offset `setView` wants without round-tripping through the GPU.
   */
  get halfExtent(): number {
    const { geom, fill } = this.o;
    return geom.kind === "annulus"
      ? geom.hi / fill
      : 0.5 * Math.max(geom.width, geom.hi - geom.lo) / fill;
  }

  get zoom(): number { return this.pf[F.zoom]; }
  get panX(): number { return this.pf[F.panX]; }
  get panY(): number { return this.pf[F.panY]; }

  /**
   * The view transform: `zoom` magnifies about the world origin, `(panX, panY)`
   * then recentres it — both in the same world units as `halfExtent`. A pure
   * uniform write, like Ra: the vertex shader reparametrises which world window
   * the fixed full-screen triangle covers, so nothing here is precomputed.
   */
  setView(zoom: number, panX: number, panY: number): void {
    this.pf[F.zoom] = zoom;
    this.pf[F.panX] = panX;
    this.pf[F.panY] = panY;
    this.syncParams();
  }

  get colormap(): ColormapName { return this.o.colormap; }

  /**
   * Swap the temperature colour map. The control points are compiled into
   * the fragment shader (see `wgsl.ts`), so this rebuilds only the render
   * pipeline — one shader module and one pipeline, sharing every buffer the
   * old one bound — rather than anything the solve depends on.
   */
  setColormap(name: ColormapName): void {
    if (name === this.o.colormap) return;
    this.o.colormap = name;
    this.buildRender(this.format);
  }

  /** Restart from the settled initial condition, clock included — and draw a fresh tracer cloud, if one is attached, since that is what pressing "reseed" means for it too. */
  reseed(amp = 0.05, wavenumber = 4): void {
    const t = new Temperature(this.o.geom, this.o.gnr, this.o.gna);
    t.reset(amp, wavenumber);
    this.particles?.seed();
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
    // The compositional half of the load, subtracted into `b` in place —
    // dispatched only when it could actually change anything, the same
    // guard `particles.project` below uses and for the same reason: no
    // cloud, or Rb = 0, means the gather chain would scatter and collapse a
    // buffer of zeros into `b` for nothing.
    if (this.particles && this.Rb !== 0) {
      this.dispatch(p, "cq", this.nrq * this.naq);
      this.dispatch(p, "gC", this.nrq * this.na);
      this.dispatch(p, "bC", this.nr * this.na);
    }
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
    this.dispatch(p, "rms", 1);      // velocity balancing the T this ψ was solved from
    this.dispatch(p, "cfl", 1);      // CFL speed of that same velocity, for the next dt
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
    if (!this.o.tackley && !this.o.tosi && !this.o.blankenbach && !this.o.vanKeken)
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
      // Tracers push on the same ψ the temperature transport above just
      // used — the Stokes solve below re-solves ψ from the T this step
      // arrived at, so pushing first (still inside this pass, dispatches
      // within it are ordered) is what keeps a tracer's path and the field
      // it is drawn over showing the same instant of the flow. The
      // composition projection only runs when it can actually change the
      // load: at Rb = 0 the buoyancy load's compositional gather is never
      // even dispatched (see `stokes`), so scattering tracers onto a grid
      // nothing samples would be pure waste.
      if (this.particles) {
        this.particles.push(p);
        if (this.Rb !== 0) this.particles.project(p);
      }
      this.dispatch(p, "nusselt", 1);
      p.end();
      // psiMax and rms run inside `stokes`, right after ψ is written; stokes
      // itself is what reads the composition grid back through `tqSource`.
      this.stokes(enc);
    });
    this.time += this.dt;
    this.steps++;
  }

  /** Draw the current temperature field straight from its storage buffer, and the tracer overlay over it, if one is attached. */
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
      // Composited over the field by draw order, in the same pass — no
      // second attachment, no second submission.
      this.particles?.draw(p);
      p.end();
    });
  }

  /**
   * Refresh `stats` from the GPU-side reductions. Off the frame's dependency
   * chain: the copy is queued behind work already submitted and the map resolves
   * whenever it resolves, so the hot loop never blocks on it. `max|ψ|` is read
   * only for the HUD — the renderer takes it from the buffer directly.
   * `maxSpeed` is what `main.ts` sizes the adaptive step from — see
   * `adaptiveDt` — so this must be called every frame, not on the HUD's
   * cadence, or the step lags further behind the flow than the hysteresis band
   * assumes.
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
    const at = this.time, atStep = this.steps;
    this.encode((enc) =>
      enc.copyBufferToBuffer(this.buf.stat, 0, this.buf.statRead, 0, 5 * S));
    void this.buf.statRead.mapAsync(GPUMapMode.READ).then(() => {
      const a = new Float32Array(this.buf.statRead.getMappedRange().slice(0));
      this.buf.statRead.unmap();
      this.stats = {
        nuInner: a[0], nuOuter: a[1], psiMax: a[2], vrms: a[3], maxSpeed: a[4],
        at, atStep,
      };
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
