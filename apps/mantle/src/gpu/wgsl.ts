/**
 * WGSL sources for the compute pipeline.
 *
 * Every kernel here is a line-by-line port of the validated CPU reference — same
 * quadrature, same monotone bicubic, same BFECC sequence, same Thomas sweep — so
 * that the parity suite compares the two implementations rather than two
 * different schemes. Where the CPU carries f64, the GPU carries f32; the two
 * places where that would have mattered (the biharmonic factorisation and the
 * diffusion factorisation) are precomputed on the CPU in f64 and uploaded, so
 * only *applications* run in f32.
 *
 * Sources are built by functions rather than kept as `.wgsl` files because two
 * things are specialised at pipeline-creation time: the FFT on its transform
 * length — stage count, workgroup size and shared-array extents are all
 * compile-time constants in the emitted code — and the **metric** on the
 * geometry. The latter is a specialisation and not a uniform on purpose: in a
 * box the transverse metric is 1 and `1/r` would be evaluated at r = 0, so a
 * branch on a uniform would have to compute an infinity and multiply it by zero.
 * Emitting `hOf`/`ihOf`/`DH` per geometry means the box's kernels simply do not
 * contain the division. Changing geometry is a rebuild either way.
 *
 * Every binding declared in a source below is statically used by its entry
 * point, which is what makes `layout: "auto"` safe: the derived bind group
 * layout then matches the declaration order exactly (see `gpu/sim.ts`).
 */

import { COLORMAPS, type ColormapName } from "../colormaps";
import type { Geometry } from "../geometry";
import { EPS_MIN, TACKLEY_TRANSITION_DEPTH, TACKLEY_STRAIN_FLOOR } from "../solver/rheology";
import { CIC_FIXED_POINT_SCALE, type TintEntry } from "../particles";

/**
 * Uniform block shared by every kernel: 16 i32 then 28 f32, padded to 176 B — a
 * multiple of 16, which is what a uniform struct's size must round to. Layout is
 * mirrored by `I` and `F` in `gpu/sim.ts`, which is what the runtime controls
 * write through — `Ra`, `dt`, `levels`, `lineW`, `mesh`, `gamma`, `nExp`, `cz`
 * and the view transform (`zoom`, `panX`, `panY`) all change from the UI
 * without touching a pipeline.
 *
 * New fields are appended rather than slotted in beside a related one, so that
 * every existing index in `F` keeps its meaning: the JS side writes this block
 * by index, and a field inserted mid-struct would silently shift `mesh` and the
 * rest by one. `Rb`, `etaLight` and `etaDense` are what the particle feature
 * adds here — see the header on `PART` for why everything else it needs lives
 * in a uniform of its own instead. `Rb` is the compositional Rayleigh number:
 * the buoyancy load reads `Ra·T − Rb·C` rather than `Ra·(T − B·C)`
 * (`tqSource`/`bSource` for the first term, `cqSource`/`bcSource` for the
 * second), which is not a cosmetic difference — it is what lets `Rb` be
 * nonzero while `Ra = 0`, the purely compositional (isothermal) limit, which a
 * single ratio `B = Rb/Ra` could never express. `etaLight`/`etaDense` are the
 * two endpoints of the composition-linear viscosity law (`vanKekenMuSource`),
 * read only by that law's own kernel. All three reuse what were previously
 * padding floats, so the struct's size does not move.
 */
export const PARAMS = /* wgsl */ `
const P: i32 = 3;
const TAU: f32 = 6.28318530717959;

struct Params {
  nr: i32, na: i32, gnr: i32, gna: i32,
  ni: i32, gni: i32, nRq: i32, nAq: i32,
  rBase: i32, rNLast: i32, aBase: i32, aNLast: i32,
  k0: i32, nRel: i32, nAel: i32, ipad2: i32,
  ri: f32, ro: f32, dr: f32, dphi: f32,
  tIn: f32, tOut: f32, Ra: f32, dt: f32,
  aLo: f32, aLen: f32, fill: f32, levels: f32,
  lineW: f32, gamma: f32, nExp: f32, mesh: f32,
  cz: f32, zoom: f32, panX: f32, panY: f32,
  sigmaY: f32, sigmaB: f32, etaStar: f32, fpad0: f32,
  Rb: f32, etaLight: f32, etaDense: f32, fpad3: f32,
};
@group(0) @binding(0) var<uniform> pp: Params;
`;

/**
 * The metric of `src/geometry.ts`, as the three quantities every kernel that
 * knows about geometry needs: `h`, `1/h` and the constant `h′`.
 *
 * `pp.ri`/`pp.ro` are the limits of the non-periodic axis in both cases — the
 * two radii, or z = 0 and z = 1 — and `pp.aLo`/`pp.aLen` already carry the
 * period of the transverse one, so nothing else about the domain needs to reach
 * the shaders separately.
 */
export const metric = (g: Geometry): string => g.kind === "annulus"
  ? /* wgsl */ `
fn hOf(r: f32) -> f32 { return r; }
fn ihOf(r: f32) -> f32 { return 1.0 / r; }
const DH: f32 = 1.0;
`
  : /* wgsl */ `
fn hOf(r: f32) -> f32 { return 1.0; }
fn ihOf(r: f32) -> f32 { return 1.0; }
const DH: f32 = 0.0;
`;

/**
 * Close a tracer's transverse coordinate onto one period, the same way
 * `metric` closes the radial one: periodic domains — the annulus, and a box
 * with periodic side walls — just wrap, exactly as `wrapPhi` already does for
 * every spline evaluation. A **free-slip** box does not: x = 0 and x = width
 * are solid walls, so nothing in the continuous flow ever reaches them with
 * outward speed, and a tracer that numerically overshoots one by the push's
 * own truncation error belongs reflected back in, not wrapped around into
 * what would be its own mirror image. Requires `wrapPhi` (see `BASIS`) in
 * scope for the periodic branch, and needs `pp.aLo`/`pp.aLen` either way — it
 * is emitted per `walls`, next to `metric`, for the reason `metric` is
 * emitted per `kind`: the two are the same kind of choice, made about the two
 * different axes.
 */
export const foldPhi = (g: Geometry): string => g.walls === "free-slip"
  ? /* wgsl */ `
fn foldPhi(phi: f32) -> f32 {
  let u = (((phi - pp.aLo) % pp.aLen) + pp.aLen) % pp.aLen;
  let w = pp.aLen * 0.5;
  return pp.aLo + select(u, pp.aLen - u, u > w);
}
`
  : /* wgsl */ `
fn foldPhi(phi: f32) -> f32 { return wrapPhi(phi); }
`;

/** Catmull–Rom clamped to its bracketing values — `cubic` in temperature.ts. */
const CUBIC = /* wgsl */ `
fn cubic(p0: f32, p1: f32, p2: f32, p3: f32, t: f32) -> f32 {
  let v = p1 + 0.5 * t * (p2 - p0 + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3
    + t * (3.0 * (p1 - p2) + p3 - p0)));
  return clamp(v, min(p1, p2), max(p1, p2));
}
`;

/**
 * A grid's shape, as the four expressions `cell()` needs — generalised so the
 * one interpolation scheme reads the (coarser) composition grid as readily as
 * it reads the temperature grid. `T_GRID` is `Temperature`'s own grid, fixed
 * by the resolution and so read straight off `pp`; `C_GRID` is the particle
 * feature's composition grid (`particles.ts`), whose node count lives in
 * `PART` rather than `PARAMS` — see that struct's header — and whose spacing
 * is not stored at all, only derived from the domain's own extent at the
 * count `PART` gives it, the same relation `pp.dr`/`pp.dphi` already encode
 * for the temperature grid.
 */
interface Grid { tag: string; nr: string; na: string; dr: string; dphi: string }

const T_GRID: Grid = { tag: "T", nr: "pp.gnr", na: "pp.gna", dr: "pp.dr", dphi: "pp.dphi" };
const C_GRID: Grid = {
  tag: "C", nr: "pt.cnr", na: "pt.cna",
  dr: "((pp.ro - pp.ri) / f32(pt.cnr - 1))", dphi: "(pp.aLen / f32(pt.cna))",
};

/** Grid cell containing (r, φ): r clamped to the domain, φ wrapped. */
const cellFn = (g: Grid) => /* wgsl */ `
struct Cell${g.tag} { i: i32, j: i32, tr: f32, tp: f32 };

fn cell${g.tag}(r: f32, phi: f32) -> Cell${g.tag} {
  let x = clamp((r - pp.ri) / ${g.dr}, 0.0, f32(${g.nr} - 1));
  let i = min(${g.nr} - 2, i32(floor(x)));
  var y = phi / ${g.dphi};
  y = y - floor(y / f32(${g.na})) * f32(${g.na});
  let j = i32(floor(y));
  return Cell${g.tag}(i, j, x - f32(i), y - f32(j));
}
`;

const CELL = cellFn(T_GRID);

/** Monotone bicubic sample of a grid buffer — `Temperature.sample`, generalised over which grid (see `Grid`). */
const sampleFn = (buf: string, g: Grid = T_GRID) => /* wgsl */ `
fn sample_${buf}(r: f32, phi: f32) -> f32 {
  let c = cell${g.tag}(r, phi);
  let na = ${g.na};
  var col: array<f32, 4>;
  for (var m = -1; m <= 2; m++) {
    let row = clamp(c.i + m, 0, ${g.nr} - 1) * na;
    col[m + 1] = cubic(${buf}[row + (c.j - 1 + na) % na], ${buf}[row + c.j % na],
                       ${buf}[row + (c.j + 1) % na], ${buf}[row + (c.j + 2) % na], c.tp);
  }
  return cubic(col[0], col[1], col[2], col[3], c.tr);
}
`;

/** Local extrema of the containing cell — the BFECC limiter bracket. Always the temperature grid: BFECC is a T-transport device only. */
const bracketFn = (buf: string) => /* wgsl */ `
fn bracket_${buf}(r: f32, phi: f32) -> vec2f {
  let c = cellT(r, phi);
  var lo = 1e30; var hi = -1e30;
  for (var m = 0; m <= 1; m++) {
    for (var l = 0; l <= 1; l++) {
      let v = ${buf}[(c.i + m) * pp.gna + (c.j + l) % pp.gna];
      lo = min(lo, v); hi = max(hi, v);
    }
  }
  return vec2f(lo, hi);
}
`;

/**
 * B-spline basis and derivatives — Piegl & Tiller A2.1/A2.3, specialised to
 * p = 3. `knots` holds both axes end to end; `base` selects one.
 *
 * The `ndu` table is built once per abscissa and the derivatives are read off
 * it, because it already holds every lower-degree basis value the derivative
 * recurrence
 *
 *   N'_{i,p} = p [ N_{i,p−1}/(u_{i+p} − u_i) − N_{i+1,p−1}/(u_{i+p+1} − u_{i+1}) ]
 *
 * needs: `m[r][j]` is the r-th nonzero degree-j value, and `m[j][r]` the knot
 * width that recurrence divides by. `val2` is that same recurrence applied a
 * second time, at degree 2 and then at degree 3 — which is why the general
 * algorithm's `a`-table loop, whose bounds collapse to the two end terms for
 * every r at these degrees, is not needed here.
 */
const BASIS = /* wgsl */ `
struct ND { span: i32, m: array<array<f32, 4>, 4> };
struct B1 { span: i32, n0: vec4f, n1: vec4f };

fn wrapPhi(u: f32) -> f32 {
  return pp.aLo + (((u - pp.aLo) % pp.aLen) + pp.aLen) % pp.aLen;
}

fn findSpan(base: i32, nLast: i32, u: f32) -> i32 {
  if (u >= knots[base + nLast + 1]) { return nLast; }
  if (u <= knots[base + P]) { return P; }
  var lo = P; var hi = nLast + 1; var mid = (lo + hi) / 2;
  loop {
    if (u >= knots[base + mid] && u < knots[base + mid + 1]) { break; }
    if (u < knots[base + mid]) { hi = mid; } else { lo = mid; }
    mid = (lo + hi) / 2;
  }
  return mid;
}

fn ndu(base: i32, nLast: i32, u: f32) -> ND {
  let span = findSpan(base, nLast, u);
  var m: array<array<f32, 4>, 4>;
  var left: array<f32, 4>;
  var right: array<f32, 4>;
  m[0][0] = 1.0;
  for (var j = 1; j <= P; j++) {
    left[j] = u - knots[base + span + 1 - j];
    right[j] = knots[base + span + j] - u;
    var saved = 0.0;
    for (var r = 0; r < j; r++) {
      m[j][r] = right[r + 1] + left[j - r];
      let tmp = m[r][j - 1] / m[j][r];
      m[r][j] = saved + right[r + 1] * tmp;
      saved = left[j - r] * tmp;
    }
    m[j][j] = saved;
  }
  return ND(span, m);
}

fn val0(d: ND) -> vec4f { return vec4f(d.m[0][3], d.m[1][3], d.m[2][3], d.m[3][3]); }

fn val1(d: ND) -> vec4f {
  var o: array<f32, 4>;
  for (var r = 0; r <= P; r++) {
    var s = 0.0;
    if (r >= 1) { s = d.m[r - 1][2] / d.m[3][r - 1]; }
    if (r <= 2) { s = s - d.m[r][2] / d.m[3][r]; }
    o[r] = f32(P) * s;
  }
  return vec4f(o[0], o[1], o[2], o[3]);
}

fn val2(d: ND) -> vec4f {
  var e: array<f32, 3>;             // d/du of the three degree-2 basis functions
  for (var b = 0; b <= 2; b++) {
    var s = 0.0;
    if (b >= 1) { s = d.m[b - 1][1] / d.m[2][b - 1]; }
    if (b <= 1) { s = s - d.m[b][1] / d.m[2][b]; }
    e[b] = 2.0 * s;
  }
  var o: array<f32, 4>;
  for (var a = 0; a <= P; a++) {
    var s = 0.0;
    if (a >= 1) { s = e[a - 1] / d.m[3][a - 1]; }
    if (a <= 2) { s = s - e[a] / d.m[3][a]; }
    o[a] = f32(P) * s;
  }
  return vec4f(o[0], o[1], o[2], o[3]);
}

fn basis1(base: i32, nLast: i32, u: f32) -> B1 {
  let d = ndu(base, nLast, u);
  return B1(d.span, val0(d), val1(d));
}
`;

/** ψ itself, for the streamline contours. */
const PSI_AT = /* wgsl */ `
fn psiAt(r: f32, phi: f32) -> f32 {
  let R = basis1(pp.rBase, pp.rNLast, r);
  let A = basis1(pp.aBase, pp.aNLast, wrapPhi(phi));
  var s = 0.0;
  for (var a = 0; a < 4; a++) {
    let row = (R.span - P + a) * pp.na;
    for (var b = 0; b < 4; b++) {
      s += R.n0[a] * A.n0[b] * psi[row + ((A.span - P + b) % pp.na + pp.na) % pp.na];
    }
  }
  return s;
}
`;

/**
 * u = ∇×(ψ ẑ) evaluated from the spline coefficients — `Field.velocity` — plus
 * the two ways of stepping a point through it: `departure` traces backward,
 * for semi-Lagrangian advection of a field; `advanceRK2` traces forward, for
 * pushing a particle along its own path. Both are RK2, and the one line that
 * tells them apart is not an oversight — see `advanceRK2`'s own comment, and
 * `particles.ts`, for why a difference that costs a one-cell hop nothing
 * matters once the same integrator runs for thousands of steps.
 */
const VELOCITY = /* wgsl */ `
fn velocity(r: f32, phi: f32) -> vec2f {
  let R = basis1(pp.rBase, pp.rNLast, r);
  let A = basis1(pp.aBase, pp.aNLast, wrapPhi(phi));
  var psi_r = 0.0; var psi_p = 0.0;
  for (var a = 0; a < 4; a++) {
    let row = (R.span - P + a) * pp.na;
    for (var b = 0; b < 4; b++) {
      let c = psi[row + ((A.span - P + b) % pp.na + pp.na) % pp.na];
      psi_r += R.n1[a] * A.n0[b] * c;
      psi_p += R.n0[a] * A.n1[b] * c;
    }
  }
  return vec2f(psi_p * ihOf(r), -psi_r);   // (u_r, u_φ)
}

/**
 * Backward RK2 characteristic trace; dt < 0 traces forward. The transverse
 * displacement is an *arc* — u_φ dt / h — which in a box is just u_x dt.
 */
fn departure(r: f32, phi: f32, dt: f32) -> vec2f {
  let ih = ihOf(r);
  let a = velocity(r, phi);
  let rm = r - 0.5 * dt * a.x;
  let pm = phi - 0.5 * dt * (a.y * ih);
  let b = velocity(clamp(rm, pp.ri, pp.ro), pm);
  return vec2f(r - dt * b.x, phi - dt * (b.y * ih));
}

/**
 * Forward explicit-midpoint RK2 trace of one tracer — the particle push.
 * Forward in the literal sense: unlike departure above, which traces
 * backward from a grid point to the value it advects from, this traces a
 * parcel of mantle forward to wherever the flow carries it next, which is
 * what a tracer's path (a pathline) actually is.
 *
 * The one line this changes from departure run at a negated dt is the final
 * transverse step, taken at the RK2 midpoint's h, not the starting radius's
 * (rm, not r, below) — free for a semi-Lagrangian hop one cell long, and not
 * for a trace integrated over thousands of steps. See particles.ts and the
 * CPU twin for the arithmetic behind why.
 */
fn advanceRK2(r: f32, phi: f32, dt: f32) -> vec2f {
  let a = velocity(r, phi);
  let rm = clamp(r + 0.5 * dt * a.x, pp.ri, pp.ro);
  let pm = phi + 0.5 * dt * (a.y * ihOf(r));
  let b = velocity(rm, pm);
  return vec2f(r + dt * b.x, phi + dt * (b.y * ihOf(rm)));   // ← rm, not r
}
`;

export const WG = 64; // workgroup size for the flat element-wise kernels

const flat = (n: string) => /* wgsl */ `
  let g = i32(gid.x);
  if (g >= ${n}) { return; }
`;

/**
 * Uniform block owned by the particle feature, not the shared `PARAMS` one:
 * the composition grid's node counts, the tracer count, the render-only
 * radius/opacity/canvas-side, and the clock the *age* colour mode reads. This
 * is `GpuParticles`' own, in the same way `Pass` (`advectSource`) is the
 * advection kernels' own rather than a field on `pp` — a tracer count or
 * render radius changing has nothing to do with the temperature/Stokes
 * pipelines, and giving it a separate small uniform is what keeps changing it
 * from ever touching theirs. `tqSource`'s buoyancy load binds it too, but only
 * for `cnr`/`cna`: the load needs to know the composition grid's *shape* to
 * sample it, not any of the render-only fields, and it reaches this same
 * struct rather than a third one so the shape can never disagree between the
 * kernel that writes C and the one that reads it.
 */
const PART = /* wgsl */ `
struct Part {
  count: i32, cnr: i32, cna: i32, ipad: i32,
  radius: f32, opacity: f32, side: f32, age: f32,
  tau: f32, fpad0: f32, fpad1: f32, fpad2: f32,
};
`;

// ---- buoyancy load: three gather passes (see solver/assembly.ts) -------------

/**
 * Tq[qr][qa] = T(r_qr, φ_qa): one bicubic sample per quadrature point.
 *
 * Pure T, always — this buffer is shared with the T-dependent viscosity laws'
 * own `muEval` (`muSource`, `tackleyMuSource`, …), which read it expecting a
 * value in [0, 1], so it must never carry anything else mixed in. The
 * compositional contribution to the buoyancy load is a *separate* gather
 * chain (`cqSource`/`bcSource`, below) for exactly that reason.
 */
export const tqSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> T: array<f32>;
@group(0) @binding(2) var<storage, read> rq: array<f32>;
@group(0) @binding(3) var<storage, read> phiq: array<f32>;
@group(0) @binding(4) var<storage, read_write> Tq: array<f32>;
` + CUBIC + CELL + sampleFn("T") + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let r = rq[g / pp.nAq]; let phi = phiq[g % pp.nAq];
  Tq[g] = sample_T(r, phi);
}
`;

/**
 * Cq[qr][qa] = C(r_qr, φ_qa) — the composition's own quadrature samples, the
 * twin of `tqSource` for the compositional half of the buoyancy load. Kept
 * out of `tqSource` itself (see that function's header) rather than folded
 * in as a second output, so nothing downstream has to know two fields are
 * being sampled from one dispatch.
 */
export const cqSource = () => PARAMS + PART + /* wgsl */ `
@group(0) @binding(1) var<storage, read> rq: array<f32>;
@group(0) @binding(2) var<storage, read> phiq: array<f32>;
@group(0) @binding(3) var<uniform> pt: Part;
@group(0) @binding(4) var<storage, read> C: array<f32>;
@group(0) @binding(5) var<storage, read_write> Cq: array<f32>;
` + CUBIC + cellFn(C_GRID) + sampleFn("C", C_GRID) + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let r = rq[g / pp.nAq]; let phi = phiq[g % pp.nAq];
  Cq[g] = sample_C(r, phi);
}
`;

/** G[qr][l] = Σ_t aVal[l][t] Tq[qr][aIdx[l][t]] — collapse the φ quadrature. */
export const gSource = (slots: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Tq: array<f32>;
@group(0) @binding(2) var<storage, read> aIdx: array<i32>;
@group(0) @binding(3) var<storage, read> aVal: array<f32>;
@group(0) @binding(4) var<storage, read_write> G: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.na")}
  let qr = g / pp.na; let l = g % pp.na;
  var s = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let b = l * ${slots} + t;
    s += aVal[b] * Tq[qr * pp.nAq + aIdx[b]];
  }
  G[g] = s;
}
`;

/** b[i][l] = Ra Σ_t rVal[i][t] G[rIdx[i][t]][l] — collapse the r quadrature. */
export const bSource = (slots: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> G: array<f32>;
@group(0) @binding(2) var<storage, read> rIdx: array<i32>;
@group(0) @binding(3) var<storage, read> rVal: array<f32>;
@group(0) @binding(4) var<storage, read_write> b: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let i = g / pp.na; let l = g % pp.na;
  // psi = const (or psi = psi_r = 0, no-slip) is essential, so margin boundary
  // DOFs per end are not in the trial space -- derived from nr/ni rather than
  // a new uniform (see geometry.ts's radialMargin). Tier 1 simply ignores
  // these rows; the Krylov tier would otherwise carry a residual component
  // its operator can never remove.
  let margin = (pp.nr - pp.ni) / 2;
  if (i < margin || i >= pp.nr - margin) { b[g] = 0.0; return; }
  var s = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let q = i * ${slots} + t;
    s += rVal[q] * G[rIdx[q] * pp.na + l];
  }
  b[g] = pp.Ra * s;
}
`;

/**
 * b[i][l] −= Rb Σ_t rVal[i][t] GC[rIdx[i][t]][l] — the compositional half of
 * the load, subtracted into the buffer `bSource` already wrote. Two
 * independently scaled dispatches into one buffer rather than one dispatch
 * reading a pre-mixed `Ra·T − Rb·C`, because `Ra` and `Rb` have to be free to
 * move independently — `Rb` alone, with `Ra = 0`, is the purely compositional
 * (isothermal) buoyancy the Rayleigh–Taylor case needs, and a single ratio
 * `B = Rb/Ra` can never express that limit. Dispatched only when a tracer
 * cloud is attached and actually coupled (`GpuSimulation.stokes`); the sign
 * is the same one `tqSource`'s retired `T − buoy·C` used, so a dense layer
 * still subtracts from the load exactly as a cold anomaly does.
 */
export const bcSource = (slots: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> GC: array<f32>;
@group(0) @binding(2) var<storage, read> rIdx: array<i32>;
@group(0) @binding(3) var<storage, read> rVal: array<f32>;
@group(0) @binding(4) var<storage, read_write> b: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let i = g / pp.na; let l = g % pp.na;
  let margin = (pp.nr - pp.ni) / 2;
  if (i < margin || i >= pp.nr - margin) { return; }
  var s = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let q = i * ${slots} + t;
    s += rVal[q] * GC[rIdx[q] * pp.na + l];
  }
  b[g] -= pp.Rb * s;
}
`;

// ---- FFT --------------------------------------------------------------------

/**
 * Stockham autosort radix-2 FFT, one workgroup per row, entirely in shared
 * memory: no bit-reversal pass and no global ping-pong. The two halves of a
 * 2N-element shared array alternate as source and destination, so `logN` stages
 * need one barrier each and nothing else.
 *
 * Convention matches `src/dft.ts`: forward carries the 1/N, inverse does not.
 */
const fftBody = (n: number, sign: number, scale: string) => {
  const half = n / 2, wg = Math.min(half, 256);
  const logn = Math.log2(n);
  if (!Number.isInteger(logn)) throw new Error(`FFT length ${n} is not a power of two`);
  return /* wgsl */ `
const N: i32 = ${n};
const HALF: i32 = ${half};
const WGS: i32 = ${wg};
var<workgroup> xr: array<f32, ${2 * n}>;
var<workgroup> xi: array<f32, ${2 * n}>;

// Runs the stages; returns which half of the shared arrays holds the result.
fn transform(t: i32) -> i32 {
  var l = HALF; var m = 1; var src = 0;
  for (var stage = 0; stage < ${logn}; stage++) {
    workgroupBarrier();
    let dst = N - src;
    for (var p = t; p < HALF; p += WGS) {
      let j = p / m; let k = p % m;
      let i0 = src + k + j * m;
      let i1 = i0 + l * m;
      let o0 = dst + k + 2 * j * m;
      let ang = ${sign === 1 ? "" : "-"}3.14159265358979 * f32(j) / f32(l);
      let w = vec2f(cos(ang), sin(ang));
      let d = vec2f(xr[i0] - xr[i1], xi[i0] - xi[i1]);
      xr[o0] = xr[i0] + xr[i1];
      xi[o0] = xi[i0] + xi[i1];
      xr[o0 + m] = w.x * d.x - w.y * d.y;
      xi[o0 + m] = w.x * d.y + w.y * d.x;
    }
    src = dst; l = l / 2; m = m * 2;
  }
  workgroupBarrier();
  return src;
}

const SCALE: f32 = ${scale};
`;
};

/** Real field (rows × N) → complex modes. */
export const fftForwardSource = (n: number) => /* wgsl */ `
@group(0) @binding(0) var<storage, read> src: array<f32>;
@group(0) @binding(1) var<storage, read_write> outRe: array<f32>;
@group(0) @binding(2) var<storage, read_write> outIm: array<f32>;
` + fftBody(n, -1, `1.0 / ${n}.0`) + `
@compute @workgroup_size(${Math.min(n / 2, 256)})
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_id) lid: vec3u) {
  let row = i32(wid.x) * N;
  let t = i32(lid.x);
  for (var i = t; i < N; i += WGS) { xr[i] = src[row + i]; xi[i] = 0.0; }
  let o = transform(t);
  for (var i = t; i < N; i += WGS) {
    outRe[row + i] = xr[o + i] * SCALE;
    outIm[row + i] = xi[o + i] * SCALE;
  }
}
`;

/** Complex modes → real field (rows × N). */
export const fftInverseSource = (n: number) => /* wgsl */ `
@group(0) @binding(0) var<storage, read> srcRe: array<f32>;
@group(0) @binding(1) var<storage, read> srcIm: array<f32>;
@group(0) @binding(2) var<storage, read_write> dst: array<f32>;
` + fftBody(n, 1, "1.0") + `
@compute @workgroup_size(${Math.min(n / 2, 256)})
fn main(@builtin(workgroup_id) wid: vec3u, @builtin(local_invocation_id) lid: vec3u) {
  let row = i32(wid.x) * N;
  let t = i32(lid.x);
  for (var i = t; i < N; i += WGS) { xr[i] = srcRe[row + i]; xi[i] = srcIm[row + i]; }
  let o = transform(t);
  for (var i = t; i < N; i += WGS) { dst[row + i] = xr[o + i] * SCALE; }
}
`;

// ---- Stokes: per-mode radial solve ------------------------------------------

/**
 * ψ̂_k = A_k⁻¹ b̂_k, a dense matvec against the f64-computed inverse.
 * Only `k = 0 … na/2` are stored (`A_k = A_{n−k}`); `margin` boundary DOFs per
 * end are zero (ψ = const, free-slip; ψ = ψ_r = 0, no-slip — `margin` derived
 * from `pp.nr`/`pp.ni`, same as `bSource`).
 *
 * The mean mode is governed by `pp.k0`. Tier 1 zeroes it: constant μ does not
 * couple the modes and the buoyancy load has no mean component (∮∂_φT dφ = 0),
 * so ψ̂₀ ≡ 0 exactly, and *solving* it would only amplify the f32 transform's
 * round-off by κ(A_0). Tier 2 solves it, because this kernel is then the
 * preconditioner and variable μ gives the residual genuine mean content. The
 * block is nonsingular either way — see `modeInverses`.
 */
export const radialSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> bRe: array<f32>;
@group(0) @binding(2) var<storage, read> bIm: array<f32>;
@group(0) @binding(3) var<storage, read> inv: array<f32>;
@group(0) @binding(4) var<storage, read_write> outRe: array<f32>;
@group(0) @binding(5) var<storage, read_write> outIm: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let i = g / pp.na; let k = g % pp.na;
  let margin = (pp.nr - pp.ni) / 2;
  if (i < margin || i >= pp.nr - margin || (k == 0 && pp.k0 == 0)) {
    outRe[g] = 0.0; outIm[g] = 0.0; return;
  }
  let base = min(k, pp.na - k) * pp.ni * pp.ni + (i - margin) * pp.ni;
  var sr = 0.0; var si = 0.0;
  for (var c = 0; c < pp.ni; c++) {
    let a = inv[base + c];
    sr += a * bRe[(c + margin) * pp.na + k];
    si += a * bIm[(c + margin) * pp.na + k];
  }
  outRe[g] = sr; outIm[g] = si;
}
`;

// ---- variable-μ operator and conjugate gradients -----------------------------

/**
 * `(∂_φA[ψ], C[ψ])` at one quadrature point, with
 * `A[ψ] = ψ_r/h − h′ψ/h²` and `C[ψ] = ψ_rr − (h′/h)ψ_r − ψ_φφ/h²` — the twin of
 * `deformation` in `solver/assembly.ts`, and shared here for the same reason: it
 * is both the first pass of the operator apply and, unweighted, the strain-rate
 * tensor the power law reads. Both kernels below name their coefficient buffer
 * `c`.
 *
 * The basis is recomputed per quadrature point rather than tabulated. That is a
 * deliberate trade: tables for `(a, g, b)` and `(N, N', N'')` would need six more
 * storage bindings, and WebGPU guarantees only eight per stage — the operator
 * kernel would not fit. Recomputing costs two `ndu` builds against a 16-term
 * tensor sum, which on a GPU is the cheap side of the trade anyway.
 */
const DEFORM = /* wgsl */ `
fn deformation(r: f32, phi: f32) -> vec2f {
  let R = ndu(pp.rBase, pp.rNLast, r);
  let A = ndu(pp.aBase, pp.aNLast, phi);
  let ih = ihOf(r); let ih2 = ih * ih;
  let r0 = val0(R); let r1 = val1(R);
  let ra = r1 * ih - DH * r0 * ih2;   // pairs with N'
  let rg = val2(R) - DH * r1 * ih;    // pairs with N
  let rb = r0 * ih2;                  // pairs with −N''
  let a0 = val0(A); let a1 = val1(A); let a2 = val2(A);

  var d = vec2f(0.0, 0.0);
  for (var p = 0; p < 4; p++) {
    let row = (R.span - P + p) * pp.na;
    for (var q = 0; q < 4; q++) {
      let v = c[row + ((A.span - P + q) % pp.na + pp.na) % pp.na];
      d.x += ra[p] * a1[q] * v;
      d.y += (rg[p] * a0[q] - rb[p] * a2[q]) * v;
    }
  }
  return d;
}
`;

/**
 * Quadrature-point stresses of the variable-μ dissipation form — pass 1 of the
 * matrix-free apply, the twin of `applyOperator` in `solver/assembly.ts`:
 *
 *   Pq = 4 μ h ∂_φA[ψ],   Qq = μ h C[ψ].
 *
 * μ is **read**, not evaluated: under μ(T) alone it could be exponentiated
 * inline, but the power law makes it a function of ψ, and this
 * kernel runs once per Krylov iteration against the *search direction*. Taking μ
 * from a buffer the frame filled once (`muSource`) is what keeps the operator
 * CG sees linear and symmetric — evaluating it here would silently make it
 * neither.
 */
export const qevalSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> knots: array<f32>;
@group(0) @binding(2) var<storage, read> c: array<f32>;
@group(0) @binding(3) var<storage, read> rq: array<f32>;
@group(0) @binding(4) var<storage, read> phiq: array<f32>;
@group(0) @binding(5) var<storage, read> mu: array<f32>;
@group(0) @binding(6) var<storage, read_write> Pq: array<f32>;
@group(0) @binding(7) var<storage, read_write> Qq: array<f32>;
` + metric(geom) + BASIS + DEFORM + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let r = rq[g / pp.nAq];
  let d = deformation(r, phiq[g % pp.nAq]);
  let m = mu[g] * hOf(r);
  Pq[g] = 4.0 * m * d.x;
  Qq[g] = m * d.y;
}
`;

/**
 * Pass 1 of the rheology update: `ε̇_II = √(ε_rr² + ε_rφ²)` from ψ, written
 * into the μ buffer that the next two kernels transform in place. `ε_rr = ∂_φA`
 * and `ε_rφ = −½C`, so this is `deformation` and one `length` — the power law
 * costs one dispatch and no new geometry.
 */
export const strainSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> knots: array<f32>;
@group(0) @binding(2) var<storage, read> c: array<f32>;
@group(0) @binding(3) var<storage, read> rq: array<f32>;
@group(0) @binding(4) var<storage, read> phiq: array<f32>;
@group(0) @binding(5) var<storage, read_write> mu: array<f32>;
` + metric(geom) + BASIS + DEFORM + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let d = deformation(rq[g / pp.nAq], phiq[g % pp.nAq]);
  mu[g] = length(vec2f(d.x, 0.5 * d.y));
}
`;

/**
 * Pass 2: the strain-rate normalisation — `δ = ε̇_min·rms(ε̇)` into `sc[4]`
 * and the geometric mean of `ε̇ + δ` into `sc[5]`. The twin of `strainScale`.
 *
 * Reducing on the GPU is not an optimisation but a requirement: the law's
 * normalisation is a global quantity, and reading it back would be exactly the
 * per-frame host round-trip the frame loop forbids.
 *
 * **One dispatch, two reductions**, because the second needs the first: δ sets
 * the offset that keeps `log ε̇` finite at a stagnation point. A workgroup
 * barrier between them is enough — splitting into two kernels would mean two
 * dispatches and a round trip through the scalar buffer for a value 256 threads
 * already hold.
 *
 * A zero field has no scale, and gives (1, 1) — normalised strain ≡ 1, the
 * unstrained case. Reachable only on the first solve after ψ is cleared.
 */
export const srefSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> e: array<f32>;
@group(0) @binding(2) var<storage, read_write> sc: array<f32>;

var<workgroup> sm: array<f32, 256>;
var<workgroup> scale: vec2f;

fn reduce(t: u32, v: f32) -> f32 {
  sm[t] = v;
  for (var k = 128u; k > 0u; k >>= 1u) {
    workgroupBarrier();
    if (t < k) { sm[t] += sm[t + k]; }
  }
  workgroupBarrier();
  return sm[0];
}

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  let n = pp.nRq * pp.nAq;

  var s = 0.0;
  for (var i = i32(t); i < n; i += 256) { s += e[i] * e[i]; }
  let rms = sqrt(reduce(t, s) / f32(n));
  if (t == 0u) { scale = select(vec2f(1.0, 0.0), vec2f(${EPS_MIN} * rms, 1.0), rms > 0.0); }
  workgroupBarrier();
  let d = scale.x;

  var l = 0.0;
  for (var i = i32(t); i < n; i += 256) { l += log(e[i] + d); }
  let gm = exp(reduce(t, l) / f32(n));
  if (t == 0u) { sc[4] = d; sc[5] = select(1.0, gm, scale.y > 0.0); }
}
`;

/**
 * Pass 3: the law itself, in place over the strain rates —
 *
 *   μ = clamp( exp(−γ(T−½) + c(d−½)) ((ε̇+δ)/G)^((1−n)/n),
 *              e^{−(γ+c)/2}, e^{+(γ+c)/2} ),
 *
 * the twin of `viscosity` in `solver/rheology.ts`. `Tq` is the bicubic sample
 * the buoyancy assembly already needed, so T costs nothing here either, and the
 * depth `d = (r_o − r)/(r_o − r_i)` comes from the radial abscissa the strain
 * pass already binds — so the depth term adds one buffer read and no dispatch.
 *
 * At `nExp = 1` the exponent is zero and the power-law factor is exactly 1,
 * while the thermal and depth terms already lie inside the clamp — so this
 * reduces to the μ(T, d) law bit for bit, and switching between the two laws is
 * a uniform write rather than a rebuild.
 */
export const muSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Tq: array<f32>;
@group(0) @binding(2) var<storage, read> sc: array<f32>;
@group(0) @binding(3) var<storage, read> rq: array<f32>;
@group(0) @binding(4) var<storage, read_write> mu: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let hi = exp(0.5 * (pp.gamma + pp.cz));
  let d = (pp.ro - rq[g / pp.nAq]) / (pp.ro - pp.ri);
  let m = exp(-pp.gamma * (clamp(Tq[g], 0.0, 1.0) - 0.5) + pp.cz * (d - 0.5))
    * pow((mu[g] + sc[4]) / sc[5], (1.0 - pp.nExp) / pp.nExp);
  mu[g] = clamp(m, 1.0 / hi, hi);
}
`;

/**
 * The Tackley pseudo-plastic law, the twin of `tackleyViscosity` in
 * `solver/rheology.ts` — an Arrhenius branch and a Bingham (yield-stress)
 * branch combined as resistors in parallel. `mu[g]` holds raw ε̇_II from
 * `strainSource`, used directly rather than normalised: the yield stress is an
 * absolute threshold, so there is no `sref` pass and no reduction to wait on.
 */
export const tackleyMuSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Tq: array<f32>;
@group(0) @binding(2) var<storage, read> rq: array<f32>;
@group(0) @binding(3) var<storage, read_write> mu: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let T = clamp(Tq[g], 0.0, 1.0);
  let d = (pp.ro - rq[g / pp.nAq]) / (pp.ro - pp.ri);
  let A0 = select(1.0, 30.0, d >= ${TACKLEY_TRANSITION_DEPTH});
  let etaLin = A0 * exp(27.631 / (T + 0.88)) * 5.86052e-13;
  let etaPlast = pp.etaStar
    + (pp.sigmaY + pp.sigmaB * d) / max(mu[g], ${TACKLEY_STRAIN_FLOOR});
  mu[g] = 1.0 / (1.0 / etaLin + 1.0 / etaPlast);
}
`;

/**
 * Blankenbach et al. (1989)'s own law, the twin of `blankenbachViscosity` in
 * `solver/rheology.ts` — `pp.gamma`/`pp.cz` read as the paper's own b, c, the
 * same two uniforms μ(T, d) and Tosi read. No strain rate anywhere in it, so
 * unlike those two this kernel does not even read `mu[g]` — `strainSource`
 * still runs before it (see `rheology()` in `gpu/sim.ts`), its result simply
 * unused, since skipping the dispatch itself would need a further branch in
 * the tier's own dispatch order for one kernel this cheap to save.
 */
export const blankenbachMuSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Tq: array<f32>;
@group(0) @binding(2) var<storage, read> rq: array<f32>;
@group(0) @binding(3) var<storage, read_write> mu: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let T = clamp(Tq[g], 0.0, 1.0);
  let d = (pp.ro - rq[g / pp.nAq]) / (pp.ro - pp.ri);
  mu[g] = exp(-pp.gamma * T + pp.cz * d);
}
`;

/**
 * van Keken et al. (1997)'s composition-linear law, the twin of
 * `vanKekenViscosity` in `solver/rheology.ts` —
 *
 *   μ(φ) = η_light + φ (η_dense − η_light),
 *
 * a plain linear interpolation between the two materials' own viscosities by
 * the composition a tracer cloud has projected onto this quadrature point.
 * Reads `C` directly (`sample_C`, the same bicubic composition read
 * `cqSource` and the buoyancy load use) rather than `Tq`: this law has no T
 * dependence at all, so it is the one muEval variant that never touches the
 * buffer the T-dependent laws share — see `tqSource`'s own header on why
 * that separation matters. `strainSource` still runs before it, unused, for
 * the same reason `blankenbachMuSource` leaves it unused: skipping the
 * dispatch would need a further branch in `rheology()` for one kernel this
 * cheap to save.
 */
export const vanKekenMuSource = () => PARAMS + PART + /* wgsl */ `
@group(0) @binding(1) var<storage, read> rq: array<f32>;
@group(0) @binding(2) var<storage, read> phiq: array<f32>;
@group(0) @binding(3) var<uniform> pt: Part;
@group(0) @binding(4) var<storage, read> C: array<f32>;
@group(0) @binding(5) var<storage, read_write> mu: array<f32>;
` + CUBIC + cellFn(C_GRID) + sampleFn("C", C_GRID) + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let r = rq[g / pp.nAq]; let phi = phiq[g % pp.nAq];
  let c = sample_C(r, phi);
  mu[g] = pp.etaLight + c * (pp.etaDense - pp.etaLight);
}
`;

/**
 * The Tosi et al. (2015) viscoplastic law, the twin of `tosiViscosity` in
 * `solver/rheology.ts` — the μ(T, d) exponential (unclamped) in parallel with
 * the same Bingham branch Tackley's law states, combined as the paper's own
 * harmonic average (hence the factor of 2, unlike Tackley's plain
 * parallel-resistor sum). `b`/`c` — γ_T/γ_z in the paper — are read from
 * `pp.gamma`/`pp.cz`, the same two contrast uniforms the power law reads,
 * since Tosi's contrast sliders are those sliders, not a second pair. As
 * with Tackley, `mu[g]` holds raw ε̇_II from `strainSource`: the yield stress
 * is an absolute threshold, so there is no `sref` pass here either.
 */
export const tosiMuSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Tq: array<f32>;
@group(0) @binding(2) var<storage, read> rq: array<f32>;
@group(0) @binding(3) var<storage, read_write> mu: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  let T = clamp(Tq[g], 0.0, 1.0);
  let d = (pp.ro - rq[g / pp.nAq]) / (pp.ro - pp.ri);
  let etaLin = exp(-pp.gamma * (T - 0.5) + pp.cz * (d - 0.5));
  let etaPlast = pp.etaStar
    + (pp.sigmaY + pp.sigmaB * d) / max(mu[g], ${TACKLEY_STRAIN_FLOOR});
  mu[g] = 2.0 / (1.0 / etaLin + 1.0 / etaPlast);
}
`;

/**
 * Pass 2: collapse the φ quadrature against N', N and N''. `aOp` interleaves the
 * three azimuthal gather tables and `Gop` holds the three results as consecutive
 * planes — packing that keeps the kernel inside the eight-storage-buffer limit
 * without a second dispatch.
 */
export const gphiOpSource = (slots: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Pq: array<f32>;
@group(0) @binding(2) var<storage, read> Qq: array<f32>;
@group(0) @binding(3) var<storage, read> aIdx: array<i32>;
@group(0) @binding(4) var<storage, read> aOp: array<f32>;
@group(0) @binding(5) var<storage, read_write> Gop: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.na")}
  let qr = g / pp.na; let l = g % pp.na;
  var sa = 0.0; var sb = 0.0; var sc = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let j = l * ${slots} + t;
    let q = qr * pp.nAq + aIdx[j];
    sa += aOp[j * 3 + 1] * Pq[q];
    sb += aOp[j * 3 + 0] * Qq[q];
    sc += aOp[j * 3 + 2] * Qq[q];
  }
  let M = pp.nRq * pp.na;
  Gop[g] = sa; Gop[M + g] = sb; Gop[2 * M + g] = sc;
}
`;

/** Pass 3: collapse the r quadrature against a, g and −b, giving A ψ. */
export const grOpSource = (slots: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> Gop: array<f32>;
@group(0) @binding(2) var<storage, read> rIdx: array<i32>;
@group(0) @binding(3) var<storage, read> rOp: array<f32>;
@group(0) @binding(4) var<storage, read_write> out: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let i = g / pp.na; let l = g % pp.na;
  // Zeroing the output rows as well as the input's is what keeps the operator
  // symmetric, which conjugate gradients requires. margin boundary rows per
  // end, derived the same way bSource/radialSource do.
  let margin = (pp.nr - pp.ni) / 2;
  if (i < margin || i >= pp.nr - margin) { out[g] = 0.0; return; }
  let M = pp.nRq * pp.na;
  var s = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let j = i * ${slots} + t;
    let q = rIdx[j] * pp.na + l;
    s += rOp[j * 3 + 0] * Gop[q] + rOp[j * 3 + 1] * Gop[M + q]
       - rOp[j * 3 + 2] * Gop[2 * M + q];
  }
  out[g] = s;
}
`;

/**
 * ⟨a, b⟩ over the ψ space, reduced in a single workgroup into `sc[slot]`.
 *
 * The scalar stays on the GPU: reading α and β back to the host would put a
 * synchronisation in the middle of every Krylov iteration, which is exactly the
 * readback the frame loop forbids. The consuming kernels divide by these slots
 * themselves — the dispatches are ordered within a compute pass, so no explicit
 * barrier is needed.
 */
export const dotSource = (slot: number) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> a: array<f32>;
@group(0) @binding(2) var<storage, read> b: array<f32>;
@group(0) @binding(3) var<storage, read_write> sc: array<f32>;

var<workgroup> sm: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  var s = 0.0;
  for (var i = i32(t); i < pp.nr * pp.na; i += 256) { s += a[i] * b[i]; }
  sm[t] = s;
  for (var k = 128u; k > 0u; k >>= 1u) {
    workgroupBarrier();
    if (t < k) { sm[t] += sm[t + k]; }
  }
  if (t == 0u) { sc[${slot}] = sm[0]; }
}
`;

/**
 * res = b − A x — the initial residual of a warm-started solve, and the only
 * *true* residual the loop forms (see `krylov` in sim.ts on why the recursive
 * one must not be reported).
 */
export const cgInitSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> b: array<f32>;
@group(0) @binding(2) var<storage, read> Ax: array<f32>;
@group(0) @binding(3) var<storage, read_write> res: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  res[g] = b[g] - Ax[g];
}
`;

/** dst = src. Sets the first search direction to the preconditioned residual. */
export const cgCopySource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> src: array<f32>;
@group(0) @binding(2) var<storage, read_write> dst: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  dst[g] = src[g];
}
`;

/**
 * x += α p, res −= α Ap, with α = ⟨r,z⟩/⟨p,Ap⟩ read from the scalar buffer.
 *
 * The guard is not defensive padding: at γ = 0 the preconditioner *is* the
 * operator, so a warm-started solve empties the Krylov space on its second
 * iteration and ⟨p,Ap⟩ is legitimately zero. Without it that iteration would
 * write NaN into ψ and the isoviscous limiting case — the one check that the two
 * tiers agree — would fail for a reason that has nothing to do with the physics.
 */
export const cgUpdateXSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> sc: array<f32>;
@group(0) @binding(2) var<storage, read> dir: array<f32>;
@group(0) @binding(3) var<storage, read> Ap: array<f32>;
@group(0) @binding(4) var<storage, read_write> x: array<f32>;
@group(0) @binding(5) var<storage, read_write> res: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let a = select(0.0, sc[0] / sc[1], sc[1] > 0.0 && sc[0] > 0.0);
  x[g] += a * dir[g];
  res[g] -= a * Ap[g];
}
`;

/** p = z + β p, with β = ⟨r,z⟩_new/⟨r,z⟩_old. */
export const cgUpdatePSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> sc: array<f32>;
@group(0) @binding(2) var<storage, read> z: array<f32>;
@group(0) @binding(3) var<storage, read_write> dir: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nr * pp.na")}
  let b = select(0.0, sc[2] / sc[0], sc[0] > 0.0);
  dir[g] = z[g] + b * dir[g];
}
`;

/**
 * ⟨r,z⟩_old ← ⟨r,z⟩_new, closing the iteration. A separate dispatch rather than a
 * line inside `cgUpdatePSource`, where the threads still reading slot 0 would
 * race the one writing it.
 */
export const cgRollSource = () => /* wgsl */ `
@group(0) @binding(0) var<storage, read_write> sc: array<f32>;

@compute @workgroup_size(1)
fn main() { sc[0] = sc[2]; }
`;

// ---- temperature transport --------------------------------------------------

/**
 * One semi-Lagrangian pass. `limiter` additionally clamps to the local extrema
 * of `brk` at the departure point (the BFECC monotone limiter); `bc` selects the
 * isothermal boundary rows (1) or zeros them (0, the reverse BFECC pass, whose
 * boundary rows the CPU reference also leaves at zero before correcting).
 */
export const advectSource = (geom: Geometry) => PARAMS + /* wgsl */ `
struct Pass { dt: f32, limiter: i32, bc: i32, pad: i32 };
@group(0) @binding(1) var<uniform> ps: Pass;
@group(0) @binding(2) var<storage, read> knots: array<f32>;
@group(0) @binding(3) var<storage, read> psi: array<f32>;
@group(0) @binding(4) var<storage, read> src: array<f32>;
@group(0) @binding(5) var<storage, read> brk: array<f32>;
@group(0) @binding(6) var<storage, read_write> dst: array<f32>;
` + metric(geom) + CUBIC + CELL + sampleFn("src") + bracketFn("brk")
  + BASIS + VELOCITY + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.gnr * pp.gna")}
  let i = g / pp.gna; let j = g % pp.gna;
  if (i == 0 || i == pp.gnr - 1) {
    dst[g] = select(0.0, select(pp.tOut, pp.tIn, i == 0), ps.bc == 1);
    return;
  }
  let d = departure(pp.ri + f32(i) * pp.dr, f32(j) * pp.dphi, ps.dt);
  var v = sample_src(d.x, d.y);
  if (ps.limiter == 1) {
    let e = bracket_brk(d.x, d.y);
    v = clamp(v, e.x, e.y);
  }
  dst[g] = v;
}
`;

/** BFECC error correction: T + ½(T − T″), with the isothermal rows restored. */
export const bfeccSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> T: array<f32>;
@group(0) @binding(2) var<storage, read> T2: array<f32>;
@group(0) @binding(3) var<storage, read_write> dst: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.gnr * pp.gna")}
  let i = g / pp.gna;
  if (i == 0) { dst[g] = pp.tIn; return; }
  if (i == pp.gnr - 1) { dst[g] = pp.tOut; return; }
  dst[g] = T[g] + 0.5 * (T[g] - T2[g]);
}
`;

/**
 * Implicit diffusion: one Thomas sweep per (mode, component). Sequential in r,
 * so parallelism is only 2·gna wide — but the whole solve is a few hundred
 * thousand flops, far below the cost of a single advection pass.
 *
 * The isothermal rows are written directly in *mode* space: a row constant in φ
 * has DFT (value, 0, 0, …), so setting k = 0 and zeroing the rest reproduces
 * `applyBC` exactly once the inverse transform runs, saving a pass.
 *
 * **Free-slip side walls are imposed here too**, and for the same reason: a
 * walled box is the mirrored domain of `geometry.ts`, "T even about x = 0" is
 * exactly "T̂ real", and this kernel is already in mode space. The imaginary
 * half of the work simply writes zeros and returns — so the wall costs a
 * *negative* number of flops, needs no mirror-indexing anywhere else in the
 * pipeline, and is the identical projection `Temperature.diffuse` makes.
 */
export const tridiagSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> tri: array<f32>;
@group(0) @binding(2) var<storage, read> inRe: array<f32>;
@group(0) @binding(3) var<storage, read> inIm: array<f32>;
@group(0) @binding(4) var<storage, read_write> outRe: array<f32>;
@group(0) @binding(5) var<storage, read_write> outIm: array<f32>;
` + metric(geom) + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.gna * 2")}
  let k = g / 2; let imag = g % 2;
  let n = pp.gni;
  let last = (pp.gnr - 1) * pp.gna + k;

  // Boundary rows, in mode space (see header).
  if (imag == 0) {
    let d = select(0.0, 1.0, k == 0);
    outRe[k] = pp.tIn * d; outRe[last] = pp.tOut * d;
  } else {
    outIm[k] = 0.0; outIm[last] = 0.0;
  }
${geom.walls === "free-slip" ? /* wgsl */ `
  // Free-slip walls: Im T̂ ≡ 0 is the mirror projection (see header).
  if (imag == 1) {
    for (var i = 0; i < n; i++) { outIm[(i + 1) * pp.gna + k] = 0.0; }
    return;
  }
` : ""}
  let f = min(k, pp.gna - k) * 3 * n;
  let r1 = pp.ri + pp.dr; let rn = pp.ri + f32(n) * pp.dr;
  // Dirichlet data is constant in φ, so it forces the k = 0 mode only. The
  // first-derivative term is the metric's, and absent in a box (DH = 0).
  var d0 = 0.0; var dn = 0.0;
  if (k == 0 && imag == 0) {
    d0 = pp.dt * (1.0 / (pp.dr * pp.dr) - DH * ihOf(r1) / (2.0 * pp.dr)) * pp.tIn;
    dn = pp.dt * (1.0 / (pp.dr * pp.dr) + DH * ihOf(rn) / (2.0 * pp.dr)) * pp.tOut;
  }

  for (var i = 0; i < n; i++) {
    let g0 = (i + 1) * pp.gna + k;
    var rhs = select(inRe[g0], inIm[g0], imag == 1);
    if (i == 0) { rhs += d0; }
    if (i == n - 1) { rhs += dn; }
    let prev = select(0.0, select(outRe[g0 - pp.gna], outIm[g0 - pp.gna], imag == 1), i > 0);
    let x = (rhs - tri[f + i] * prev) / tri[f + 2 * n + i];
    if (imag == 0) { outRe[g0] = x; } else { outIm[g0] = x; }
  }
  for (var i = n - 2; i >= 0; i--) {
    let g0 = (i + 1) * pp.gna + k;
    let up = tri[f + n + i];
    if (imag == 0) { outRe[g0] -= up * outRe[g0 + pp.gna]; }
    else { outIm[g0] -= up * outIm[g0 + pp.gna]; }
  }
}
`;

/**
 * max|ψ| over the spline coefficients, which sets the streamline contour
 * spacing. Reducing the *coefficients* rather than the field is not an
 * approximation to apologise for: B-splines are non-negative and form a
 * partition of unity, so the surface lies inside the convex hull of its control
 * values and `max|c|` bounds `max|ψ|` — a valid, and very cheap, scale.
 *
 * Doing this on the GPU is what keeps the contour density meaningful as Ra and
 * the flow change without the host ever seeing ψ.
 */
export const psiMaxSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> psi: array<f32>;
@group(0) @binding(2) var<storage, read_write> stat: array<f32>;

var<workgroup> sm: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  var m = 0.0;
  for (var i = i32(t); i < pp.nr * pp.na; i += 256) { m = max(m, abs(psi[i])); }
  sm[t] = m;
  for (var s = 128u; s > 0u; s >>= 1u) {
    workgroupBarrier();
    if (t < s) { sm[t] = max(sm[t], sm[t + s]); }
  }
  if (t == 0u) { stat[2] = sm[0]; }
}
`;

/**
 * Nusselt number at both boundaries, reduced in a single workgroup so the only
 * host-visible traffic is a handful of floats — read back asynchronously for
 * the HUD, never inside the frame's dependency chain.
 *
 * The normalisation is `Geometry.nuScale`, written out of the uniform block
 * rather than baked as a constant so that it cannot disagree with the domain the
 * rest of the kernels were built for: `r ln(r_o/r_i)/2π` at a radius, `1/L`
 * across a box. Both make pure conduction read Nu = 1 at each end.
 */
export const nusseltSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> T: array<f32>;
@group(0) @binding(2) var<storage, read_write> out: array<f32>;

const NT: u32 = 256u;
var<workgroup> si: array<f32, 256>;
var<workgroup> so: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  let c = vec4f(48.0, -36.0, 16.0, -3.0) / (12.0 * pp.dr);
  let c0 = -25.0 / (12.0 * pp.dr);
  var qi = 0.0; var qo = 0.0;
  for (var j = i32(t); j < pp.gna; j += i32(NT)) {
    let a = pp.gna;
    qi += -(c0 * T[j] + c[0] * T[a + j] + c[1] * T[2 * a + j]
            + c[2] * T[3 * a + j] + c[3] * T[4 * a + j]);
    let e = (pp.gnr - 1) * a + j;
    qo += c0 * T[e] + c[0] * T[e - a] + c[1] * T[e - 2 * a]
          + c[2] * T[e - 3 * a] + c[3] * T[e - 4 * a];
  }
  si[t] = qi; so[t] = qo;
  for (var s = NT / 2u; s > 0u; s >>= 1u) {
    workgroupBarrier();
    if (t < s) { si[t] += si[t + s]; so[t] += so[t + s]; }
  }
  if (t == 0u) {
${geom.kind === "annulus" ? /* wgsl */ `
    let k = log(pp.ro / pp.ri) / TAU * pp.dphi;
    out[0] = k * pp.ri * si[0];
    out[1] = k * pp.ro * so[0];` : /* wgsl */ `
    let k = pp.dphi / pp.aLen;
    out[0] = k * si[0];
    out[1] = k * so[0];`}
  }
}
`;

/**
 * √⟨|u|²⟩ over the domain — area-weighted by the metric, so the number reads
 * the same in both geometries (`Temperature.rmsVelocity` is the CPU twin).
 * Reduced over the same `gnr × gna` grid `nusselt` walks, trapezoidal on r and
 * exact on the uniform periodic φ — `dφ` cancels between the numerator and
 * the area and so is written nowhere, here or there.
 *
 * Not folded into the `nusselt` kernel despite sharing that grid and a
 * workgroup shape: that reduction sums two *boundary* rows, `gna` points each,
 * while this one walks all `gnr` of them, and giving both jobs to one
 * workgroup would size the smaller one for the larger one's work.
 *
 * `out[3]` is the `stat` buffer's fourth float — see `gpu/sim.ts`.
 */
export const rmsSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> knots: array<f32>;
@group(0) @binding(2) var<storage, read> psi: array<f32>;
@group(0) @binding(3) var<storage, read_write> out: array<f32>;
` + metric(geom) + BASIS + VELOCITY + `
const NR: u32 = 256u;
var<workgroup> sn: array<f32, 256>;
var<workgroup> sa: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  let n = pp.gnr * pp.gna;
  var num = 0.0; var area = 0.0;
  for (var g = i32(t); g < n; g += i32(NR)) {
    let i = g / pp.gna; let j = g % pp.gna;
    let r = pp.ri + f32(i) * pp.dr;
    var w = pp.dr * hOf(r);
    if (i == 0 || i == pp.gnr - 1) { w *= 0.5; }
    let v = velocity(r, f32(j) * pp.dphi);
    num += w * dot(v, v);
    area += w;
  }
  sn[t] = num; sa[t] = area;
  for (var s = NR / 2u; s > 0u; s >>= 1u) {
    workgroupBarrier();
    if (t < s) { sn[t] += sn[t + s]; sa[t] += sa[t + s]; }
  }
  if (t == 0u) { out[3] = sqrt(sn[0] / sa[0]); }
}
`;

/**
 * `max(|u_r|/dr, |u_φ|/(h(r)dφ))` over the interior rows of the T grid — the
 * advective CFL measure `Temperature.maxSpeed` sizes the CPU reference's step
 * with. Reduced here for the same reason `rmsSource` is: reading `u` back to
 * size dt on the host would be the readback the frame loop forbids.
 *
 * Dispatched inside `stokes()`, after ψ is written, so it measures the flow
 * the *next* step advects with — not the one that just ran. `main.ts` reads it
 * back through `pollStats` and turns it into a step through `adaptiveDt`, which
 * lags by design (see its header): a few frames stale is far inside the safety
 * margin a Courant number already carries.
 *
 * `out[4]` is the `stat` buffer's fifth float.
 */
export const cflSource = (geom: Geometry) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> knots: array<f32>;
@group(0) @binding(2) var<storage, read> psi: array<f32>;
@group(0) @binding(3) var<storage, read_write> out: array<f32>;
` + metric(geom) + BASIS + VELOCITY + `
const NC: u32 = 256u;
var<workgroup> sc: array<f32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid: vec3u) {
  let t = lid.x;
  let n = (pp.gnr - 2) * pp.gna;   // interior rows only — the boundary carries no u
  var m = 0.0;
  for (var g = i32(t); g < n; g += i32(NC)) {
    let i = g / pp.gna + 1; let j = g % pp.gna;
    let r = pp.ri + f32(i) * pp.dr;
    let v = velocity(r, f32(j) * pp.dphi);
    m = max(m, max(abs(v.x) / pp.dr, abs(v.y) / (hOf(r) * pp.dphi)));
  }
  sc[t] = m;
  for (var s = NC / 2u; s > 0u; s >>= 1u) {
    workgroupBarrier();
    if (t < s) { sc[t] = max(sc[t], sc[t + s]); }
  }
  if (t == 0u) { out[4] = sc[0]; }
}
`;

/**
 * Screen position → the solver's (r, φ), and the half-extent of world the
 * viewport spans. The *only* part of the render pass that knows which domain it
 * is drawing: everything after it — the colour map, the contours, the mesh — is
 * written against (r, φ) and is identical in both.
 *
 * The annulus is polar, with the disk sized so that r_o covers `fill` of the
 * half-viewport. A box is an offset of the same square viewport, sized by its
 * *longer* side, so a long box fills the width and a unit one fills the frame;
 * depth increases upward on screen, which puts z = 0 — the hot boundary — at the
 * bottom, as the literature draws it.
 *
 * A **walled** box is drawn over `[0, width]`, which is half the period it is
 * solved on. The other half is its mirror image and carries no information, so
 * showing it would be showing the same cell twice. The ratio is emitted here
 * rather than passed as a uniform because the wall setting is a rebuild anyway.
 */
const domainFn = (g: Geometry): string => g.kind === "annulus"
  ? /* wgsl */ `
struct Dom { r: f32, phi: f32, inside: bool };

fn halfExtent() -> f32 { return pp.ro / pp.fill; }

fn domain(p: vec2f) -> Dom {
  let r = length(p);
  var phi = atan2(p.y, p.x);
  if (phi < 0.0) { phi += TAU; }
  return Dom(r, phi, r >= pp.ri && r <= pp.ro);
}
`
  : /* wgsl */ `
struct Dom { r: f32, phi: f32, inside: bool };

// Drawn width over solved period: 1 when x wraps, ½ when it is mirrored.
const SHOWN: f32 = ${g.walls === "free-slip" ? "0.5" : "1.0"};

fn width() -> f32 { return pp.aLen * SHOWN; }

fn halfExtent() -> f32 {
  return 0.5 * max(width(), pp.ro - pp.ri) / pp.fill;
}

fn domain(p: vec2f) -> Dom {
  let w = width();
  let x = p.x + 0.5 * w;
  let z = p.y + 0.5 * (pp.ri + pp.ro);
  return Dom(z, x, x >= 0.0 && x <= w && z >= pp.ri && z <= pp.ro);
}
`;

/**
 * (r, φ) → screen world position — the *inverse* of `domain` above, and the
 * particle vertex shader's only geometry-specific line: everywhere a pixel's
 * position tells `domain` where in the solver it is looking, this tells a
 * tracer's solver position where its dot belongs on screen. On a walled box
 * `width()` is already the drawn half of the solved period (see `domain`'s
 * own header), which is what keeps a tracer's dot appearing exactly where
 * `domain` would read the pixel under it back to the same (r, φ) — the
 * belt-and-braces fold `particles.ts`'s push already applies is what keeps a
 * tracer's φ inside that drawn half to begin with.
 */
const worldOfFn = (g: Geometry): string => g.kind === "annulus"
  ? /* wgsl */ `
fn worldOf(r: f32, phi: f32) -> vec2f { return r * vec2f(cos(phi), sin(phi)); }
`
  : /* wgsl */ `
fn worldOf(r: f32, phi: f32) -> vec2f {
  return vec2f(phi - 0.5 * width(), r - 0.5 * (pp.ri + pp.ro));
}
`;

/**
 * Render pass: screen → (r, φ) → the same monotone bicubic sample of T that the
 * solver uses, through a perceptually ordered map, overlaid with **ψ isocontours
 * — which are exactly the streamlines**, since `u = ∇×(ψ ẑ)` is tangent to level
 * sets of ψ. No particle tracing, no geometry, no second buffer: one spline
 * evaluation per pixel, and screen-space derivatives give the lines a constant
 * width at any zoom.
 *
 * The optional mesh overlay is drawn the same way — as a distance field, not as
 * geometry — from the element widths alone, since both discretisations are
 * uniform in (r, φ). It goes *under* the streamlines: the contours are the
 * reading, the mesh is what they are resolved on.
 *
 * Both fields are read straight from the solver's storage buffers, so a frame
 * never leaves the GPU.
 */
/** Colour-map control points as WGSL `vec3f` literals, in shader source order. */
const stopsWgsl = (colormap: ColormapName): string =>
  COLORMAPS[colormap].map(([r, g, b]) =>
    `vec3f(${r.toFixed(6)}, ${g.toFixed(6)}, ${b.toFixed(6)})`).join(", ");

export const renderSource = (geom: Geometry, colormap: ColormapName) => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> T: array<f32>;
@group(0) @binding(2) var<storage, read> knots: array<f32>;
@group(0) @binding(3) var<storage, read> psi: array<f32>;
@group(0) @binding(4) var<storage, read> stat: array<f32>;
` + metric(geom) + domainFn(geom) + CUBIC + CELL + sampleFn("T")
  + BASIS + PSI_AT + /* wgsl */ `
// One family of mesh lines: t is the distance to the nearest line and gap the
// spacing between neighbours, both in pixels.
//
// The gap term fades a family out as it approaches one line per pixel, the same
// policy the contours follow: past that the mesh is not a mesh any more, it is
// moiré, and it would read as texture in the temperature field.
fn meshLine(t: f32, gap: f32, w: f32) -> f32 {
  return (1.0 - smoothstep(0.0, w, t)) * smoothstep(1.5, 3.5, gap);
}

struct VSOut { @builtin(position) pos: vec4f, @location(0) p: vec2f };

@vertex fn vs(@builtin(vertex_index) i: u32) -> VSOut {
  var v = array(vec2f(-1, -1), vec2f(3, -1), vec2f(-1, 3));
  var o: VSOut;
  o.pos = vec4f(v[i], 0, 1);
  // fill = fraction of the half-viewport the domain's longest half-extent spans.
  // zoom/pan reparametrise which world window that half-viewport shows; every
  // downstream fwidth() is taken from the interpolated result, so contour and
  // mesh line widths stay a constant number of *screen* pixels at any zoom.
  o.p = v[i] * halfExtent() / pp.zoom + vec2f(pp.panX, pp.panY);
  return o;
}

@fragment fn fs(in: VSOut) -> @location(0) vec4f {
  let dom = domain(in.p);
  let r = dom.r; let phi = dom.phi;

  // Contour coordinate, evaluated for *every* pixel: fwidth is only defined in
  // uniform control flow, so it cannot sit behind the in-domain test below.
  // stat[2] = max|ψ| from the GPU reduction, so the spacing tracks the flow.
  let spacing = 2.0 * max(stat[2], 1e-20) / max(pp.levels, 1.0);
  let f = psiAt(clamp(r, pp.ri, pp.ro), phi) / spacing;
  let df = fwidth(f);
  // World units per pixel, for the mesh. Taken from the position rather than
  // from fwidth of the mesh coordinates themselves because φ has a branch cut on
  // the annulus: a screen-space derivative across it is enormous and
  // meaningless, and would erase the radial spokes along that one ray. p.x is
  // linear in screen x, so this is the exact scale, and it is the same for both
  // axes and both geometries.
  let px = max(fwidth(in.p.x), 1e-20);

  if (!dom.inside) { return vec4f(0.02, 0.02, 0.047, 1); }

  // Colour map control points (see colormaps.ts) — the shared table also
  // draws the pane's colour-bar legend, so the two can never disagree.
  var cm = array<vec3f, ${COLORMAPS[colormap].length}>(${stopsWgsl(colormap)});
  let u = clamp(sample_T(r, phi), 0.0, 1.0) * ${(COLORMAPS[colormap].length - 1).toFixed(1)};
  let i = min(${COLORMAPS[colormap].length - 2}, i32(u));
  var col = mix(cm[i], cm[i + 1], u - f32(i));

  // Element boundaries. Both discretisations are uniform in (r, φ) — clamped
  // uniform knots for ψ, a uniform grid for T — so a line family is just the
  // distance to the nearest multiple of its width, and the two meshes differ
  // only in how many elements they divide the domain into. Drawn at a lower
  // weight than the contours, and beneath them.
  if (pp.mesh > 0.0) {
    let spline = pp.mesh < 1.5;
    let hr = (pp.ro - pp.ri) / select(f32(pp.gnr - 1), f32(pp.nRel), spline);
    let ha = pp.aLen / select(f32(pp.gna), f32(pp.nAel), spline);
    // The transverse family's spacing on screen is an *arc*: on the annulus its
    // lines are spokes, so the gap grows with r and they fade from the inside
    // out rather than all at once. In a box h = 1 and the gap is uniform.
    let arc = ha * hOf(r);
    let w = 0.5 * pp.lineW;
    let a = 0.45 * max(
      meshLine(abs(fract((r - pp.ri) / hr - 0.5) - 0.5) * hr / px, hr / px, w),
      meshLine(abs(fract(phi / ha - 0.5) - 0.5) * arc / px, arc / px, w));
    let lum = dot(col, vec3f(0.30, 0.59, 0.11));
    col = mix(col, select(vec3f(1.0), vec3f(0.0), lum > 0.45), a);
  }

  if (pp.levels > 0.0) {
    let d = abs(fract(f - 0.5) - 0.5) / max(df, 1e-8);   // distance to a level, in pixels
    // Fade where the contours approach one per pixel: beyond Nyquist they would
    // alias into moiré, and an honest blank reads better than a false texture.
    let a = (1.0 - smoothstep(0.0, pp.lineW, d)) * 0.7 * (1.0 - smoothstep(0.3, 0.65, df));
    // Inferno spans black to near-white, so a fixed line colour disappears at
    // one end; pick against the local luminance instead.
    let lum = dot(col, vec3f(0.30, 0.59, 0.11));
    col = mix(col, select(vec3f(1.0), vec3f(0.0), lum > 0.45), a);
  }
  return vec4f(col, 1);
}
`;

// ---- particles: push, chemical projection, and their own render pass --------
//
// One tracer is 16 bytes, `r`/`phi`/`a`/`c` in that order — its position, its
// screen colour (always in [0, 1], whatever `particles.ts`'s tint registry
// currently fills it with), and its composition (`particles.ts` §5's
// two-species field). A push is a pointwise read-modify-write of one element:
// unlike the temperature grid, a tracer's next position depends on no
// neighbour's, so there is nothing here to double-buffer and nothing to race.

const PARTICLE = /* wgsl */ `
struct Particle { r: f32, phi: f32, a: f32, c: f32 };
`;

/**
 * Advance every tracer one step and, for the colour modes that track a
 * tracer's present surroundings rather than where it started, refresh `a`.
 *
 * `tint.wgsl` is `null` for a mode that needs nothing done here at all —
 * either the colour is fixed (*uniform*) or it was written once, at seeding,
 * on the host (*initial depth*, *initial φ*) — otherwise it is spliced in
 * verbatim as the new `a`, evaluated at the tracer's position *after* the
 * push, which is what makes *temperature* and *speed* read the flow a tracer
 * is in now rather than the one it left. Which extra buffer that expression
 * needs is read off its own text rather than carried as a separate flag,
 * so the registry entry stays the one place a colour mode is described.
 *
 * Binds `ps: Pass` — the *same* struct `advectSource` reads `dt` from, and at
 * runtime the same buffer: `passA` already carries exactly the dt the
 * temperature this step advects with, and rebinding it here is what keeps a
 * tracer's path and the field it is drawn over advancing through the same
 * instant of time, with no second dt uniform that could fall out of step
 * with the first (`particles.ts` §1 "dt").
 */
export const pushSource = (geom: Geometry, tint: TintEntry) => {
  const wantsT = (tint.wgsl ?? "").includes("sample_T");
  const wantsStat = (tint.wgsl ?? "").includes("stat[");
  const extraBind = wantsT
    ? "@group(0) @binding(6) var<storage, read> T: array<f32>;\n"
    : wantsStat
    ? "@group(0) @binding(6) var<storage, read> stat: array<f32>;\n"
    : "";
  const extraFns = wantsT ? CUBIC + CELL + sampleFn("T") : "";
  return PARAMS + PART + PARTICLE + /* wgsl */ `
struct Pass { dt: f32, limiter: i32, bc: i32, pad: i32 };
@group(0) @binding(1) var<uniform> pt: Part;
@group(0) @binding(2) var<uniform> ps: Pass;
@group(0) @binding(3) var<storage, read> knots: array<f32>;
@group(0) @binding(4) var<storage, read> psi: array<f32>;
@group(0) @binding(5) var<storage, read_write> Prt: array<Particle>;
${extraBind}` + metric(geom) + BASIS + VELOCITY + foldPhi(geom) + extraFns + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pt.count")}
  var p = Prt[g];
  let d = advanceRK2(p.r, p.phi, ps.dt);
  p.r = clamp(d.x, pp.ri, pp.ro);
  p.phi = foldPhi(d.y);
  let r = p.r; let phi = p.phi;
${tint.wgsl ? `  p.a = ${tint.wgsl};\n` : ""}  Prt[g] = p;
}
`;
};

/**
 * Zero the composition projection's two fixed-point accumulators
 * (`particles.ts`'s `CIC_FIXED_POINT_SCALE`) before each step's scatter —
 * `cnr·cna` threads, one accumulator pair each.
 *
 * No `PARAMS` here — this kernel touches nothing of the shared block, and a
 * binding declared but never read would be dropped from the `layout: "auto"`
 * bind group `sim.ts` derives from it, silently shifting every index after
 * it (see this file's own header on why every binding here is one its entry
 * point actually uses).
 */
export const clearSource = () => PART + /* wgsl */ `
@group(0) @binding(0) var<uniform> pt: Part;
@group(0) @binding(1) var<storage, read_write> num: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> den: array<atomic<u32>>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pt.cnr * pt.cna")}
  atomicStore(&num[g], 0u);
  atomicStore(&den[g], 0u);
}
`;

/**
 * Deposit every tracer's composition onto the four composition-grid nodes
 * bracketing it — cloud-in-cell, the same bilinear weights `sample_C` reads
 * back with — one thread per tracer, four `atomicAdd` pairs. WGSL has integer
 * atomics only, so each weighted contribution is scaled by
 * `CIC_FIXED_POINT_SCALE` and rounded before it is added; `particles.ts` has
 * the overflow argument for why that scale leaves headroom rather than a
 * limit meant to be approached.
 */
export const scatterSource = () => PARAMS + PART + PARTICLE + /* wgsl */ `
@group(0) @binding(1) var<uniform> pt: Part;
@group(0) @binding(2) var<storage, read> Prt: array<Particle>;
@group(0) @binding(3) var<storage, read_write> num: array<atomic<u32>>;
@group(0) @binding(4) var<storage, read_write> den: array<atomic<u32>>;
` + cellFn(C_GRID) + `
const SCALE: f32 = ${CIC_FIXED_POINT_SCALE.toFixed(1)};

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pt.count")}
  let p = Prt[g];
  let cell = cellC(p.r, p.phi);
  let na = pt.cna;
  let j1 = (cell.j + 1) % na;
  let w00 = (1.0 - cell.tr) * (1.0 - cell.tp);
  let w01 = (1.0 - cell.tr) * cell.tp;
  let w10 = cell.tr * (1.0 - cell.tp);
  let w11 = cell.tr * cell.tp;
  let i00 = cell.i * na + cell.j;
  let i01 = cell.i * na + j1;
  let i10 = (cell.i + 1) * na + cell.j;
  let i11 = (cell.i + 1) * na + j1;
  atomicAdd(&num[i00], u32(w00 * p.c * SCALE));
  atomicAdd(&den[i00], u32(w00 * SCALE));
  atomicAdd(&num[i01], u32(w01 * p.c * SCALE));
  atomicAdd(&den[i01], u32(w01 * SCALE));
  atomicAdd(&num[i10], u32(w10 * p.c * SCALE));
  atomicAdd(&den[i10], u32(w10 * SCALE));
  atomicAdd(&num[i11], u32(w11 * p.c * SCALE));
  atomicAdd(&den[i11], u32(w11 * SCALE));
}
`;

/**
 * num/den → C, `cnr·cna` threads, one composition-grid node each — the
 * "leave it alone" empty-cell policy `particles.ts` argues for: a node no
 * tracer reached this step keeps whatever composition it last held rather
 * than being reset to zero, since zero would read as "no dense material"
 * rather than "no data yet" and would spike the buoyancy load as tracers
 * wander through.
 *
 * **The free-slip mirror is imposed here too**, exactly as `tridiagSource`
 * imposes it for T: tracers only ever occupy `[0, width]` (`particles.ts`
 * §1–2), so a node in the drawn-nowhere upper half reads its accumulators
 * from its mirror partner in the lower half instead of its own — never from
 * that partner's *already-normalised* `C`, which a thread in a different
 * workgroup has no guarantee of having written yet within the same dispatch.
 * Reading the raw accumulators, which the prior `scatterSource` dispatch left
 * fully resolved before this one began, sidesteps that ordering question
 * rather than relying on one.
 *
 * No `PARAMS` here, for the reason `clearSource` has none: this kernel never
 * reads `pp`, and `layout: "auto"` drops a binding nothing reads.
 */
export const normaliseSource = (geom: Geometry) => PART + /* wgsl */ `
@group(0) @binding(0) var<uniform> pt: Part;
@group(0) @binding(1) var<storage, read> num: array<u32>;
@group(0) @binding(2) var<storage, read> den: array<u32>;
@group(0) @binding(3) var<storage, read_write> C: array<f32>;

@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pt.cnr * pt.cna")}
  let na = pt.cna;
  var idx = g;
${geom.walls === "free-slip" ? /* wgsl */ `
  let i = g / na; let j = g % na; let half = na / 2;
  if (j > half) { idx = i * na + (na - j); }
` : ""}
  let d = f32(den[idx]);
  if (d > 0.0) { C[g] = f32(num[idx]) / d; }
}
`;

/**
 * Draw every tracer as a soft-edged disc, `triangle-strip`/`draw(4, N)`
 * instanced over the particle buffer — the vertex shader indexes it directly
 * by `instance_index`, so there is no second copy of the positions and no
 * vertex layout to declare. Issued as a second pipeline inside the field's
 * own render pass (`gpu/sim.ts`'s `draw`), after it, so particles composite
 * over the field by draw order rather than needing a second attachment.
 *
 * The renderer is mode-agnostic: it only ever reads `a` and never asks which
 * colour mode filled it. `discrete` is the one exception, and it is a render
 * *pipeline* choice rather than a per-pixel one — the species map reads as
 * two bands, not a gradient, so `a` is thresholded at ½ instead of
 * interpolated across the colour map's control points.
 */
export const particleRenderSource = (
  geom: Geometry, colormap: ColormapName, discrete: boolean,
) => PARAMS + PART + PARTICLE + /* wgsl */ `
@group(0) @binding(1) var<uniform> pt: Part;
@group(0) @binding(2) var<storage, read> Prt: array<Particle>;
` + domainFn(geom) + worldOfFn(geom) + /* wgsl */ `
struct VSOut { @builtin(position) pos: vec4f, @location(0) uv: vec2f, @location(1) a: f32 };

@vertex fn vs(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VSOut {
  let p = Prt[ii];
  var corners = array(vec2f(-1.0, -1.0), vec2f(1.0, -1.0), vec2f(-1.0, 1.0), vec2f(1.0, 1.0));
  let uv = corners[vi];
  let world = worldOf(p.r, p.phi);
  let center = (world - vec2f(pp.panX, pp.panY)) * pp.zoom / halfExtent();
  var o: VSOut;
  o.pos = vec4f(center + uv * pt.radius * 2.0 / pt.side, 0.0, 1.0);
  o.uv = uv;
  o.a = p.a;
  return o;
}

@fragment fn fs(in: VSOut) -> @location(0) vec4f {
  let d2 = dot(in.uv, in.uv);
  if (d2 > 1.0) { discard; }
  var cm = array<vec3f, ${COLORMAPS[colormap].length}>(${stopsWgsl(colormap)});
  let a = clamp(in.a, 0.0, 1.0);
${discrete ? /* wgsl */ `
  let col = select(cm[0], cm[${COLORMAPS[colormap].length - 1}], a > 0.5);` : /* wgsl */ `
  let u = a * ${(COLORMAPS[colormap].length - 1).toFixed(1)};
  let i = min(${COLORMAPS[colormap].length - 2}, i32(u));
  let col = mix(cm[i], cm[i + 1], u - f32(i));`}
  // Soft edge: a disc, not a lit sphere — feathered so a dot survives at 2–3 px
  // without aliasing, not shaded as if it had depth.
  let edge = 1.0 - smoothstep(0.7, 1.0, sqrt(d2));
  return vec4f(col, pt.opacity * edge);
}
`;

// ---- the 3D cutaway-globe view: a purely cosmetic alternate render ----------
//
// See PLAN-3d-view.md. This is not a second simulation: the field is the same
// (r, φ) data `renderSource` draws, and the only new geometry is decorative
// (an Earth-textured exterior shell and a core sphere) around a wedge cut out
// of it that exposes two flat faces — each showing the *whole* field, exactly
// as `worldOf`'s annulus branch already lays it out, just embedded as a plane
// in 3D instead of filling the screen.
//
// Drawn by one fullscreen triangle, like `renderSource`, but the fragment
// shader casts a ray into a small analytic scene instead of mapping a pixel
// straight to (r, φ). Every intersection below is closed-form — a quadratic
// for each sphere, one division for a plane — so there is no marching loop,
// in the same spirit as the mesh overlay and the contours being distance
// fields rather than geometry. Contours and the mesh overlay are not drawn on
// the cut faces here — a documented simplification, not an oversight; see the
// plan.
//
// **The wedge is anchored at azimuth 0, on purpose.** `planeBasis(0)` is
// `(1, 0, 0)`, so the cut face at that azimuth already shares its in-plane
// axes with the flat 2D view's own (x, y) — camera azimuth/elevation 0 looks
// straight down −Z at exactly that plane with the same (right, up) the flat
// canvas uses. That is what lets `Globe3D`'s transition start at a camera
// pose that reproduces the 2D view pixel-for-pixel (its orthographic branch,
// `rayOrigin` below) rather than needing a second, unrelated "matching" pose
// to be found by hand.

/**
 * Camera and transition state, uniform like `Part` — nothing here is a
 * simulation quantity, so it has no business on the shared `Params` block.
 * `binding` is a parameter rather than fixed because the two pipelines that
 * use this struct disagree on what else shares group 0: the scene pass also
 * binds `Params`/`T` and so needs `gp` at binding 1, while the particle pass
 * needs nothing from the solver at all and puts it at binding 0.
 *
 * `persp` is the ortho→perspective blend `Globe3D`'s tween drives from 0 to
 * 1, `orthoHalf`/`panX`/`panY` are the 2D view's own `halfExtent`/zoom and
 * pan carried over so the ortho branch reproduces it exactly, and `reveal`
 * fades the exterior shell and core in without touching the cut faces —
 * which stay at full strength throughout, since they are the one thing that
 * must read as continuous with the 2D view a reader just came from.
 */
const globeStruct = (binding: number): string => /* wgsl */ `
struct Globe {
  eyeAz: f32, eyeEl: f32, eyeDist: f32, fovY: f32,
  wedgeW: f32, persp: f32, orthoHalf: f32, panX: f32,
  panY: f32, reveal: f32, gpad0: f32, gpad1: f32,
};
@group(0) @binding(${binding}) var<uniform> gp: Globe;
`;

/**
 * Cheap 3D value noise — a hash and its trilinear blend, then five fBm
 * octaves. Not physical, just enough irregularity that the exterior shell
 * reads as a planet rather than a lit sphere.
 */
const NOISE3 = /* wgsl */ `
fn hash3(p: vec3f) -> f32 {
  var p3 = fract(p * 0.3183099 + vec3f(0.10, 0.20, 0.30));
  p3 += dot(p3, p3.yzx + 19.19);
  return fract((p3.x + p3.y) * p3.z);
}

fn noise3(p: vec3f) -> f32 {
  let i = floor(p); let f = fract(p);
  let u = f * f * (3.0 - 2.0 * f);
  return mix(
    mix(mix(hash3(i + vec3f(0, 0, 0)), hash3(i + vec3f(1, 0, 0)), u.x),
        mix(hash3(i + vec3f(0, 1, 0)), hash3(i + vec3f(1, 1, 0)), u.x), u.y),
    mix(mix(hash3(i + vec3f(0, 0, 1)), hash3(i + vec3f(1, 0, 1)), u.x),
        mix(hash3(i + vec3f(0, 1, 1)), hash3(i + vec3f(1, 1, 1)), u.x), u.y),
    u.z);
}

fn fbm3(p: vec3f) -> f32 {
  var v = 0.0; var amp = 0.55; var freq = p;
  for (var o = 0; o < 5; o++) {
    v += amp * (noise3(freq) * 2.0 - 1.0);
    freq *= 2.02;
    amp *= 0.55;
  }
  return v;
}
`;

/**
 * The camera, as an orbit around the origin, and the two things every
 * fragment or vertex needs from it: a ray (the scene pass raymarches) or an
 * inverse projection (the particle pass rasterises). Shared by both sources
 * — `rayOrigin`/`rayDir` go unused by the particle pass and `projectNdc` by
 * the scene pass, but all four only ever read `gp`, which every caller binds
 * regardless.
 *
 * `rayOrigin`'s orthographic branch and `projectNdc`'s both solve the *same*
 * equation `Globe3D` writes into `orthoHalf`/`panX`/`panY` — one screen↔world
 * map, read in opposite directions, which is what keeps a tracer's dot
 * appearing exactly where the field pixel under it does.
 */
const CAMERA = /* wgsl */ `
struct Cam { eye: vec3f, fwd: vec3f, right: vec3f, up: vec3f };

fn camera() -> Cam {
  let ca = cos(gp.eyeAz); let sa = sin(gp.eyeAz);
  let ce = cos(gp.eyeEl); let se = sin(gp.eyeEl);
  let dir = vec3f(ce * sa, se, ce * ca);   // target → eye, unit
  let fwd = -dir;
  let worldUp = vec3f(0.0, 1.0, 0.0);
  let right = normalize(cross(fwd, worldUp));
  let up = cross(right, fwd);
  return Cam(dir * gp.eyeDist, fwd, right, up);
}

fn rayOrigin(ndc: vec2f, cam: Cam) -> vec3f {
  let ortho = (ndc.x * gp.orthoHalf + gp.panX) * cam.right
            + (ndc.y * gp.orthoHalf + gp.panY) * cam.up;
  return mix(ortho, cam.eye, gp.persp);
}

fn rayDir(ndc: vec2f, cam: Cam) -> vec3f {
  let t = tan(gp.fovY * 0.5);
  let persp = normalize(cam.fwd + ndc.x * t * cam.right + ndc.y * t * cam.up);
  return normalize(mix(cam.fwd, persp, gp.persp));
}

/** World point → NDC, the inverse of the two functions above — a rasterised vertex's own screen position rather than a fragment's ray. */
fn projectNdc(P: vec3f, cam: Cam) -> vec2f {
  let ortho = vec2f((dot(P, cam.right) - gp.panX) / gp.orthoHalf,
                     (dot(P, cam.up) - gp.panY) / gp.orthoHalf);
  let rel = P - cam.eye;
  let cz = max(dot(rel, cam.fwd), 1e-3);
  let t = tan(gp.fovY * 0.5);
  let persp = vec2f(dot(rel, cam.right) / (cz * t), dot(rel, cam.up) / (cz * t));
  return mix(ortho, persp, gp.persp);
}

/** Camera-space depth of a world point — monotone, not perspective-correct, and computed identically by the scene and particle passes so the shared depth buffer sorts them against each other correctly. */
const FAR: f32 = 32.0;
fn camDepth(P: vec3f, cam: Cam) -> f32 { return clamp(dot(P - cam.eye, cam.fwd) / FAR, 0.0, 1.0); }
`;

/**
 * The scene: an exterior shell of radius `ro`, a core of radius `ri`, and the
 * two flat faces the wedge exposes — every one of them a closed-form ray
 * intersection, and the *nearest valid* one wins. "Valid" for the shell is
 * *outside* the wedge's azimuth range (θ ∈ [0, wedgeW), the anchor the header
 * above explains); for a cut face it is *inside* the field's own radii,
 * ri ≤ ρ ≤ ro. The shell's far root is kept too, not just the near one — a
 * ray that enters through the wedge mouth and threads past the core without
 * meeting either cut face is looking at the *inside* of the shell on the far
 * side of the cavity, which reads better shaded as rock than as a hole into
 * nothing.
 *
 * No candidate depends on another's result, so the five below are independent
 * "compute, compare, keep" blocks — the same shape `meshLine`'s distance
 * fields already use, extended from two dimensions to three.
 */
const SCENE_TRACE = /* wgsl */ `
const KIND_BG: i32 = 0;
const KIND_OUTER: i32 = 1;
const KIND_CORE: i32 = 2;
const KIND_PLANE: i32 = 3;

/** The direction at plane azimuth θ — also that plane's in-plane x-axis; its in-plane y-axis is always world up. */
fn planeBasis(theta: f32) -> vec3f { return vec3f(cos(theta), 0.0, sin(theta)); }

struct Hit { t: f32, kind: i32, theta: f32 };

fn traceScene(O: vec3f, D: vec3f) -> Hit {
  var best = Hit(1e30, KIND_BG, 0.0);

  // exterior shell, radius ro — both roots (see header on why the far one matters too)
  {
    let b = dot(O, D); let c = dot(O, O) - pp.ro * pp.ro; let disc = b * b - c;
    if (disc > 0.0) {
      let s = sqrt(disc);
      for (var k = 0; k < 2; k++) {
        let t = select(-b - s, -b + s, k == 1);
        if (t > 1e-4 && t < best.t) {
          let P = O + t * D;
          let theta = atan2(P.z, P.x);
          if (theta < 0.0 || theta > gp.wedgeW) { best = Hit(t, KIND_OUTER, 0.0); }
        }
      }
    }
  }
  // core, radius ri — no angular restriction; the wedge only removes shell.
  {
    let b = dot(O, D); let c = dot(O, O) - pp.ri * pp.ri; let disc = b * b - c;
    if (disc > 0.0) {
      let s = sqrt(disc);
      var t = -b - s;
      if (t < 1e-4) { t = -b + s; }
      if (t > 1e-4 && t < best.t) { best = Hit(t, KIND_CORE, 0.0); }
    }
  }
  // the two cut faces, at θ = 0 and θ = wedgeW
  for (var k = 0; k < 2; k++) {
    let theta = select(0.0, gp.wedgeW, k == 1);
    let e1 = planeBasis(theta);
    let n = vec3f(-sin(theta), 0.0, cos(theta));
    let denom = dot(D, n);
    if (abs(denom) > 1e-6) {
      let t = -dot(O, n) / denom;
      if (t > 1e-4 && t < best.t) {
        let P = O + t * D;
        let rho = length(vec2f(dot(P, e1), P.y));
        if (rho >= pp.ri && rho <= pp.ro) { best = Hit(t, KIND_PLANE, theta); }
      }
    }
  }
  return best;
}

const LIGHT: vec3f = vec3f(0.5477, 0.6086, 0.5744);   // fixed in globe space, not the camera's

fn shadeOuter(P: vec3f) -> vec3f {
  let n = normalize(P);
  let land = smoothstep(-0.02, 0.05, fbm3(n * 2.2) - 0.02);
  let base = mix(vec3f(0.10, 0.28, 0.55), vec3f(0.20, 0.42, 0.16), land);
  let ice = smoothstep(0.72, 0.92, abs(n.y));
  let surf = mix(base, vec3f(0.92, 0.95, 0.97), ice);
  let diff = max(dot(n, LIGHT), 0.0);
  return surf * (0.28 + 0.72 * diff);
}

fn shadeCore(P: vec3f) -> vec3f {
  let t = clamp(fbm3(normalize(P) * 3.0 + vec3f(7.3, 1.1, 4.6)) + 0.35, 0.0, 1.0);
  return mix(vec3f(1.0, 0.74, 0.20), vec3f(1.0, 0.98, 0.88), t);
}
`;

/**
 * The scene pass: one fullscreen triangle, one ray per pixel, `frag_depth`
 * written explicitly so the particle pass below can be depth-tested against
 * it. `T`/`sample_T`/the colour map are exactly `renderSource`'s — see this
 * section's header on why the cut faces read the field rather than a copy
 * of it.
 */
export const globeSource = (colormap: ColormapName) => PARAMS + globeStruct(1) + /* wgsl */ `
@group(0) @binding(2) var<storage, read> T: array<f32>;
` + CUBIC + CELL + sampleFn("T") + NOISE3 + CAMERA + SCENE_TRACE + /* wgsl */ `
fn shadePlane(theta: f32, P: vec3f) -> vec3f {
  let u = dot(P, planeBasis(theta)); let v = P.y;
  var phi = atan2(v, u);
  if (phi < 0.0) { phi += TAU; }
  var cm = array<vec3f, ${COLORMAPS[colormap].length}>(${stopsWgsl(colormap)});
  let s = clamp(sample_T(length(vec2f(u, v)), phi), 0.0, 1.0) * ${(COLORMAPS[colormap].length - 1).toFixed(1)};
  let i = min(${COLORMAPS[colormap].length - 2}, i32(s));
  return mix(cm[i], cm[i + 1], s - f32(i));
}

struct VSOut { @builtin(position) pos: vec4f, @location(0) ndc: vec2f };

@vertex fn vs(@builtin(vertex_index) i: u32) -> VSOut {
  var v = array(vec2f(-1, -1), vec2f(3, -1), vec2f(-1, 3));
  var o: VSOut;
  o.pos = vec4f(v[i], 0, 1);
  o.ndc = v[i];
  return o;
}

struct FSOut { @location(0) col: vec4f, @builtin(frag_depth) depth: f32 };

@fragment fn fs(in: VSOut) -> FSOut {
  let cam = camera();
  let O = rayOrigin(in.ndc, cam);
  let D = rayDir(in.ndc, cam);
  let hit = traceScene(O, D);
  if (hit.kind == KIND_BG) { return FSOut(vec4f(0.02, 0.02, 0.047, 1.0), 1.0); }
  let P = O + hit.t * D;
  let depth = camDepth(P, cam);
  if (hit.kind == KIND_PLANE) { return FSOut(vec4f(shadePlane(hit.theta, P), 1.0), depth); }
  let shaded = select(shadeCore(P), shadeOuter(P), hit.kind == KIND_OUTER);
  let col = mix(vec3f(0.02, 0.02, 0.047), shaded, gp.reveal);
  return FSOut(vec4f(col, 1.0), depth);
}
`;

/**
 * Particles in the globe scene: the same `Prt`/`Part` a `GpuParticles`
 * already owns, borrowed the way `Globe3D` borrows `T` from the host
 * simulation — no second cloud, no second push. Every tracer is placed on
 * the θ = 0 cut face (`worldOf`'s own map, `r·(cos φ, sin φ)`, padded to
 * `z = 0`) rather than duplicated onto both, which is what keeps the overlay
 * reading as one layer instead of two copies of it.
 *
 * Rasterised, not raymarched — `projectNdc` is `rayOrigin`/`rayDir`'s
 * inverse — but depth-tested against the scene pass's `frag_depth` all the
 * same, via the identical `camDepth`, so a tracer is correctly hidden behind
 * the exterior shell once the camera orbits round to the solid side.
 */
export const globeParticleSource = (colormap: ColormapName, discrete: boolean) =>
  globeStruct(0) + PART + PARTICLE + /* wgsl */ `
@group(0) @binding(1) var<uniform> pt: Part;
@group(0) @binding(2) var<storage, read> Prt: array<Particle>;
` + CAMERA + /* wgsl */ `
struct VSOut {
  @builtin(position) pos: vec4f, @location(0) uv: vec2f,
  @location(1) a: f32, @location(2) depth: f32,
};

@vertex fn vs(@builtin(vertex_index) vi: u32, @builtin(instance_index) ii: u32) -> VSOut {
  let p = Prt[ii];
  var corners = array(vec2f(-1.0, -1.0), vec2f(1.0, -1.0), vec2f(-1.0, 1.0), vec2f(1.0, 1.0));
  let uv = corners[vi];
  let P = vec3f(p.r * cos(p.phi), p.r * sin(p.phi), 0.0);
  let cam = camera();
  let center = projectNdc(P, cam);
  let radiusNdc = pt.radius * 2.0 / pt.side;
  var o: VSOut;
  o.pos = vec4f(center + uv * radiusNdc, 0.0, 1.0);
  o.uv = uv;
  o.a = p.a;
  o.depth = camDepth(P, cam);
  return o;
}

struct FSOut { @location(0) col: vec4f, @builtin(frag_depth) depth: f32 };

@fragment fn fs(in: VSOut) -> FSOut {
  let d2 = dot(in.uv, in.uv);
  if (d2 > 1.0) { discard; }
  var cm = array<vec3f, ${COLORMAPS[colormap].length}>(${stopsWgsl(colormap)});
  let a = clamp(in.a, 0.0, 1.0);
${discrete ? /* wgsl */ `
  let col = select(cm[0], cm[${COLORMAPS[colormap].length - 1}], a > 0.5);` : /* wgsl */ `
  let u = a * ${(COLORMAPS[colormap].length - 1).toFixed(1)};
  let i = min(${COLORMAPS[colormap].length - 2}, i32(u));
  let col = mix(cm[i], cm[i + 1], u - f32(i));`}
  let edge = 1.0 - smoothstep(0.7, 1.0, sqrt(d2));
  return FSOut(vec4f(col, pt.opacity * edge), in.depth);
}
`;
