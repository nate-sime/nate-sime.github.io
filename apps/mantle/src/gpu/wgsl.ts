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

/**
 * Uniform block shared by every kernel: 16 i32 then 24 f32, padded to 160 B — a
 * multiple of 16, which is what a uniform struct's size must round to. Layout is
 * mirrored by `I` and `F` in `gpu/sim.ts`, which is what the runtime controls
 * write through — `Ra`, `dt`, `levels`, `lineW`, `mesh`, `gamma`, `nExp`, `cz`
 * and the view transform (`zoom`, `panX`, `panY`) all change from the UI
 * without touching a pipeline.
 *
 * New fields are appended rather than slotted in beside a related one, so that
 * every existing index in `F` keeps its meaning: the JS side writes this block
 * by index, and a field inserted mid-struct would silently shift `mesh` and the
 * rest by one.
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

/** Catmull–Rom clamped to its bracketing values — `cubic` in temperature.ts. */
const CUBIC = /* wgsl */ `
fn cubic(p0: f32, p1: f32, p2: f32, p3: f32, t: f32) -> f32 {
  let v = p1 + 0.5 * t * (p2 - p0 + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3
    + t * (3.0 * (p1 - p2) + p3 - p0)));
  return clamp(v, min(p1, p2), max(p1, p2));
}
`;

/** Grid cell containing (r, φ): r clamped to the domain, φ wrapped. */
const CELL = /* wgsl */ `
struct Cell { i: i32, j: i32, tr: f32, tp: f32 };

fn cell(r: f32, phi: f32) -> Cell {
  let x = clamp((r - pp.ri) / pp.dr, 0.0, f32(pp.gnr - 1));
  let i = min(pp.gnr - 2, i32(floor(x)));
  var y = phi / pp.dphi;
  y = y - floor(y / f32(pp.gna)) * f32(pp.gna);
  let j = i32(floor(y));
  return Cell(i, j, x - f32(i), y - f32(j));
}
`;

/** Monotone bicubic sample of a grid buffer — `Temperature.sample`. */
const sampleFn = (buf: string) => /* wgsl */ `
fn sample_${buf}(r: f32, phi: f32) -> f32 {
  let c = cell(r, phi);
  let na = pp.gna;
  var col: array<f32, 4>;
  for (var m = -1; m <= 2; m++) {
    let row = clamp(c.i + m, 0, pp.gnr - 1) * na;
    col[m + 1] = cubic(${buf}[row + (c.j - 1 + na) % na], ${buf}[row + c.j % na],
                       ${buf}[row + (c.j + 1) % na], ${buf}[row + (c.j + 2) % na], c.tp);
  }
  return cubic(col[0], col[1], col[2], col[3], c.tr);
}
`;

/** Local extrema of the containing cell — the BFECC limiter bracket. */
const bracketFn = (buf: string) => /* wgsl */ `
fn bracket_${buf}(r: f32, phi: f32) -> vec2f {
  let c = cell(r, phi);
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

/** u = ∇×(ψ ẑ) evaluated from the spline coefficients — `Field.velocity`. */
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
`;

export const WG = 64; // workgroup size for the flat element-wise kernels

const flat = (n: string) => /* wgsl */ `
  let g = i32(gid.x);
  if (g >= ${n}) { return; }
`;

// ---- buoyancy load: three gather passes (see solver/assembly.ts) -------------

/** Tq[qr][qa] = T(r_qr, φ_qa): one bicubic sample per quadrature point. */
export const tqSource = () => PARAMS + /* wgsl */ `
@group(0) @binding(1) var<storage, read> T: array<f32>;
@group(0) @binding(2) var<storage, read> rq: array<f32>;
@group(0) @binding(3) var<storage, read> phiq: array<f32>;
@group(0) @binding(4) var<storage, read_write> Tq: array<f32>;
` + CUBIC + CELL + sampleFn("T") + `
@compute @workgroup_size(${WG})
fn main(@builtin(global_invocation_id) gid: vec3u) {
${flat("pp.nRq * pp.nAq")}
  Tq[g] = sample_T(rq[g / pp.nAq], phiq[g % pp.nAq]);
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
  // ψ = const is essential, so the boundary DOFs are not in the trial space.
  // Tier 1 simply ignores these rows; the Krylov tier would otherwise carry a
  // residual component its operator can never remove.
  if (i == 0 || i == pp.nr - 1) { b[g] = 0.0; return; }
  var s = 0.0;
  for (var t = 0; t < ${slots}; t++) {
    let q = i * ${slots} + t;
    s += rVal[q] * G[rIdx[q] * pp.na + l];
  }
  b[g] = pp.Ra * s;
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
 * Only `k = 0 … na/2` are stored (`A_k = A_{n−k}`); the two boundary DOFs are
 * zero (ψ = const).
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
  if (i == 0 || i == pp.nr - 1 || (k == 0 && pp.k0 == 0)) {
    outRe[g] = 0.0; outIm[g] = 0.0; return;
  }
  let base = min(k, pp.na - k) * pp.ni * pp.ni + (i - 1) * pp.ni;
  var sr = 0.0; var si = 0.0;
  for (var c = 0; c < pp.ni; c++) {
    let a = inv[base + c];
    sr += a * bRe[(c + 1) * pp.na + k];
    si += a * bIm[(c + 1) * pp.na + k];
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
  // symmetric, which conjugate gradients requires.
  if (i == 0 || i == pp.nr - 1) { out[g] = 0.0; return; }
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
