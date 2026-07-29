/**
 * Quadrature gather tables for the buoyancy load and for the variable-μ
 * operator.
 *
 * `buoyancyLoad` in `step.ts` assembles ℓ_{il} = Ra ∫∫ T B_i N_l' dr dφ by
 * *scattering* from element quadrature points into DOFs — the natural loop for a
 * CPU, and useless on a GPU, where concurrent scatters need atomics. The same
 * integral is re-expressed here as a per-DOF **gather**: each DOF stores the
 * quadrature points in its support together with the weight `w · N^(d)`. Because
 * the geometry never moves, those tables are static and computed once in f64.
 *
 * The integrand does not separate (T is a bicubic sample of the grid, and the
 * monotone clamp makes it nonlinear in T), but the *basis* does, so the assembly
 * factors into three cheap passes rather than one quadratic one:
 *
 *   Tq[qr][qa] = T(r_qr, φ_qa)                        one sample per quad point
 *   G [qr][l ] = Σ_t aVal[l][t] · Tq[qr][ aIdx[l][t] ]        collapse φ
 *   b [i ][l ] = Ra Σ_t rVal[i][t] · G [ rIdx[i][t] ][l]      collapse r
 *
 * A degree-3 DOF is supported on at most P+1 elements, so a table row is at most
 * (P+1)·NQ long; shorter rows are padded with zero weight, which costs a few
 * wasted lanes and buys a branch-free inner loop.
 *
 * `applyOperator` reuses exactly that shape for the variable-μ stiffness apply —
 * the only structural difference is that its quadrature-point field is computed
 * from ψ rather than read from a grid, so a *fourth* pass sits in front.
 */

import { ANNULUS, type Geometry } from "../geometry";
import { Axis, P } from "../spline";
import { mat } from "../linalg";
import { gauss } from "../quad";

export const NQ = 5;                 // Gauss points per element (see src/quad.ts)
export const SLOTS = (P + 1) * NQ;   // quadrature points in one DOF's support

/**
 * The scalar a basis function contributes at a quadrature point: `N` are its
 * derivatives (rows 0..2), `p` the local index, `x` the abscissa. Everything the
 * two bilinear forms need is one of the six below.
 */
export type Weight = (N: Float64Array[], p: number, x: number) => number;

export const W0: Weight = (N, p) => N[0][p];
export const W1: Weight = (N, p) => N[1][p];
export const W2: Weight = (N, p) => N[2][p];

// Radial factors of the dissipation form (see `operators.ts` for the derivation):
// ε_rr = ∂_φ A[ψ] pairs a(r) with N'(φ); ε_rφ = −½C[ψ] pairs g(r) with N and
// −b(r) with N''. All three carry the metric, so they are built per geometry
// rather than being constants — in a box they are N′, N″ and N.
export const wA = (g: Geometry): Weight => (N, p, r) => {
  const ih = 1 / g.h(r);
  return N[1][p] * ih - g.dh * N[0][p] * ih * ih;
};
export const wG = (g: Geometry): Weight => (N, p, r) =>
  N[2][p] - g.dh * N[1][p] / g.h(r);
export const wB = (g: Geometry): Weight => (N, p, r) =>
  N[0][p] / (g.h(r) * g.h(r));

export interface QuadTable {
  x: Float64Array;    // abscissae, element-major: q = e·NQ + g
  idx: Int32Array;    // [n · SLOTS] quadrature point index
  val: Float64Array;  // [n · SLOTS] w · f at that point
}

/**
 * Gather table for one weight on `ax`. `idx` depends only on the axis, so tables
 * built with different weights are index-compatible and the GPU uploads a single
 * copy of it per axis.
 */
export function quadTable(ax: Axis, f: Weight): QuadTable {
  const els = ax.elements();
  const x = new Float64Array(els.length * NQ);
  const idx = new Int32Array(ax.n * SLOTS);
  const val = new Float64Array(ax.n * SLOTS);
  const fill = new Int32Array(ax.n);
  els.forEach(([a, b], e) => {
    gauss(a, b).forEach(({ x: u, w }, g) => {
      const q = e * NQ + g;
      x[q] = u;
      const { span, N } = ax.ders(u, 2);
      for (let p = 0; p <= P; p++) {
        const i = ax.dof(span, p), t = i * SLOTS + fill[i]++;
        idx[t] = q;
        val[t] = w * f(N, p, u);
      }
    });
  });
  return { x, idx, val };
}

export interface EvalTable {
  dof: Int32Array;    // [nq · (P+1)] DOF index of each locally supported basis fn
  val: Float64Array;  // [nq · (P+1)] the weight there, *without* the quadrature w
}

/**
 * The transpose of `quadTable`: per quadrature point, the P+1 basis functions
 * that are nonzero there. This is the scatter direction, which is exactly right
 * for *evaluating* a field at quadrature points (a read, not an accumulate) and
 * costs one basis evaluation per abscissa rather than per (r, φ) pair.
 */
export function evalTable(ax: Axis, f: Weight): EvalTable {
  const els = ax.elements(), nq = els.length * NQ, w = P + 1;
  const dof = new Int32Array(nq * w), val = new Float64Array(nq * w);
  els.forEach(([a, b], e) => {
    gauss(a, b).forEach(({ x: u }, g) => {
      const q = e * NQ + g;
      const { span, N } = ax.ders(u, 2);
      for (let p = 0; p <= P; p++) {
        dof[q * w + p] = ax.dof(span, p);
        val[q * w + p] = f(N, p, u);
      }
    });
  });
  return { dof, val };
}

/**
 * The gather form of `buoyancyLoad`, given T pre-sampled at the tensor grid of
 * quadrature points (`Tq[qr · nAq + qa]`). This is the CPU twin of the GPU
 * kernels and exists so the tables themselves can be regression-tested in f64
 * against the reference scatter assembly.
 */
export function tabulatedLoad(
  rt: QuadTable, at: QuadTable, nr: number, na: number, Tq: Float64Array, Ra: number,
): Float64Array[] {
  const nRq = rt.x.length, nAq = at.x.length;
  const G = mat(nRq, na);
  for (let q = 0; q < nRq; q++)
    for (let l = 0; l < na; l++) {
      let s = 0;
      for (let t = 0; t < SLOTS; t++)
        s += at.val[l * SLOTS + t] * Tq[q * nAq + at.idx[l * SLOTS + t]];
      G[q][l] = s;
    }

  const b = mat(nr, na);
  for (let i = 0; i < nr; i++)
    for (let l = 0; l < na; l++) {
      let s = 0;
      for (let t = 0; t < SLOTS; t++)
        s += rt.val[i * SLOTS + t] * G[rt.idx[i * SLOTS + t]][l];
      b[i][l] = Ra * s;
    }
  return b;
}

// ---- variable-μ operator ----------------------------------------------------

/**
 * Everything `applyOperator` needs, built once from the geometry. The three
 * radial and three azimuthal gather tables share one index array per axis; the
 * eval tables are per-axis, so the basis is evaluated `nRq + nAq` times at init
 * rather than `nRq · nAq` times per apply.
 */
export interface OperatorTables {
  rx: Float64Array; ax: Float64Array;                     // abscissae
  rIdx: Int32Array; rA: Float64Array; rG: Float64Array; rB: Float64Array;
  aIdx: Int32Array; aN0: Float64Array; aN1: Float64Array; aN2: Float64Array;
  rDof: Int32Array; eA: Float64Array; eG: Float64Array; eB: Float64Array;
  aDof: Int32Array; eN0: Float64Array; eN1: Float64Array; eN2: Float64Array;
}

export function operatorTables(
  rAx: Axis, aAx: Axis, geom: Geometry = ANNULUS,
): OperatorTables {
  const A = quadTable(rAx, wA(geom)), G = quadTable(rAx, wG(geom));
  const B = quadTable(rAx, wB(geom));
  const n0 = quadTable(aAx, W0), n1 = quadTable(aAx, W1), n2 = quadTable(aAx, W2);
  const eA = evalTable(rAx, wA(geom)), eG = evalTable(rAx, wG(geom));
  const eB = evalTable(rAx, wB(geom));
  const e0 = evalTable(aAx, W0), e1 = evalTable(aAx, W1), e2 = evalTable(aAx, W2);
  return {
    rx: A.x, ax: n0.x,
    rIdx: A.idx, rA: A.val, rG: G.val, rB: B.val,
    aIdx: n0.idx, aN0: n0.val, aN1: n1.val, aN2: n2.val,
    rDof: eA.dof, eA: eA.val, eG: eG.val, eB: eB.val,
    aDof: e0.dof, eN0: e0.val, eN1: e1.val, eN2: e2.val,
  };
}

/**
 * `∂_φA[ψ]` and `C[ψ]` at every quadrature point, from the spline coefficients.
 *
 * This is the first pass of the operator apply and, unweighted by μ, it is also
 * exactly the strain-rate tensor — so the power law shares it rather than
 * evaluating the basis a second time. `A[ψ] = ψ_r/h − h′ψ/h²` and
 * `C[ψ] = ψ_rr − (h′/h)ψ_r − ψ_φφ/h²` as in `operators.ts` — the metric is
 * already inside the tables, so this pass is the same in both geometries.
 */
function deformation(
  t: OperatorTables, c: Float64Array[],
): [Float64Array, Float64Array] {
  const nRq = t.rx.length, nAq = t.ax.length, W = P + 1;
  const d1 = new Float64Array(nRq * nAq), d2 = new Float64Array(nRq * nAq);
  for (let qr = 0; qr < nRq; qr++)
    for (let qa = 0; qa < nAq; qa++) {
      let s1 = 0, s2 = 0;
      for (let p = 0; p < W; p++) {
        const e = qr * W + p, row = c[t.rDof[e]];
        const a = t.eA[e], g = t.eG[e], b = t.eB[e];
        for (let q = 0; q < W; q++) {
          const f = qa * W + q, v = row[t.aDof[f]];
          s1 += a * t.eN1[f] * v;
          s2 += (g * t.eN0[f] - b * t.eN2[f]) * v;
        }
      }
      const o = qr * nAq + qa;
      d1[o] = s1; d2[o] = s2;
    }
  return [d1, d2];
}

/**
 * ε̇_II = √(ε_rr² + ε_rφ²) at the quadrature points, with `ε_rr = ∂_φA[ψ]`
 * and `ε_rφ = −½C[ψ]` — the deformation pass and one `hypot`. The power law
 * therefore costs the CPU one extra sweep of the basis and the GPU one extra
 * dispatch; no new tables, and nothing the operator did not already need.
 */
export function strainRate(t: OperatorTables, c: Float64Array[]): Float64Array {
  const [d1, d2] = deformation(t, c);
  return d1.map((v, i) => Math.hypot(v, 0.5 * d2[i]));
}

/**
 * Matrix-free application of the variable-μ dissipation form
 *
 *   a(ψ, v) = ∫ μ [ 4 ∂_φA[ψ] ∂_φA[v] + C[ψ] C[v] ] h dr dφ,
 *
 * with `A[ψ] = ψ_r/h − h′ψ/h²` and `C[ψ] = ψ_rr − (h′/h)ψ_r − ψ_φφ/h²` as in
 * `operators.ts`. `mu` is μ sampled at the tensor grid of quadrature points, so
 * it may vary in φ — which is precisely why the operator no longer separates and
 * the DFT solve of tier 1 no longer applies. Nothing else changes: the form, and
 * therefore the free-slip natural condition, is the same one, with
 * μ moved inside the integral.
 *
 * Three passes, mirroring the load assembly:
 *
 *   1. per quadrature point, the two stress components of ψ, weighted by μ r;
 *   2. collapse φ against N', N and N'';
 *   3. collapse r against a, g and b.
 *
 * The first and last radial DOFs are excluded on **both** sides — the input's
 * boundary rows never enter (the eval tables read them, but callers keep them
 * zero) and the output's are forced to zero. That is the discrete statement of
 * ψ = const, and doing it symmetrically is what keeps the operator SPD, which
 * conjugate gradients requires.
 */
export function applyOperator(
  t: OperatorTables, nr: number, na: number, c: Float64Array[], mu: Float64Array,
  geom: Geometry = ANNULUS,
): Float64Array[] {
  const nRq = t.rx.length, nAq = t.ax.length;

  // 1. Quadrature-point stresses: 4μh ∂_φA[ψ] and μh C[ψ].
  const [Pq, Qq] = deformation(t, c);
  for (let qr = 0; qr < nRq; qr++) {
    const h = geom.h(t.rx[qr]);
    for (let qa = 0; qa < nAq; qa++) {
      const o = qr * nAq + qa, m = mu[o] * h;
      Pq[o] *= 4 * m;
      Qq[o] *= m;
    }
  }

  // 2. Collapse the φ quadrature against N', N, N''.
  const GA = mat(nRq, na), GB = mat(nRq, na), GC = mat(nRq, na);
  for (let qr = 0; qr < nRq; qr++)
    for (let l = 0; l < na; l++) {
      let sa = 0, sb = 0, sc = 0;
      for (let s = 0; s < SLOTS; s++) {
        const j = l * SLOTS + s, q = qr * nAq + t.aIdx[j];
        sa += t.aN1[j] * Pq[q];
        sb += t.aN0[j] * Qq[q];
        sc += t.aN2[j] * Qq[q];
      }
      GA[qr][l] = sa; GB[qr][l] = sb; GC[qr][l] = sc;
    }

  // 3. Collapse the r quadrature; boundary rows stay zero (ψ = const, §3.1).
  const out = mat(nr, na);
  for (let i = 1; i < nr - 1; i++)
    for (let l = 0; l < na; l++) {
      let s = 0;
      for (let k = 0; k < SLOTS; k++) {
        const j = i * SLOTS + k, q = t.rIdx[j];
        s += t.rA[j] * GA[q][l] + t.rG[j] * GB[q][l] - t.rB[j] * GC[q][l];
      }
      out[i][l] = s;
    }
  return out;
}

/**
 * μ at the tensor grid of quadrature points, in the layout `applyOperator`
 * wants. `strain` supplies the second argument of the law — normally the
 * normalised array built from `strainRate` and `strainScale`, and absent (hence
 * 1, hence irrelevant at n = 1) for the μ(T, d) tier. The third is the radial
 * abscissa, which is what a depth-dependent law needs and what a law without one
 * simply does not declare.
 */
export function viscosityAt(
  t: OperatorTables, T: (r: number, phi: number) => number,
  mu: (T: number, strain: number, r: number) => number, strain?: Float64Array,
): Float64Array {
  const nRq = t.rx.length, nAq = t.ax.length;
  const out = new Float64Array(nRq * nAq);
  for (let qr = 0; qr < nRq; qr++)
    for (let qa = 0; qa < nAq; qa++) {
      const o = qr * nAq + qa;
      out[o] = mu(T(t.rx[qr], t.ax[qa]), strain ? strain[o] : 1, t.rx[qr]);
    }
  return out;
}
