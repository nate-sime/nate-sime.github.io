/**
 * Degree-p tensor-product B-spline space on (r, φ). CPU reference.
 *
 *   radial axis     : clamped (open) knot vector on [r_i, r_o] — or [0, 1] in a box,
 *   azimuthal axis  : uniform periodic knot vector on [0, span).
 *
 * A scalar stream function ψ = Σ_{ij} c_{ij} B_i(r) B_j(φ) yields, via
 * u = ∇×(ψ ẑ), the velocity
 *
 *   u_r = ψ_φ / h(r),     u_φ = −ψ_r,
 *
 * which is divergence-free pointwise and *exactly* — the property is structural
 * (div curl ≡ 0), independent of the coefficients and of the metric `h` that
 * `geometry.ts` supplies. The space is C^{p−1}; with p = 3 the velocity gradient
 * (hence the strain rate) is continuous.
 *
 * Basis evaluation follows Piegl & Tiller, "The NURBS Book", algorithms
 * A2.1 (knot span) and A2.3 (basis functions and derivatives).
 */

import { ANNULUS, type Geometry } from "./geometry";
import { mat, lu, solve } from "./linalg";

export const P = 3; // spline degree

function findSpan(U: Float64Array, nLast: number, u: number): number {
  if (u >= U[nLast + 1]) return nLast;
  if (u <= U[P]) return P;
  let lo = P, hi = nLast + 1, mid = (lo + hi) >> 1;
  while (u < U[mid] || u >= U[mid + 1]) {
    u < U[mid] ? (hi = mid) : (lo = mid);
    mid = (lo + hi) >> 1;
  }
  return mid;
}

/** Nonzero basis functions and derivatives up to order `d`: rows 0..d, cols 0..P. */
function basisDers(U: Float64Array, span: number, u: number, d: number): Float64Array[] {
  const ndu = mat(P + 1, P + 1), left = new Float64Array(P + 1), right = new Float64Array(P + 1);
  ndu[0][0] = 1;
  for (let j = 1; j <= P; j++) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;
    let saved = 0;
    for (let r = 0; r < j; r++) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const tmp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * tmp;
      saved = left[j - r] * tmp;
    }
    ndu[j][j] = saved;
  }
  const ders = mat(d + 1, P + 1);
  for (let j = 0; j <= P; j++) ders[0][j] = ndu[j][P];
  const a = mat(2, P + 1);
  for (let r = 0; r <= P; r++) {
    let s1 = 0, s2 = 1;
    a[0][0] = 1;
    for (let k = 1; k <= d; k++) {
      let der = 0;
      const rk = r - k, pk = P - k;
      if (r >= k) { a[s2][0] = a[s1][0] / ndu[pk + 1][rk]; der = a[s2][0] * ndu[rk][pk]; }
      const j1 = rk >= -1 ? 1 : -rk;
      const j2 = r - 1 <= pk ? k - 1 : P - r;
      for (let j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
        der += a[s2][j] * ndu[rk + j][pk];
      }
      if (r <= pk) { a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]; der += a[s2][k] * ndu[r][pk]; }
      ders[k][r] = der;
      [s1, s2] = [s2, s1];
    }
  }
  let f = P;
  for (let k = 1; k <= d; k++) { for (let j = 0; j <= P; j++) ders[k][j] *= f; f *= P - k; }
  return ders;
}

/** One parametric direction: knot vector, DOF indexing (mod n if periodic). */
export class Axis {
  readonly U: Float64Array;
  readonly nLast: number;
  constructor(readonly n: number, U: number[], readonly periodic: boolean) {
    this.U = Float64Array.from(U);
    this.nLast = this.U.length - P - 2;
  }
  dof(span: number, a: number): number {
    const i = span - P + a;
    return this.periodic ? ((i % this.n) + this.n) % this.n : i;
  }
  /**
   * On a periodic axis the parameter is wrapped first. `findSpan` clamps out-of-
   * range arguments, so without this a φ that strays outside [0, 2π) — as the
   * semi-Lagrangian trace routinely does near φ = 0 — would be *extrapolated*
   * off the end knot instead of read from the other side of the seam.
   */
  ders(u: number, d: number) {
    if (this.periodic) {
      const a = this.U[P], L = this.U[this.nLast + 1] - a;
      u = a + ((((u - a) % L) + L) % L);
    }
    const span = findSpan(this.U, this.nLast, u);
    return { span, N: basisDers(this.U, span, u, d) };
  }
  /** Non-degenerate knot intervals (elements) covering the domain. */
  elements(): [number, number][] {
    const out: [number, number][] = [];
    for (let i = P; i < this.U.length - P - 1; i++)
      if (this.U[i + 1] > this.U[i]) out.push([this.U[i], this.U[i + 1]]);
    return out;
  }
}

/** Clamped uniform axis with `n` basis functions on [a, b]. */
export function clampedAxis(n: number, a: number, b: number): Axis {
  const U: number[] = [];
  for (let i = 0; i <= P; i++) U.push(a);
  const nInt = n - P - 1;
  for (let i = 1; i <= nInt; i++) U.push(a + ((b - a) * i) / (nInt + 1));
  for (let i = 0; i <= P; i++) U.push(b);
  return new Axis(n, U, false);
}

/** Uniform periodic axis with `n` basis functions on [0, L). */
export function periodicAxis(n: number, L = 2 * Math.PI): Axis {
  const d = L / n, U: number[] = [];
  for (let k = 0; k <= n + 2 * P; k++) U.push((k - P) * d);
  return new Axis(n, U, true);
}

/** Interpolation sites: Greville abscissae (clamped) / knots (periodic). */
export function sites(ax: Axis): number[] {
  if (ax.periodic) {
    const d = ax.U[P + 1] - ax.U[P];
    return Array.from({ length: ax.n }, (_, b) => b * d);
  }
  return Array.from({ length: ax.n }, (_, i) => {
    let s = 0;
    for (let k = 1; k <= P; k++) s += ax.U[i + k];
    return s / P;
  });
}

function collocation(ax: Axis, s: number[]): Float64Array[] {
  const M = mat(ax.n, ax.n);
  s.forEach((u, row) => {
    const { span, N } = ax.ders(u, 0);
    for (let a = 0; a <= P; a++) M[row][ax.dof(span, a)] += N[0][a];
  });
  return M;
}

// ---- tensor-product stream function -----------------------------------------

export class Field {
  readonly c: Float64Array[]; // control coefficients [nr][nφ]
  /**
   * `geom` enters only through `h`, and defaults to the annulus because that is
   * the geometry the manufactured solutions and the convergence suite are
   * written on — see `src/mms.ts` and `src/verify.ts`.
   */
  constructor(readonly r: Axis, readonly a: Axis, readonly geom: Geometry = ANNULUS) {
    this.c = mat(r.n, a.n);
  }

  /** Interpolate g(r, φ) at the tensor grid of interpolation sites. */
  interpolate(g: (r: number, phi: number) => number): void {
    const rs = sites(this.r), as = sites(this.a);
    const Br = lu(collocation(this.r, rs)), Ba = lu(collocation(this.a, as));
    const nr = this.r.n, na = this.a.n;
    const X = mat(nr, na);
    for (let j = 0; j < na; j++) {
      const col = new Float64Array(nr);
      for (let i = 0; i < nr; i++) col[i] = g(rs[i], as[j]);
      const y = solve(Br, col);
      for (let i = 0; i < nr; i++) X[i][j] = y[i];
    }
    for (let i = 0; i < nr; i++) this.c[i] = solve(Ba, X[i]);
  }

  /** ψ and the derivatives needed for velocity, divergence and strain rate. */
  eval(r: number, phi: number) {
    const R = this.r.ders(r, 2), A = this.a.ders(phi, 2);
    let psi = 0, psi_r = 0, psi_p = 0, psi_rp = 0, psi_rr = 0, psi_pp = 0;
    for (let a = 0; a <= P; a++) {
      const i = this.r.dof(R.span, a);
      for (let b = 0; b <= P; b++) {
        const c = this.c[i][this.a.dof(A.span, b)];
        psi += R.N[0][a] * A.N[0][b] * c;
        psi_r += R.N[1][a] * A.N[0][b] * c;
        psi_p += R.N[0][a] * A.N[1][b] * c;
        psi_rp += R.N[1][a] * A.N[1][b] * c;
        psi_rr += R.N[2][a] * A.N[0][b] * c;
        psi_pp += R.N[0][a] * A.N[2][b] * c;
      }
    }
    return { psi, psi_r, psi_p, psi_rp, psi_rr, psi_pp };
  }

  /**
   * Shear strain rate ε_rφ = ½[ψ_φφ/h² − ψ_rr + (h′/h) ψ_r].
   * Vanishing at both boundaries is exactly the free-slip condition — the
   * `−u_φ h′/h` curvature term is what distinguishes this from `ω = 0`, and it
   * is the term a box does not have.
   */
  shearRate(r: number, phi: number): number {
    const { psi_r, psi_rr, psi_pp } = this.eval(r, phi);
    const ih = 1 / this.geom.h(r);
    return 0.5 * (psi_pp * ih * ih - psi_rr + this.geom.dh * psi_r * ih);
  }

  /**
   * Second invariant ε̇_II = √(ε_rr² + ε_rφ²), the strain-rate measure the power
   * law is a function of. Incompressibility gives ε_φφ = −ε_rr, so
   * √(½ ε:ε) collapses to these two components with no separate ε_φφ term.
   *
   * This is the definition; `strainRate` in `solver/assembly.ts` computes the
   * same quantity from the operator's own gather tables, and is checked against
   * this one.
   */
  strainRate(r: number, phi: number): number {
    const { psi_r, psi_p, psi_rp, psi_rr, psi_pp } = this.eval(r, phi);
    const ih = 1 / this.geom.h(r), dh = this.geom.dh;
    const e_rr = psi_rp * ih - dh * psi_p * ih * ih;
    return Math.hypot(e_rr, 0.5 * (psi_pp * ih * ih - psi_rr + dh * psi_r * ih));
  }

  velocity(r: number, phi: number) {
    const { psi_r, psi_p } = this.eval(r, phi);
    return { ur: psi_p / this.geom.h(r), up: -psi_r };
  }

  /** ∇·u assembled from independent velocity-gradient terms; ≈ 0 to round-off. */
  divergence(r: number, phi: number): number {
    const { psi_p, psi_rp } = this.eval(r, phi);
    const ih = 1 / this.geom.h(r), dh = this.geom.dh;
    const curvature = dh * psi_p * ih * ih;        // (h′/h) u_r
    const dur_dr = psi_rp * ih - dh * psi_p * ih * ih; // ∂u_r/∂r
    const duphi_dphi_over_h = -psi_rp * ih;        // (1/h) ∂u_φ/∂φ
    return curvature + dur_dr + duphi_dphi_over_h;
  }
}
