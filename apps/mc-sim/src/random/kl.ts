// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Karhunen–Loève expansions of a stationary Gaussian field on [0, 1].
 *
 * A zero-mean field with unit variance and correlation C(|x − y|) is
 *
 *   g(x, ω) = Σ_j √λ_j φ_j(x) ξ_j(ω),   ξ_j iid N(0, 1),
 *
 * where (λ_j, φ_j) are the eigenpairs of the covariance operator,
 * ∫₀¹ C(|x − y|) φ(y) dy = λ φ(x), L²-orthonormal. Truncating after M terms
 * keeps Σ_{j<M} λ_j of the total variance Σ λ_j = ∫₀¹ C(0) = 1.
 *
 * Kernels, with correlation length ℓ:
 *
 *   exponential        e^{−r/ℓ}                  paths continuous, nowhere differentiable
 *   Matérn ν = 3/2     (1 + √3 r/ℓ) e^{−√3 r/ℓ}   paths once differentiable
 *   squared exponential e^{−r²/2ℓ²}               paths analytic
 *
 * and the smoothness of the paths is what sets how fast λ_j decays (j⁻², j⁻⁴,
 * faster than any power) — and, later, how fast the discretisation error of a
 * sample can fall.
 *
 * Nyström. The operator is discretised by composite Gauss–Legendre (nodes x_i,
 * weights w_i): the symmetric matrix A_ik = √w_i C(|x_i − x_k|) √w_k has
 * eigenvectors y_j, and φ_j(x_i) = y_ij / √w_i. Off the nodes the eigenfunction
 * is the Nyström interpolant
 *
 *   φ_j(x) = λ_j⁻¹ Σ_i w_i C(|x − x_i|) φ_j(x_i),
 *
 * which is exactly the eigen-equation with its integral done by the same rule —
 * so a field sample can be evaluated at any point, and in particular at the
 * quadrature points of any level of the hierarchy, from the same ξ.
 */

import { symEig } from "../eig";
import { gaussLegendre } from "../quad";

export type Kernel = "exponential" | "matern32" | "gaussian";

export function correlation(kernel: Kernel, ell: number): (r: number) => number {
  if (kernel === "exponential") return (r) => Math.exp(-r / ell);
  if (kernel === "matern32") {
    const a = Math.sqrt(3) / ell;
    return (r) => (1 + a * r) * Math.exp(-a * r);
  }
  return (r) => Math.exp((-0.5 * r * r) / (ell * ell));
}

/** The serialisable part of an expansion — what is posted to a worker. */
export interface KLData {
  readonly kernel: Kernel;
  readonly ell: number;
  /** Nyström nodes and weights on [0, 1]. */
  readonly x: Float64Array;
  readonly w: Float64Array;
  /** Every eigenvalue the solve resolved, descending (also those past the truncation). */
  readonly spectrum: Float64Array;
  /** The retained eigenvalues, descending; M = values.length. */
  readonly values: Float64Array;
  /** φ_j at the nodes, flat [j·N + i], L²-orthonormal under the rule. */
  readonly vectors: Float64Array;
}

/** Most terms an expansion keeps; the pane offers up to this. */
export const MAX_TERMS = 64;

/**
 * Eigenvalues below this fraction of the largest are dropped: the Nyström
 * interpolant divides by λ_j, so a λ_j at round-off would hand back an
 * eigenfunction made of round-off. Smooth kernels reach it within a few dozen
 * terms; the effective M is then smaller than the M asked for, and says so.
 */
const FLOOR = 1e-12;

/**
 * The expansion, truncated after `terms`, on N = 96 × 4 Gauss nodes by default.
 * Against the exponential kernel's exact eigenvalues the j-th is good to about
 * (j/N)² — 3·10⁻² at j = 64, 7·10⁻⁴ at j = 10: the kink of e^{−r/ℓ} at r = 0 caps
 * the rule at second order whatever its degree. The smoother kernels do better.
 * The solve is dense, O(N³), ~0.2 s; it is done once per kernel and length.
 */
export function karhunenLoeve(kernel: Kernel, ell: number, terms: number, panels = 96, g = 4): KLData {
  if (!(ell > 0)) throw new Error(`correlation length must be positive, got ${ell}`);
  const rule = gaussLegendre(g), N = panels * g;
  const x = new Float64Array(N), w = new Float64Array(N);
  for (let e = 0; e < panels; e++)
    for (let q = 0; q < g; q++) {
      x[e * g + q] = (e + rule.x[q]) / panels;
      w[e * g + q] = rule.w[q] / panels;
    }
  const C = correlation(kernel, ell), sw = w.map(Math.sqrt);
  const A = Array.from({ length: N }, (_, i) => Float64Array.from({ length: N }, (_, k) => sw[i] * C(Math.abs(x[i] - x[k])) * sw[k]));
  const { values, vectors } = symEig(A, true);
  const order = Array.from({ length: N }, (_, j) => N - 1 - j); // descending
  const spectrum = Float64Array.from(order, (j) => values[j]);
  let M = 0;
  while (M < Math.min(terms, N) && spectrum[M] > FLOOR * spectrum[0]) M++;
  const V = new Float64Array(M * N);
  for (let j = 0; j < M; j++) {
    const y = vectors![order[j]];
    // A fixed sign — positive mean, or positive slope at 0 for an odd mode — so
    // the j-th mode, and with it what ξ_j does, does not flip between kernels.
    let s = 0;
    for (let i = 0; i < N; i++) s += sw[i] * y[i];
    if (Math.abs(s) < 1e-8) s = y[1] / sw[1] - y[0] / sw[0];
    const sign = s < 0 ? -1 : 1;
    for (let i = 0; i < N; i++) V[j * N + i] = (sign * y[i]) / sw[i];
  }
  return { kernel, ell, x, w, spectrum, values: spectrum.slice(0, M), vectors: V };
}

/** Number of retained terms. */
export const termsOf = (kl: KLData): number => kl.values.length;

/** Σ_{j<M} λ_j: the share of the field's variance the truncation keeps. */
export const captured = (kl: KLData): number => kl.values.reduce((s, v) => s + v, 0);

/**
 * The scaled modes √λ_j φ_j at arbitrary points, flat [q·M + j] — what one
 * level tabulates once and every sample on it reuses: g(x_q) = Σ_j Φ[q, j] ξ_j.
 */
export function modesAt(kl: KLData, points: ArrayLike<number>): Float64Array {
  const N = kl.x.length, M = termsOf(kl), C = correlation(kl.kernel, kl.ell);
  const out = new Float64Array(points.length * M), c = new Float64Array(N);
  for (let q = 0; q < points.length; q++) {
    const xq = points[q];
    for (let i = 0; i < N; i++) c[i] = kl.w[i] * C(Math.abs(xq - kl.x[i]));
    for (let j = 0; j < M; j++) {
      let s = 0;
      const o = j * N;
      for (let i = 0; i < N; i++) s += c[i] * kl.vectors[o + i];
      // √λ φ(x) = λ⁻¹ᐟ² Σ_i w_i C φ(x_i).
      out[q * M + j] = s / Math.sqrt(kl.values[j]);
    }
  }
  return out;
}

/**
 * The exact eigenvalues of the exponential kernel on [0, 1], the first `count`
 * (Ghanem & Spanos, Stochastic Finite Elements, §2.3.3). With c = 1/ℓ, a = ½
 * and θ = ωa, they are λ = 2c / (ω² + c²) at the roots of
 *
 *   ca cos θ − θ sin θ = 0  (even modes),   θ cos θ + ca sin θ = 0  (odd modes),
 *
 * which alternate: the n-th lies in (nπ/2, (n + 1)π/2), so bisection needs no guess.
 */
export function exponentialEigenvalues(ell: number, count: number): Float64Array {
  const c = 1 / ell, a = 0.5, ca = c * a;
  const even = (t: number) => ca * Math.cos(t) - t * Math.sin(t);
  const odd = (t: number) => t * Math.cos(t) + ca * Math.sin(t);
  return Float64Array.from({ length: count }, (_, n) => {
    const f = n % 2 === 0 ? even : odd;
    let lo = (n * Math.PI) / 2, hi = ((n + 1) * Math.PI) / 2;
    let flo = f(lo);
    for (let it = 0; it < 200 && hi - lo > 1e-15 * hi; it++) {
      const mid = 0.5 * (lo + hi), fm = f(mid);
      if (Math.sign(fm) === Math.sign(flo)) { lo = mid; flo = fm; } else hi = mid;
    }
    const om = (0.5 * (lo + hi)) / a;
    return (2 * c) / (om * om + c * c);
  });
}

/** The first M terms of an expansion (fewer if it resolved fewer). */
export function truncate(kl: KLData, M: number): KLData {
  const m = Math.max(0, Math.min(M, termsOf(kl))), N = kl.x.length;
  return { ...kl, values: kl.values.slice(0, m), vectors: kl.vectors.slice(0, m * N) };
}
