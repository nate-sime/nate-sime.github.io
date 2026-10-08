// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Closed-form answers for the uniform beam (e = μ = 1), against which every
 * discretisation error the app reports is measured.
 *
 * Frequencies. The modes of a uniform beam are combinations of sin, cos, sinh
 * and cosh of βx, and the supports reduce to one transcendental equation in β
 * (L = 1), with λ_n = ω̂_n² = β_n⁴. Each is written divided through by cosh β,
 * so it stays O(1) for every mode instead of overflowing past n ≈ 230:
 *
 *   clamped–free     cos β + sech β = 0
 *   clamped–clamped  cos β − sech β = 0
 *   pinned–pinned    sin β = 0                       (β_n = nπ)
 *   clamped–pinned   sin β − tanh β cos β = 0
 *
 * Every one changes sign exactly once on [nπ, (n + 1)π] for n ≥ 1 (and the
 * cantilever's first root lies in (0, π)), so bisection on those brackets finds
 * the n-th root without a starting guess.
 *
 * Statics. Under a uniform load q and point forces P_j at interior a_j,
 *
 *   w(x) = c₀ + c₁x + c₂x² + c₃x³ + q x⁴/24 + Σ_j P_j (x − a_j)₊³ / 6,
 *
 * and the four cᵢ follow from the two conditions at each end. Forces applied at
 * an end enter there as shear (they do no work on a support).
 */

import type { End, Supports } from "./beam";

type Pair = "CF" | "CC" | "PP" | "CP";

function pairOf({ left, right }: Supports): Pair | null {
  const key = [left, right].sort().join("-");
  const map: Record<string, Pair> = {
    "clamped-free": "CF", "clamped-clamped": "CC", "pinned-pinned": "PP", "clamped-pinned": "CP",
  };
  return map[key] ?? null;
}

const sech = (x: number) => 1 / Math.cosh(x);
const EQUATION: Record<Pair, (b: number) => number> = {
  CF: (b) => Math.cos(b) + sech(b),
  CC: (b) => Math.cos(b) - sech(b),
  PP: (b) => Math.sin(b),
  CP: (b) => Math.sin(b) - Math.tanh(b) * Math.cos(b),
};

/** β_n L for the n-th mode (n ≥ 1), or null where the supports have no closed form here. */
export function betaL(supports: Supports, n: number): number | null {
  const pair = pairOf(supports);
  if (!pair) return null;
  if (pair === "PP") return n * Math.PI;
  const f = EQUATION[pair];
  // The cantilever's n-th root is in the n-th bracket from 0; the others skip β = 0.
  const j = pair === "CF" ? n - 1 : n;
  let a = j * Math.PI + (j === 0 ? 1e-9 : 0), b = (j + 1) * Math.PI;
  let fa = f(a);
  for (let it = 0; it < 200 && b - a > 1e-15 * b; it++) {
    const c = 0.5 * (a + b), fc = f(c);
    if (fc === 0) return c;
    if (Math.sign(fc) === Math.sign(fa)) { a = c; fa = fc; } else b = c;
  }
  return 0.5 * (a + b);
}

/** The first `count` exact eigenvalues λ_n = (β_n L)⁴, or null. */
export function exactEigenvalues(supports: Supports, count: number): Float64Array | null {
  if (!pairOf(supports)) return null;
  return Float64Array.from({ length: count }, (_, i) => betaL(supports, i + 1)! ** 4);
}

export interface StaticLoad {
  /** Uniform distributed load. */
  readonly q: number;
  readonly forces: readonly { readonly x: number; readonly P: number }[];
}

export interface StaticSolution {
  /** r-th derivative of w at x, r ≤ 4 away from the a_j. */
  w(x: number, r?: number): number;
  /** ℓ(w) = ∫ q w + Σ P_j w(a_j): the work done by the load, twice the strain energy. */
  readonly compliance: number;
}

/** d^r/dx^r of x^i. */
const dpow = (i: number, r: number, x: number): number => {
  if (r > i) return 0;
  let c = 1;
  for (let k = 0; k < r; k++) c *= i - k;
  return c * x ** (i - r);
};

/** d^r/dx^r of (x − a)₊³ / 6. */
const dmac = (a: number, r: number, x: number): number => {
  if (x <= a || r > 3) return 0;
  const t = x - a;
  return [t ** 3 / 6, t ** 2 / 2, t, 1][r];
};

export function exactStatic(supports: Supports, load: StaticLoad): StaticSolution | null {
  if (!pairOf(supports)) return null;
  const interior = load.forces.filter((f) => f.x > 0 && f.x < 1);
  const endForce = (x: number) => load.forces.filter((f) => f.x === x).reduce((s, f) => s + f.P, 0);
  // Everything but the cubic: q x⁴/24 + Σ P (x − a)₊³/6.
  const part = (x: number, r: number) =>
    load.q * dpow(4, r, x) / 24 + interior.reduce((s, f) => s + f.P * dmac(f.x, r, x), 0);

  // Each row: Σ cᵢ (d^r x^i)(x₀) · sign = rhs − sign · part^(r)(x₀).
  const rows: number[][] = [], rhs: number[] = [];
  const cond = (x0: number, r: number, sign: number, value: number) => {
    rows.push([0, 1, 2, 3].map((i) => sign * dpow(i, r, x0)));
    rhs.push(value - sign * part(x0, r));
  };
  const at = (end: End, x0: number) => {
    const s0 = x0 === 0 ? 1 : -1; // w‴(0) = P₀, −w‴(1) = P₁ ; −w″(0) = 0, w″(1) = 0
    if (end === "clamped") { cond(x0, 0, 1, 0); cond(x0, 1, 1, 0); }
    if (end === "pinned") { cond(x0, 0, 1, 0); cond(x0, 2, 1, 0); }
    if (end === "free") { cond(x0, 2, 1, 0); cond(x0, 3, s0, endForce(x0)); }
  };
  at(supports.left, 0);
  at(supports.right, 1);
  const c = solve4(rows, rhs);

  const w = (x: number, r = 0) => c.reduce((s, ci, i) => s + ci * dpow(i, r, x), 0) + part(x, r);
  // ∫₀¹ w: the cubic, q x⁴/24 → q/120, and (x − a)₊³/6 → (1 − a)⁴/24.
  const integral = c.reduce((s, ci, i) => s + ci / (i + 1), 0) + load.q / 120 +
    interior.reduce((s, f) => s + (f.P * (1 - f.x) ** 4) / 24, 0);
  const compliance = load.q * integral + load.forces.reduce((s, f) => s + f.P * w(f.x), 0);
  return { w, compliance };
}

function solve4(A: number[][], b: number[]): number[] {
  const M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < 4; c++) {
    let piv = c;
    for (let r = c + 1; r < 4; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    for (let r = 0; r < 4; r++) {
      if (r === c) continue;
      const f = M[r][c] / M[c][c];
      for (let j = c; j <= 4; j++) M[r][j] -= f * M[c][j];
    }
  }
  return M.map((r, i) => r[4] / r[i]);
}

/**
 * The midspan amplitude |w(½)| of the uniform pinned–pinned beam under its
 * load applied harmonically at Ω, Rayleigh-damped (C = a M + b K) — by modal
 * superposition, which is exact here: the modes φ_n = √2 sin nπx are known and
 * Rayleigh damping does not couple them. With λ_n = (nπ)⁴,
 *
 *   w(½) = Σ_n φ_n(½) f_n / (λ_n − Ω² + iΩ(a + bλ_n)),
 *
 * f_n = ∫ φ_n = √2 (1 − cos nπ)/(nπ) for the uniform load and φ_n(½) for the
 * point load at midspan. Only odd n contribute; the terms fall as n⁻⁵ and n⁻⁴,
 * so 10⁵ of them leave a tail below 10⁻¹⁶ of the sum.
 */
export function pinnedResponse(load: "uniform" | "point", h: { Omega: number; a: number; b: number }): number {
  const W = h.Omega;
  let re = 0, im = 0;
  for (let n = 1; n <= 200001; n += 2) {
    const npi = n * Math.PI, lam = npi ** 4, s = n % 4 === 1 ? 1 : -1; // sin(nπ/2)
    const num = load === "uniform" ? (2 * s * 2) / npi : 2; // φ_n(½) f_n
    const dr = lam - W * W, di = W * (h.a + h.b * lam), d2 = dr * dr + di * di;
    re += (num * dr) / d2;
    im -= (num * di) / d2;
  }
  return Math.hypot(re, im);
}
