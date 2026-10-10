// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * A Gaussian field on the plate a × 1, by a separable Karhunen–Loève expansion.
 *
 * The covariance is a product, C(x, x′) = C₁(|x − x′|) C₁(|y − y′|), with C₁ the
 * beam's one-dimensional kernel. The covariance operator is then the tensor
 * product of two one-dimensional ones, and its eigenpairs are products of
 * theirs:
 *
 *   λ_ij = λ_i^x λ_j^y,   φ_ij(x, y) = φ_i^x(x) φ_j^y(y),
 *
 * so the two-dimensional expansion costs two one-dimensional Nyström solves —
 * no N² × N² eigenproblem — and keeps the M largest products. For the squared
 * exponential the product *is* the isotropic kernel, e^{−r²/2ℓ²}; for the
 * exponential and Matérn kernels it is the standard separable stand-in (the
 * exponential product is the textbook one: Ghanem & Spanos, §2.3.3). Along
 * either axis its smoothness is the one-dimensional kernel's.
 *
 * Along x the plate is a long, not 1: the expansion there is the one on [0, 1]
 * with length ℓ/a, evaluated at x/a — the same kernel, rescaled.
 *
 * Grid evaluation. A level's quadrature points are a tensor grid X × Y, and a
 * sample there is G = A Ξ Bᵀ — A[x, i] = √λ_i φ_i(x), B[y, j] likewise, Ξ the
 * normals placed at the kept (i, j) — in O(NX·M + NX·NY·J) for J distinct
 * j-indices, not the O(NX·NY·M) a flat table would take (and the megabytes it
 * would hold per level). Output is flat [gy·NX + gx], the layout `plate.ts`
 * reads its coefficients in.
 */

import type { PointField } from "./field";
import { modesAt, termsOf, type KLData } from "./kl";

/** The (i, j) of the M largest products λ_i^x λ_j^y, descending. */
export function productTerms(klx: KLData, kly: KLData, M: number): { i: Int32Array; j: Int32Array; values: Float64Array } {
  const all: [number, number, number][] = [];
  for (let i = 0; i < termsOf(klx); i++)
    for (let j = 0; j < termsOf(kly); j++) all.push([i, j, klx.values[i] * kly.values[j]]);
  all.sort((a, b) => b[2] - a[2] || a[0] + a[1] - (b[0] + b[1]) || a[0] - b[0]);
  const kept = all.slice(0, M);
  return {
    i: Int32Array.from(kept, (t) => t[0]),
    j: Int32Array.from(kept, (t) => t[1]),
    values: Float64Array.from(kept, (t) => t[2]),
  };
}

export class FieldOnGrid implements PointField {
  readonly n: number;
  readonly M: number;
  readonly NX: number;
  readonly NY: number;
  /** √λ_i φ_i at the x points, [gx·Mx + i]; likewise in y. */
  private readonly A: Float64Array;
  private readonly B: Float64Array;
  private readonly Mx: number;
  private readonly My: number;
  readonly terms: { i: Int32Array; j: Int32Array; values: Float64Array };
  /** Distinct j of the kept terms, and for each the kept terms with that j. */
  private readonly js: number[];
  private readonly byJ: number[][];
  /** s_M = Σ_t λ_t φ_t², the truncated pointwise variance. */
  readonly s: Float64Array;

  /** `klx` is the expansion on [0, 1] for length ℓ/a; x is scaled by `aspect` before it reads it. */
  constructor(klx: KLData, kly: KLData, M: number, readonly xs: ArrayLike<number>, readonly ys: ArrayLike<number>, aspect = 1) {
    this.terms = productTerms(klx, kly, M);
    this.M = this.terms.values.length;
    this.NX = xs.length;
    this.NY = ys.length;
    this.n = this.NX * this.NY;
    this.Mx = termsOf(klx);
    this.My = termsOf(kly);
    this.A = modesAt(klx, Float64Array.from(xs, (x) => x / aspect));
    this.B = modesAt(kly, ys);
    const groups = new Map<number, number[]>();
    for (let t = 0; t < this.M; t++) {
      const j = this.terms.j[t];
      if (!groups.has(j)) groups.set(j, []);
      groups.get(j)!.push(t);
    }
    this.js = [...groups.keys()];
    this.byJ = [...groups.values()];
    // Σ_t (A[x, i_t] B[y, j_t])²: the same contraction with every factor squared.
    this.s = this.contract(new Float64Array(this.M).fill(1), true);
  }

  /** Σ_t A[x, i_t] B[y, j_t] ξ_t, optionally with A and B squared. */
  private contract(xi: ArrayLike<number>, squared: boolean, out = new Float64Array(this.n)): Float64Array {
    const { NX, NY, A, B, Mx, My, terms } = this;
    const sq = (v: number) => (squared ? v * v : v);
    out.fill(0);
    const c = new Float64Array(NX);
    this.js.forEach((j, g) => {
      c.fill(0);
      for (const t of this.byJ[g]) {
        const i = terms.i[t], x = xi[t];
        for (let gx = 0; gx < NX; gx++) c[gx] += sq(A[gx * Mx + i]) * x;
      }
      for (let gy = 0; gy < NY; gy++) {
        const b = sq(B[gy * My + j]), o = gy * NX;
        if (b === 0) continue;
        for (let gx = 0; gx < NX; gx++) out[o + gx] += b * c[gx];
      }
    });
    return out;
  }

  gaussian(xi: ArrayLike<number>, out = new Float64Array(this.n)): Float64Array {
    return this.contract(xi, false, out);
  }

  lognormal(sigma: number, xi: ArrayLike<number>, out = new Float64Array(this.n)): Float64Array {
    this.gaussian(xi, out);
    for (let q = 0; q < this.n; q++) out[q] = Math.exp(sigma * out[q] - 0.5 * sigma * sigma * this.s[q]);
    return out;
  }
}
