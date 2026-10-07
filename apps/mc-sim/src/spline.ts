// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Univariate B-spline spaces of any degree and any inter-element continuity.
 *
 * A space is fixed by three numbers on a uniform mesh of `ne` elements:
 *
 *   p  degree,
 *   k  continuity across interior knots, C^k with 0 ≤ k ≤ p − 1,
 *
 * realised as an open (clamped) knot vector whose end knots repeat p + 1 times
 * and whose interior knots repeat m = p − k times. That one multiplicity is what
 * spans the three refinements the app teaches: h (more elements), p (degree up,
 * continuity held), and k (degree up *with* continuity, m = 1 — the smooth
 * spline spaces). p = 3, k = 1 is the classical Hermite cubic beam element;
 * p = 3, k = 2 is the cubic spline space on the same mesh.
 *
 *   n = p + 1 + (ne − 1)·m          basis functions,
 *   first[e] = e·m                  the first of the p + 1 nonzero on element e.
 *
 * Basis evaluation follows Piegl & Tiller, "The NURBS Book", A2.3, with the
 * degree a parameter rather than the constant the mantle app fixes it to.
 *
 * Layout. Everything a solver needs is precomputed by `tabulate` into flat
 * Float64Arrays indexed (element, quadrature point, derivative, local function)
 * — no per-point objects — so that the same tables can be uploaded unchanged as
 * GPU storage buffers when assembly moves to WebGPU.
 */

import type { Rule } from "./quad";

export interface SplineSpace {
  readonly p: number;
  readonly k: number;
  readonly ne: number;
  /** Interior knot multiplicity, p − k. */
  readonly m: number;
  /** Number of basis functions. */
  readonly n: number;
  readonly a: number;
  readonly b: number;
  readonly U: Float64Array;
  /** Element boundaries, ne + 1 of them. */
  readonly breaks: Float64Array;
}

export function splineSpace(p: number, ne: number, k = p - 1, a = 0, b = 1): SplineSpace {
  if (!Number.isInteger(p) || p < 1) throw new Error(`degree must be ≥ 1, got ${p}`);
  if (!Number.isInteger(ne) || ne < 1) throw new Error(`need ≥ 1 element, got ${ne}`);
  if (!Number.isInteger(k) || k < 0 || k > p - 1)
    throw new Error(`continuity C^${k} is not available at degree ${p} (0 ≤ k ≤ ${p - 1})`);
  const m = p - k;
  const breaks = Float64Array.from({ length: ne + 1 }, (_, i) => (i === ne ? b : a + ((b - a) * i) / ne));
  const U: number[] = [];
  for (let i = 0; i <= p; i++) U.push(a);
  for (let i = 1; i < ne; i++) for (let r = 0; r < m; r++) U.push(breaks[i]);
  for (let i = 0; i <= p; i++) U.push(b);
  return { p, k, ne, m, n: p + 1 + (ne - 1) * m, a, b, U: Float64Array.from(U), breaks };
}

/** Index of the first basis function nonzero on element `e`. */
export const firstDof = (s: SplineSpace, e: number): number => e * s.m;

/** The element containing x; a break belongs to the element on its right, except b. */
export function locate(s: SplineSpace, x: number): number {
  let e = Math.floor(((x - s.a) / (s.b - s.a)) * s.ne);
  e = Math.min(s.ne - 1, Math.max(0, e));
  // Floor of a rounded quotient can be one off at a break.
  if (e > 0 && x < s.breaks[e]) e--;
  else if (e < s.ne - 1 && x >= s.breaks[e + 1]) e++;
  return e;
}

/**
 * Values and derivatives up to order d of the p + 1 functions nonzero on
 * element `e`, at x. Row-major (d + 1) × (p + 1): entry [r·(p + 1) + j] is the
 * r-th derivative of function firstDof(e) + j. Derivatives above p are zero.
 *
 * x is not required to lie inside the element — evaluating element e's
 * polynomial pieces at its own end point is how one-sided limits at a knot are
 * read (`tests/spline.test.ts` uses that to measure continuity).
 */
export function basisOnElement(s: SplineSpace, e: number, x: number, d: number): Float64Array {
  const { p, U } = s;
  const span = p + e * s.m;
  const P1 = p + 1;
  const ndu = new Float64Array(P1 * P1); // ndu[j·P1 + r]
  const left = new Float64Array(P1), right = new Float64Array(P1);
  ndu[0] = 1;
  for (let j = 1; j <= p; j++) {
    left[j] = x - U[span + 1 - j];
    right[j] = U[span + j] - x;
    let saved = 0;
    for (let r = 0; r < j; r++) {
      ndu[j * P1 + r] = right[r + 1] + left[j - r];
      const tmp = ndu[r * P1 + j - 1] / ndu[j * P1 + r];
      ndu[r * P1 + j] = saved + right[r + 1] * tmp;
      saved = left[j - r] * tmp;
    }
    ndu[j * P1 + j] = saved;
  }
  const out = new Float64Array((d + 1) * P1);
  for (let j = 0; j <= p; j++) out[j] = ndu[j * P1 + p];
  const dd = Math.min(d, p);
  const a = [new Float64Array(P1), new Float64Array(P1)];
  for (let r = 0; r <= p; r++) {
    let s1 = 0, s2 = 1;
    a[0][0] = 1;
    for (let kk = 1; kk <= dd; kk++) {
      let der = 0;
      const rk = r - kk, pk = p - kk;
      if (r >= kk) {
        a[s2][0] = a[s1][0] / ndu[(pk + 1) * P1 + rk];
        der = a[s2][0] * ndu[rk * P1 + pk];
      }
      const j1 = rk >= -1 ? 1 : -rk;
      const j2 = r - 1 <= pk ? kk - 1 : p - r;
      for (let j = j1; j <= j2; j++) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[(pk + 1) * P1 + rk + j];
        der += a[s2][j] * ndu[(rk + j) * P1 + pk];
      }
      if (r <= pk) {
        a[s2][kk] = -a[s1][kk - 1] / ndu[(pk + 1) * P1 + r];
        der += a[s2][kk] * ndu[r * P1 + pk];
      }
      out[kk * P1 + r] = der;
      [s1, s2] = [s2, s1];
    }
  }
  let f = p;
  for (let kk = 1; kk <= dd; kk++) {
    for (let j = 0; j <= p; j++) out[kk * P1 + j] *= f;
    f *= p - kk;
  }
  return out;
}

/** Derivatives 0..d of the field Σ c_i B_i at x. */
export function evalField(s: SplineSpace, c: ArrayLike<number>, x: number, d = 0): Float64Array {
  const e = locate(s, x);
  const N = basisOnElement(s, e, x, d);
  const f0 = firstDof(s, e), P1 = s.p + 1;
  const out = new Float64Array(d + 1);
  for (let r = 0; r <= d; r++) {
    let v = 0;
    for (let j = 0; j < P1; j++) v += c[f0 + j] * N[r * P1 + j];
    out[r] = v;
  }
  return out;
}

/** Every basis function's r-th derivative at x, as a dense length-n vector. */
export function basisRow(s: SplineSpace, x: number, r = 0): Float64Array {
  const e = locate(s, x);
  const N = basisOnElement(s, e, x, r);
  const row = new Float64Array(s.n), f0 = firstDof(s, e), P1 = s.p + 1;
  for (let j = 0; j < P1; j++) row[f0 + j] = N[r * P1 + j];
  return row;
}

/**
 * Greville abscissae ξ_i = (U_{i+1} + … + U_{i+p}) / p. Plotting the control
 * coefficients c_i at ξ_i gives the control polygon, and Σ ξ_i B_i(x) = x.
 */
export function greville(s: SplineSpace): Float64Array {
  return Float64Array.from({ length: s.n }, (_, i) => {
    let g = 0;
    for (let r = 1; r <= s.p; r++) g += s.U[i + r];
    return g / s.p;
  });
}

/**
 * Basis tables at every quadrature point of every element.
 *
 *   x[e·nq + q]                       physical point
 *   w[e·nq + q]                       weight × Jacobian
 *   B[((e·nq + q)·(d + 1) + r)·(p + 1) + j]   r-th derivative of local function j
 */
export interface Tabulation {
  readonly space: SplineSpace;
  readonly nq: number;
  readonly d: number;
  readonly x: Float64Array;
  readonly w: Float64Array;
  readonly B: Float64Array;
}

export function tabulate(s: SplineSpace, rule: Rule, d: number): Tabulation {
  const nq = rule.x.length, P1 = s.p + 1, stride = (d + 1) * P1;
  const x = new Float64Array(s.ne * nq), w = new Float64Array(s.ne * nq);
  const B = new Float64Array(s.ne * nq * stride);
  for (let e = 0; e < s.ne; e++) {
    const x0 = s.breaks[e], h = s.breaks[e + 1] - x0;
    for (let q = 0; q < nq; q++) {
      const i = e * nq + q;
      x[i] = x0 + h * rule.x[q];
      w[i] = h * rule.w[q];
      B.set(basisOnElement(s, e, x[i], d), i * stride);
    }
  }
  return { space: s, nq, d, x, w, B };
}

/** A function sampled at a tabulation's quadrature points — the coefficient layout assembly reads. */
export function sampleAt(t: Tabulation, f: (x: number) => number): Float64Array {
  return Float64Array.from(t.x, f);
}
