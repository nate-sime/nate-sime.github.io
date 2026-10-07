// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { gaussLegendre } from "../src/quad";
import {
  basisOnElement, basisRow, evalField, firstDof, greville, splineSpace, tabulate,
} from "../src/spline";

const SPACES = [
  [1, 4, 0], [2, 5, 1], [2, 5, 0], [3, 6, 2], [3, 6, 1], [4, 3, 3], [4, 7, 1], [5, 4, 4], [6, 5, 2],
] as const;

const samples = (n = 97) => Array.from({ length: n }, (_, i) => (i + 0.37) / n);

describe("Gauss–Legendre", () => {
  it("integrates x^q exactly on [0, 1] for q ≤ 2n − 1", () => {
    for (let n = 1; n <= 12; n++) {
      const { x, w } = gaussLegendre(n);
      expect(w.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 14);
      for (let q = 0; q <= 2 * n - 1; q++) {
        let s = 0;
        x.forEach((xi, i) => (s += w[i] * xi ** q));
        expect(Math.abs(s - 1 / (q + 1))).toBeLessThan(1e-14);
      }
    }
  });
});

describe("spline space", () => {
  it("has n = p + 1 + (ne − 1)(p − k) functions", () => {
    for (const [p, ne, k] of SPACES) {
      const s = splineSpace(p, ne, k);
      expect(s.U.length).toBe(s.n + p + 1);
      expect(s.n).toBe(p + 1 + (ne - 1) * (p - k));
    }
  });

  it("rejects continuity the degree cannot carry", () => {
    expect(() => splineSpace(3, 4, 3)).toThrow();
    expect(() => splineSpace(3, 4, -1)).toThrow();
  });

  it("is a partition of unity, with derivatives summing to zero", () => {
    for (const [p, ne, k] of SPACES) {
      const s = splineSpace(p, ne, k);
      for (const x of [0, 1, ...samples()]) {
        expect(basisRow(s, x, 0).reduce((a, b) => a + b, 0)).toBeCloseTo(1, 13);
        for (let r = 1; r <= Math.min(p, 3); r++) {
          const sum = basisRow(s, x, r).reduce((a, b) => a + b, 0);
          expect(Math.abs(sum)).toBeLessThan(1e-9 * ne ** r);
        }
      }
    }
  });

  it("differentiates in agreement with central differences", () => {
    for (const [p, ne, k] of SPACES) {
      const s = splineSpace(p, ne, k);
      const c = Float64Array.from({ length: s.n }, (_, i) => Math.sin(1.3 * i + 0.4));
      const eps = 1e-6;
      for (const x of samples(23)) {
        // Keep the stencil inside one element: the derivative jumps at a knot when k is low.
        const e = Math.floor(x * ne);
        if (x - eps < e / ne || x + eps > (e + 1) / ne) continue;
        const f = evalField(s, c, x, 2);
        const fp = evalField(s, c, x + eps, 1), fm = evalField(s, c, x - eps, 1);
        expect(f[1]).toBeCloseTo((fp[0] - fm[0]) / (2 * eps), 5);
        if (p >= 2) expect(Math.abs(f[2] - (fp[1] - fm[1]) / (2 * eps))).toBeLessThan(1e-4 * ne ** 2);
      }
    }
  });

  it("reproduces every polynomial of degree ≤ p (L² projection is exact)", () => {
    for (const [p, ne, k] of SPACES) {
      const s = splineSpace(p, ne, k);
      const t = tabulate(s, gaussLegendre(p + 2), 0);
      const P1 = p + 1;
      for (let q = 0; q <= p; q++) {
        // Assemble the dense mass matrix and load, solve, compare pointwise.
        const n = s.n, M = Array.from({ length: n }, () => new Float64Array(n)), b = new Float64Array(n);
        for (let e = 0; e < ne; e++)
          for (let iq = 0; iq < t.nq; iq++) {
            const i = e * t.nq + iq, f0 = firstDof(s, e);
            for (let a = 0; a < P1; a++) {
              b[f0 + a] += t.w[i] * t.B[i * P1 + a] * t.x[i] ** q;
              for (let bb = 0; bb < P1; bb++) M[f0 + a][f0 + bb] += t.w[i] * t.B[i * P1 + a] * t.B[i * P1 + bb];
            }
          }
        const c = gaussSolve(M, b);
        for (const x of samples(31)) expect(evalField(s, c, x)[0]).toBeCloseTo(x ** q, 10);
      }
    }
  });

  it("has linear precision at the Greville abscissae", () => {
    for (const [p, ne, k] of SPACES) {
      const s = splineSpace(p, ne, k), g = greville(s);
      for (const x of samples(19)) expect(evalField(s, g, x)[0]).toBeCloseTo(x, 13);
    }
  });

  it("is exactly C^k across interior knots — and no smoother", () => {
    for (const [p, ne, k] of SPACES) {
      if (ne < 2) continue;
      const s = splineSpace(p, ne, k);
      const c = Float64Array.from({ length: s.n }, (_, i) => Math.cos(2.1 * i) + 0.3 * i);
      for (let e = 1; e < ne; e++) {
        const x = s.breaks[e];
        const L = sideValues(s, c, e - 1, x, k + 1), R = sideValues(s, c, e, x, k + 1);
        for (let r = 0; r <= k; r++) expect(Math.abs(L[r] - R[r])).toBeLessThan(1e-9 * ne ** r * (1 + Math.abs(L[r])));
        // A generic coefficient vector jumps in the (k+1)-th derivative.
        expect(Math.abs(L[k + 1] - R[k + 1])).toBeGreaterThan(1e-6);
      }
    }
  });
});

function sideValues(s: ReturnType<typeof splineSpace>, c: Float64Array, e: number, x: number, d: number): Float64Array {
  const N = basisOnElement(s, e, x, d), P1 = s.p + 1, f0 = firstDof(s, e);
  return Float64Array.from({ length: d + 1 }, (_, r) => {
    let v = 0;
    for (let j = 0; j < P1; j++) v += c[f0 + j] * N[r * P1 + j];
    return v;
  });
}

/** Dense Gaussian elimination with partial pivoting — test-only. */
function gaussSolve(A: Float64Array[], b: Float64Array): Float64Array {
  const n = b.length, M = A.map((r) => Float64Array.from(r)), x = Float64Array.from(b);
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    [M[c], M[piv]] = [M[piv], M[c]];
    [x[c], x[piv]] = [x[piv], x[c]];
    for (let r = c + 1; r < n; r++) {
      const f = M[r][c] / M[c][c];
      for (let j = c; j < n; j++) M[r][j] -= f * M[c][j];
      x[r] -= f * x[c];
    }
  }
  for (let r = n - 1; r >= 0; r--) {
    let s = x[r];
    for (let j = r + 1; j < n; j++) s -= M[r][j] * x[j];
    x[r] = s / M[r][r];
  }
  return x;
}
