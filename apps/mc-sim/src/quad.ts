// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Gauss–Legendre quadrature of any order, on the reference interval [0, 1].
 *
 * Nodes are the roots of the Legendre polynomial P_n, found by Newton's method
 * from the Tricomi estimate cos(π(i + ¾)/(n + ½)), which lands close enough that
 * a handful of iterations reach round-off for every n this app uses. An n-point
 * rule integrates polynomials of degree 2n − 1 exactly.
 */

export interface Rule {
  /** Nodes in [0, 1], ascending. */
  readonly x: Float64Array;
  /** Weights, summing to 1. */
  readonly w: Float64Array;
}

const cache = new Map<number, Rule>();

export function gaussLegendre(n: number): Rule {
  if (!Number.isInteger(n) || n < 1) throw new Error(`Gauss rule needs n ≥ 1, got ${n}`);
  const hit = cache.get(n);
  if (hit) return hit;
  const x = new Float64Array(n), w = new Float64Array(n);
  const m = (n + 1) >> 1;
  for (let i = 0; i < m; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let pp = 0;
    for (let it = 0; it < 100; it++) {
      // P_n(z) by the three-term recurrence, and P_n′ from P_n and P_{n−1}.
      let p1 = 1, p2 = 0;
      for (let j = 1; j <= n; j++) {
        const p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      pp = (n * (z * p1 - p2)) / (z * z - 1);
      const dz = p1 / pp;
      z -= dz;
      if (Math.abs(dz) < 1e-16) break;
    }
    // ξ = ∓z on [−1, 1], mapped to [0, 1]; the weight halves with the Jacobian.
    const wi = 1 / ((1 - z * z) * pp * pp);
    x[i] = (1 - z) / 2;
    x[n - 1 - i] = (1 + z) / 2;
    w[i] = w[n - 1 - i] = wi;
  }
  const rule = { x, w };
  cache.set(n, rule);
  return rule;
}
