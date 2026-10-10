// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { Beam, SUPPORTS, admissible, type Load, type Supports } from "../src/beam/beam";
import { betaL, exactEigenvalues, exactStatic } from "../src/beam/exact";
import { gaussLegendre } from "../src/quad";

const ALL = Object.values(SUPPORTS) as Supports[];
const rate = (e: number[]) => Math.log2(e[e.length - 2] / e[e.length - 1]);

/** ‖(w_h − w)^(r)‖_L²: Gauss on every element at well above the spline degree. */
function errorNorm(beam: Beam, c: Float64Array, w: (x: number, r: number) => number, r: number): number {
  const { x: gx, w: gw } = gaussLegendre(beam.spec.p + 4), { breaks } = beam.space;
  let s = 0;
  for (let e = 0; e < beam.space.ne; e++) {
    const h = breaks[e + 1] - breaks[e];
    gx.forEach((xi, q) => {
      const x = breaks[e] + h * xi;
      s += h * gw[q] * (beam.evaluate(c, x, r)[r] - w(x, r)) ** 2;
    });
  }
  return Math.sqrt(s);
}

describe("exact uniform-beam frequencies", () => {
  // Blevins, "Formulas for Natural Frequency and Mode Shape", table 8-1.
  it("matches the tabulated β_n L", () => {
    const cases: [Supports, number[]][] = [
      [SUPPORTS.cantilever, [1.87510407, 4.69409113, 7.85475744, 10.99554073]],
      [SUPPORTS["clamped–clamped"], [4.73004074, 7.85320462, 10.99560784, 14.13716549]],
      [SUPPORTS["clamped–pinned"], [3.92660231, 7.06858275, 10.21017612, 13.35176878]],
      [SUPPORTS["pinned–pinned"], [Math.PI, 2 * Math.PI, 3 * Math.PI, 4 * Math.PI]],
    ];
    for (const [s, roots] of cases) roots.forEach((r, i) => expect(betaL(s, i + 1)).toBeCloseTo(r, 7));
    // Far up the spectrum the roots approach their asymptotes without overflow.
    expect(betaL(SUPPORTS.cantilever, 400)! - 399.5 * Math.PI).toBeCloseTo(0, 10);
  });
});

describe("admissibility", () => {
  it("refuses C⁰ spaces and mechanisms", () => {
    expect(admissible({ p: 3, ne: 8, k: 0, supports: SUPPORTS.cantilever })).toMatch(/C¹/);
    expect(admissible({ p: 1, ne: 8, k: 0, supports: SUPPORTS.cantilever })).toMatch(/C¹/);
    expect(admissible({ p: 3, ne: 8, supports: { left: "pinned", right: "free" } })).toMatch(/mechanism/);
    expect(admissible({ p: 3, ne: 8, k: 1, supports: SUPPORTS["clamped–clamped"] })).toBeNull();
  });
});

describe("statics", () => {
  it("is exact when the solution lies in the space", () => {
    // A uniform load gives a quartic: in every space of degree ≥ 4.
    for (const supports of ALL)
      for (const [p, k] of [[4, 3], [4, 1], [5, 2]]) {
        const beam = new Beam({ p, k, ne: 5, supports });
        const ex = exactStatic(supports, { q: 1, forces: [] })!;
        const c = beam.solve({ q: () => 1 });
        expect(errorNorm(beam, c, ex.w, 0)).toBeLessThan(1e-12);
        expect(errorNorm(beam, c, ex.w, 2)).toBeLessThan(1e-10);
      }
    // A point force at a knot gives a piecewise cubic, C² at the force: in the
    // cubic space whenever the force sits on a knot.
    for (const supports of ALL)
      for (const k of [1, 2]) {
        const tip = supports === SUPPORTS.cantilever;
        const forces = [{ x: tip ? 1 : 0.5, P: 1 }];
        const beam = new Beam({ p: 3, k, ne: 4, supports });
        const ex = exactStatic(supports, { q: 0, forces })!;
        const c = beam.solve({ forces });
        expect(errorNorm(beam, c, ex.w, 0)).toBeLessThan(1e-12);
        const compliance = beam.load({ forces }).reduce((s, f, i) => s + f * c[i], 0);
        expect(Math.abs(compliance / ex.compliance - 1)).toBeLessThan(1e-12);
      }
  });

  it("reproduces the textbook deflections", () => {
    const at = (s: Supports, x: number, load: Parameters<typeof exactStatic>[1]) => exactStatic(s, load)!.w(x);
    const q = { q: 1, forces: [] }, P = (x: number) => ({ q: 0, forces: [{ x, P: 1 }] });
    expect(at(SUPPORTS.cantilever, 1, q)).toBeCloseTo(1 / 8, 14);
    expect(at(SUPPORTS["pinned–pinned"], 0.5, q)).toBeCloseTo(5 / 384, 14);
    expect(at(SUPPORTS["clamped–clamped"], 0.5, q)).toBeCloseTo(1 / 384, 14);
    expect(at(SUPPORTS.cantilever, 1, P(1))).toBeCloseTo(1 / 3, 14);
    expect(at(SUPPORTS["pinned–pinned"], 0.5, P(0.5))).toBeCloseTo(1 / 48, 14);
    expect(at(SUPPORTS["clamped–clamped"], 0.5, P(0.5))).toBeCloseTo(1 / 192, 14);
    expect(at(SUPPORTS["clamped–pinned"], 0.5, P(0.5))).toBeCloseTo(7 / 768, 14);
  });

  it("converges at the a-priori rates (manufactured, variable stiffness, end loads)", () => {
    // Cantilever, e = eˣ, w = 1 − cos x: then (e w″)″ = −2eˣ sin x, and the free
    // end carries the moment e w″(1) and the force −(e w″)′(1) it implies.
    const wEx = (x: number, r: number) => [1 - Math.cos(x), Math.sin(x), Math.cos(x), -Math.sin(x)][r];
    const E = Math.E, c1 = Math.cos(1), s1 = Math.sin(1);
    const load: Load = {
      q: (x) => -2 * Math.exp(x) * Math.sin(x),
      moments: [{ x: 1, M: E * c1 }],
      forces: [{ x: 1, P: -E * (c1 - s1) }],
    };
    for (const p of [2, 3, 4])
      for (const k of [1, p - 1]) {
        // Coarser meshes at p = 4: by 32 elements the L² error has reached the
        // round-off floor (~10⁻¹⁰ here) that K's conditioning sets.
        const L2: number[] = [], H2: number[] = [];
        for (const ne of p < 4 ? [4, 8, 16, 32] : [4, 8, 16]) {
          const beam = new Beam({ p, k, ne, supports: SUPPORTS.cantilever, stiffness: Math.exp, nq: p + 3 });
          const c = beam.solve(load);
          L2.push(errorNorm(beam, c, wEx, 0));
          H2.push(errorNorm(beam, c, wEx, 2));
        }
        // ‖e″‖ ~ h^(p−1); ‖e‖ ~ h^(p+1), except that Aubin–Nitsche gains only
        // 2(p − 1) over it, which at p = 2 is the binding one.
        expect(rate(H2)).toBeCloseTo(p - 1, 0);
        expect(rate(L2)).toBeCloseTo(Math.min(p + 1, 2 * (p - 1)), 0);
      }
  });
});

describe("vibration", () => {
  it("bounds every exact eigenvalue from above, and converges at h^(2(p−1))", () => {
    for (const supports of ALL)
      for (const p of [2, 3, 4]) {
        const ex = exactEigenvalues(supports, 2)!;
        // p = 4 stops at 16 elements: by 32 the error (~10⁻¹¹) is at the floor
        // that solving with K (cond ~ h⁻⁴) leaves, and can land either side of λ.
        const err: number[] = [];
        for (const ne of p < 4 ? [4, 8, 16, 32] : [4, 8, 16]) {
          const { values } = new Beam({ p, ne, supports }).modes(2);
          // Rayleigh–Ritz in a conforming space can only overestimate.
          values.forEach((v, i) => expect(v).toBeGreaterThanOrEqual(ex[i] * (1 - 1e-13)));
          err.push((values[0] - ex[0]) / ex[0]);
        }
        expect(rate(err)).toBeCloseTo(2 * (p - 1), 0);
      }
  });

  it("returns M-orthonormal modes that agree with the full spectrum", () => {
    const beam = new Beam({ p: 3, k: 1, ne: 12, supports: SUPPORTS["clamped–pinned"] });
    const all = beam.spectrum(), { values, vectors } = beam.modes(4);
    for (let i = 0; i < 4; i++) {
      expect(values[i]).toBeCloseTo(all[i], 6);
      expect(Math.abs(values[i] - all[i]) / all[i]).toBeLessThan(1e-11);
      const Mv = beam.Mfull.matvec(vectors[i]);
      for (let j = 0; j <= i; j++) expect(vectors[j].reduce((s, v, r) => s + v * Mv[r], 0)).toBeCloseTo(i === j ? 1 : 0, 10);
    }
    // Every discrete eigenvalue lies above its exact counterpart.
    const ex = exactEigenvalues(beam.spec.supports, all.length)!;
    all.forEach((v, i) => expect(v).toBeGreaterThanOrEqual(ex[i] * (1 - 1e-11)));
  });
});
