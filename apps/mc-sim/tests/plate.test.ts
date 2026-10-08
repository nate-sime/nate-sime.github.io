// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { denseGeneralizedEig, symEig } from "../src/eig";
import { FIELD_POINTS, Sampler, type McSpec } from "../src/mc/sampler";
import { Accumulator } from "../src/mc/stats";
import {
  REFERENCE, exactPlateEigenvalues, levyEigenvalues, navierEigenvalues, navierResponse, navierStatic,
} from "../src/plate/exact";
import { Plate, admissible, bandwidth, freeDofs } from "../src/plate/plate";
import {
  evaluatePlateQoI, plateHarmonicOf, plateLoadOf, plateQoiPoint, plateTheoryRate, runPlateHierarchy, type PlateCase,
} from "../src/plate/qoi";
import type { FieldSpec } from "../src/random/field";
import { FieldOnGrid } from "../src/random/field2d";
import { correlation, karhunenLoeve, type KLData } from "../src/random/kl";

const rate = (e: number[]) => Math.log2(e[e.length - 2] / e[e.length - 1]);
const square = (edges: PlateCase["edges"], load: PlateCase["load"] = "uniform"): PlateCase =>
  ({ edges, load, aspect: 1, nu: 0.3 });

describe("plate assembly", () => {
  it("numbers the free coefficients in a band of half-width p·nfx + p", () => {
    const pl = new Plate({ p: 3, ne: 4, edges: "CSFS" });
    // n = 7 per side; C removes two columns at x = 0, S one row at each of y = 0, 1.
    expect([pl.nfx, pl.nfy]).toEqual([5, 5]);
    expect(pl.dofs).toBe(freeDofs({ p: 3, ne: 4, edges: "CSFS" }));
    expect(pl.K.bw).toBe(bandwidth({ p: 3, ne: 4, edges: "CSFS" }));
    let reach = 0;
    pl.K.toDense().forEach((row, i) => row.forEach((v, j) => { if (v !== 0) reach = Math.max(reach, Math.abs(i - j)); }));
    expect(reach).toBe(pl.K.bw);
  });

  it("has exactly three rigid modes when free, and refuses mechanisms for statics", () => {
    const lam = symEig(new Plate({ p: 3, ne: 3, edges: "FFFF" }).K.toDense(), false).values;
    const top = lam[lam.length - 1];
    expect(Math.max(...lam.slice(0, 3).map(Math.abs)) / top).toBeLessThan(1e-13);
    expect(lam[3] / top).toBeGreaterThan(1e-6);
    expect(admissible({ p: 3, ne: 3, edges: "FFFF" })).toMatch(/mechanism/);
    expect(admissible({ p: 3, ne: 3, edges: "FFFF" }, true)).toBeNull();
    expect(admissible({ p: 3, ne: 3, edges: "SFFF" })).toMatch(/mechanism/);
    expect(admissible({ p: 3, ne: 3, edges: "SSFF" })).toBeNull();
    expect(admissible({ p: 3, ne: 3, edges: "CFFF" })).toBeNull();
    expect(admissible({ p: 3, ne: 3, k: 0, edges: "CCCC" })).toMatch(/C¹/);
  });

  it("reproduces a solution that lies in the space, for any ν", () => {
    // w = f(x) f(y), f = t²(1 − t)²: clamped on every edge, biquartic. Then
    // ∇⁴w = f⁗(x) f(y) + 2 f″(x) f″(y) + f(x) f⁗(y) is the load, and Galerkin
    // in a space containing w returns w itself.
    const f = (t: number) => t * t * (1 - t) ** 2, f2 = (t: number) => 2 - 12 * t + 12 * t * t;
    const q = (x: number, y: number) => 24 * f(y) + 2 * f2(x) * f2(y) + 24 * f(x);
    for (const nu of [0, 0.3])
      for (const [p, k] of [[4, 3], [4, 1], [5, 4]]) {
        const pl = new Plate({ p, k, ne: 3, edges: "CCCC", nu, stiffness: () => 1 });
        const c = pl.solve({ q });
        for (const [x, y] of [[0.5, 0.5], [0.2, 0.7], [0.91, 0.13]])
          expect(Math.abs(pl.evaluate(c, x, y) - f(x) * f(y)) / f(0.5) ** 2).toBeLessThan(1e-11);
      }
  });
});

describe("closed forms", () => {
  it("Navier: the textbook coefficients", () => {
    // Timoshenko & Woinowsky-Krieger, tables 8 and 23: 0.00406 q b⁴/D, 0.01160 P b²/D at the centre.
    expect(navierStatic(1, "uniform").w(0.5, 0.5)).toBeCloseTo(0.00406235, 7);
    expect(navierStatic(1, { x: 0.5, y: 0.5 }).w(0.5, 0.5)).toBeCloseTo(0.01160, 4);
    expect(navierEigenvalues(2, 3)).toEqual(Float64Array.from([1.25 ** 2, 2 ** 2, 3.25 ** 2], (v) => Math.PI ** 4 * v));
    // The grid evaluation is the pointwise one.
    const ns = navierStatic(1.5, "uniform"), g = ns.grid!([0.3, 0.75], [0.2, 0.5]);
    expect(g[3]).toBeCloseTo(ns.w(0.75, 0.5), 14);
    expect(g[0]).toBeCloseTo(ns.w(0.3, 0.2), 14);
  });

  it("Lévy: the CSCS frequencies, any aspect ratio, to the splines' last digits", () => {
    for (const aspect of [1, 2]) {
      const want = levyEigenvalues(aspect, 5);
      const got = new Plate({ p: 5, ne: 16, edges: "CSCS", aspect }).modes(5).values;
      want.forEach((w, i) => expect(Math.abs(got[i] / w - 1)).toBeLessThan(1e-7));
    }
    expect(Math.sqrt(levyEigenvalues(1, 1)[0])).toBeCloseTo(28.95085, 5);
  });

  it("Leissa: the splines converge to his CCCC and FFFF frequencies, and below his CFFF", () => {
    for (const edges of ["CCCC", "FFFF"] as const) {
      const ref = REFERENCE[edges]!.omega;
      const got = new Plate({ p: 4, ne: 16, edges }).modes(ref.length).values.map((v) => Math.sqrt(Math.max(v, 0)));
      ref.forEach((w, i) => (w === 0 ? expect(got[i]).toBeLessThan(1e-5) : expect(Math.abs(got[i] / w - 1)).toBeLessThan(6e-5)));
    }
    // CFFF from above, toward the converged values; Leissa's 1969 3.4917 sits 0.6% higher.
    const cfff = [8, 16, 32].map((ne) => Math.sqrt(new Plate({ p: 4, ne, edges: "CFFF" }).modes(1).values[0]));
    expect(cfff[0]).toBeGreaterThan(cfff[1]);
    expect(cfff[1]).toBeGreaterThan(cfff[2]);
    expect(cfff[2]).toBeCloseTo(REFERENCE.CFFF!.omega[0], 4);
    expect(3.4917 / cfff[2] - 1).toBeGreaterThan(0.005);
  });
});

describe("plate hierarchy", () => {
  it("converges at 2(p − 1) for ω₁ and the compliance of the simply supported plate", () => {
    for (const p of [2, 3, 4])
      for (const qoi of ["omega1", "compliance"] as const) {
        const H = runPlateHierarchy({ p, k: p - 1, plate: square("SSSS"), qoi, ne0: 2, levels: 4 });
        expect(H.alphaExact!).toBeCloseTo(plateTheoryRate(qoi, p, square("SSSS"))!, 0);
        expect(H.alpha!).toBeCloseTo(2 * (p - 1), 0);
      }
  });

  it("costs h⁻⁴ per banded solve: γ = 4, once ne is well past p", () => {
    const H = runPlateHierarchy({ p: 2, k: 1, plate: square("SSSS"), qoi: "compliance", ne0: 8, levels: 4 });
    // (ne + p − 2)² unknowns times a band of p(ne + p − 2): γ → 4 from below.
    expect(H.gamma!).toBeGreaterThan(3.8);
    expect(H.gamma!).toBeLessThan(4);
  });

  it("measures the field's error against the Navier series", () => {
    const H = runPlateHierarchy({ p: 3, k: 2, plate: square("SSSS"), qoi: "field", ne0: 2, levels: 4 });
    expect(H.alphaExact!).toBeCloseTo(4, 0);
    expect(H.exact!).toBeCloseTo(navierStatic(1, "uniform").rms, 12);
  });

  it("finds a clamped–free corner capping the cantilever plate's rate", () => {
    const H = runPlateHierarchy({ p: 4, k: 3, plate: square("CFFF"), qoi: "omega1", ne0: 2, levels: 5 });
    expect(H.exact).toBeNull();
    expect(H.alpha!).toBeLessThan(4); // against 2(p − 1) = 6 for a smooth mode
  });
});

describe("forced response of the plate", () => {
  it("is the modal superposition of the discrete modes, Rayleigh damping keeping them uncoupled", () => {
    const pc: PlateCase = { ...square("CSCS"), forcing: { ratio: 1.3, zeta: 0.03 } };
    const stiffness = (x: number, y: number) => 1 + 0.5 * x * y, mass = (x: number) => 1 + 0.3 * x;
    const pl = new Plate({ p: 3, ne: 4, edges: "CSCS", stiffness, mass });
    const h = plateHarmonicOf(pc), at = plateQoiPoint(pc), load = plateLoadOf(pc);
    const { Q } = evaluatePlateQoI(pl, "response", load, at, h);
    const F = pl.load(load), { values, vectors } = denseGeneralizedEig(pl.K.toDense(), pl.M.toDense());
    let re = 0, im = 0;
    values.forEach((lam, n) => {
      const phi = pl.expand(vectors[n]), pf = vectors[n].reduce((s, v, i) => s + v * F[i], 0);
      const v = pl.evaluate(phi, at.x, at.y) * pf;
      const dr = lam - h.Omega ** 2, di = h.Omega * (h.a + h.b * lam), d2 = dr * dr + di * di;
      re += (v * dr) / d2;
      im -= (v * di) / d2;
    });
    expect(Math.abs(Q / Math.hypot(re, im) - 1)).toBeLessThan(1e-10);
  });

  it("converges to the Navier modal series on the simply supported plate", () => {
    const pc = square("SSSS"), h = plateHarmonicOf(pc), at = plateQoiPoint(pc);
    const exact = navierResponse(1, "uniform", at, h);
    const err = [4, 8, 16].map((ne) => {
      const pl = new Plate({ p: 3, ne, edges: "SSSS" });
      return Math.abs(evaluatePlateQoI(pl, "response", plateLoadOf(pc), at, h).Q / exact - 1);
    });
    expect(err[2]).toBeLessThan(1e-5);
    expect(rate(err)).toBeGreaterThan(3.5);
    // Static limit: Ω → 0 gives back the Navier deflection.
    expect(navierResponse(1, "uniform", at, { Omega: 0, a: 0, b: 0 })).toBeCloseTo(navierStatic(1, "uniform").w(0.5, 0.5), 10);
  });
});

describe("the separable random field", () => {
  const kl = karhunenLoeve("matern32", 0.4, 64);

  it("reproduces the product covariance from its kept terms", () => {
    const xs = [0.1, 0.45, 0.8], ys = [0.2, 0.6], M = 600;
    const f = new FieldOnGrid(kl, kl, M, xs, ys);
    const G = Array.from({ length: f.M }, (_, t) => f.gaussian(Float64Array.from({ length: f.M }, (_, u) => (u === t ? 1 : 0))));
    const C = correlation("matern32", 0.4);
    for (const [a, b] of [[0, 4], [1, 5], [2, 3], [0, 0]]) {
      const cov = G.reduce((s, g) => s + g[a] * g[b], 0);
      const ax = a % 3, ay = Math.floor(a / 3), bx = b % 3, by = Math.floor(b / 3);
      // 600 of the products keep all but ~10⁻⁵ of the variance (the 1D spectrum falls as j⁻⁴).
      expect(cov).toBeCloseTo(C(Math.abs(xs[ax] - xs[bx])) * C(Math.abs(ys[ay] - ys[by])), 4);
    }
    // s is that covariance on the diagonal.
    expect(f.s[4]).toBeCloseTo(G.reduce((s, g) => s + g[4] ** 2, 0), 12);
  });

  it("is the same function of ω on every grid — the coupling MLMC needs", () => {
    const coarse = new FieldOnGrid(kl, kl, 40, [0, 0.5, 1], [0.25, 0.75]);
    const fine = new FieldOnGrid(kl, kl, 40, [0, 0.25, 0.5, 0.75, 1], [0.25, 0.5, 0.75]);
    const xi = Float64Array.from({ length: 40 }, (_, i) => Math.sin(3.7 * i + 1));
    const gc = coarse.gaussian(xi), gf = fine.gaussian(xi);
    expect(gf[0 * 5 + 2]).toBeCloseTo(gc[0 * 3 + 1], 13);
    expect(gf[2 * 5 + 4]).toBeCloseTo(gc[1 * 3 + 2], 13);
  });

  it("stretches along x with the aspect ratio", () => {
    // A 2 × 1 plate: along x the expansion for ℓ/2 at x/2 has the correlation of length ℓ.
    const C = correlation("gaussian", 0.5), klx = karhunenLoeve("gaussian", 0.25, 64), kly = karhunenLoeve("gaussian", 0.5, 64);
    const f = new FieldOnGrid(klx, kly, 400, [0.2, 1.1], [0.5], 2);
    const G = Array.from({ length: f.M }, (_, t) => f.gaussian(Float64Array.from({ length: f.M }, (_, u) => (u === t ? 1 : 0))));
    expect(G.reduce((s, g) => s + g[0] * g[1], 0)).toBeCloseTo(C(0.9), 6);
  });
});

describe("Monte Carlo on the plate", () => {
  // As montecarlo.test.ts: a perfectly correlated field scales K by one number per sample.
  const constantKL: KLData = {
    kernel: "gaussian", ell: 1e9, x: Float64Array.of(0.5), w: Float64Array.of(1),
    spectrum: Float64Array.of(1), values: Float64Array.of(1), vectors: Float64Array.of(1),
  };
  const field: FieldSpec = { kernel: "gaussian", ell: 1e9, sigma: 0.4, terms: 1, massFollows: false, loadSigma: 0 };
  const spec = (qoi: McSpec["qoi"]): McSpec => ({
    p: 3, k: 2, beam: { supports: "pinned–pinned", load: "uniform", section: "uniform" },
    plate: square("SSSS"), qoi, field, seed: 3,
  });

  it("E[w] = w₀ e^{σ²} and E[ω₁] = ω₁ₕ e^{−σ²/8}, on the discrete plate", () => {
    const sigma = 0.4, n = 3000;
    for (const qoi of ["deflection", "omega1"] as const) {
      const s = new Sampler(spec(qoi), constantKL), acc = new Accumulator(FIELD_POINTS);
      const Q = new Float64Array(n), Qc = new Float64Array(n).fill(NaN);
      for (let i = 0; i < n; i++) Q[i] = s.solve(i, 4).Q;
      acc.push({ from: 0, count: n, Q, Qc, W: null, ms: 0 });
      const pl = new Plate({ p: 3, ne: 4, edges: "SSSS" });
      const det = qoi === "deflection" ? pl.evaluate(pl.solve({ q: () => 1 }), 0.5, 0.5) : Math.sqrt(pl.modes(1).values[0]);
      const mean = qoi === "deflection" ? det * Math.exp(sigma ** 2) : det * Math.exp(-(sigma ** 2) / 8);
      expect(Math.abs(acc.q.mean - mean)).toBeLessThan(4 * acc.q.se);
    }
  });

  it("solves the same ω on a level and its parent, and records the midline", () => {
    const kl = karhunenLoeve("matern32", 0.3, 16);
    const s = new Sampler({ ...spec("deflection"), field: { ...field, kernel: "matern32", ell: 0.3, terms: 16 } }, kl);
    const r = s.sample(5, 8, true, true);
    expect(r.Qc).toBe(s.solve(5, 4).Q);
    expect(Math.abs(r.Q - r.Qc) / r.Q).toBeLessThan(1e-2);
    expect(r.w![(FIELD_POINTS - 1) / 2]).toBeCloseTo(r.Q, 14);
  });
});

describe("exact eigenvalues", () => {
  it("are offered exactly where a closed form exists", () => {
    expect(exactPlateEigenvalues("SSSS", 1, 1)).not.toBeNull();
    expect(exactPlateEigenvalues("CSCS", 1.5, 1)).not.toBeNull();
    expect(exactPlateEigenvalues("CCCC", 1, 1)).toBeNull();
  });
});
