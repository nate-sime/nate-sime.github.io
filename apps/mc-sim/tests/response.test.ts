// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { complexSolve } from "../src/band";
import { Beam, SUPPORTS } from "../src/beam/beam";
import { pinnedResponse } from "../src/beam/exact";
import { denseGeneralizedEig } from "../src/eig";
import { evaluateQoI, harmonicOf, runHierarchy, type BeamCase } from "../src/hierarchy";

const pinned = (load: BeamCase["load"], ratio = 0.8, zeta = 0.02): BeamCase =>
  ({ supports: "pinned–pinned", load, section: "uniform", forcing: { ratio, zeta } });

describe("the forced, damped response", () => {
  it("solves the complex symmetric system exactly as modal superposition does", () => {
    // Rayleigh damping is diagonal in the undamped modes, so on the discrete
    // problem u = Σ φ_n (φ_nᵀF) / (λ_n − Ω² + iΩ(a + bλ_n)) is an identity.
    for (const supports of ["cantilever", "clamped–pinned"] as const) {
      const beam = new Beam({ p: 3, ne: 12, supports: SUPPORTS[supports], stiffness: (x) => 1 + 0.5 * x, mass: (x) => 2 - x });
      const h = harmonicOf({ supports, load: "uniform", section: "uniform", forcing: { ratio: 1.3, zeta: 0.05 } });
      const F = beam.load({ q: (x) => Math.cos(3 * x) }).subarray(beam.lo, beam.hi);
      const W = h.Omega;
      const u = complexSolve(beam.K, beam.M, [1, W * h.b], [-W * W, W * h.a], F);
      const { values, vectors } = denseGeneralizedEig(beam.K.toDense(), beam.M.toDense());
      const re = new Float64Array(beam.dofs), im = new Float64Array(beam.dofs);
      values.forEach((lam, n) => {
        const f = vectors[n].reduce((s, v, i) => s + v * F[i], 0);
        const dr = lam - W * W, di = W * (h.a + h.b * lam), d2 = dr * dr + di * di;
        vectors[n].forEach((v, i) => { re[i] += (v * f * dr) / d2; im[i] -= (v * f * di) / d2; });
      });
      // Both carry round-off of order ε·cond(K): measured ~10⁻¹¹ at 12 cubic elements.
      const scale = Math.max(...Array.from(re, Math.abs), ...Array.from(im, Math.abs));
      for (let i = 0; i < beam.dofs; i++) {
        expect(Math.abs(u.re[i] - re[i]) / scale).toBeLessThan(1e-10);
        expect(Math.abs(u.im[i] - im[i]) / scale).toBeLessThan(1e-10);
      }
    }
  });

  it("fits Rayleigh damping to ζ on the first two modes", () => {
    for (const supports of Object.keys(SUPPORTS) as (keyof typeof SUPPORTS)[]) {
      const bc: BeamCase = { supports, load: "uniform", section: "uniform", forcing: { ratio: 0.5, zeta: 0.03 } };
      const h = harmonicOf(bc), beam = new Beam({ p: 4, ne: 32, supports: SUPPORTS[supports] });
      const w = beam.modes(2).values.map(Math.sqrt);
      for (const wn of w) expect(h.a / (2 * wn) + (h.b * wn) / 2).toBeCloseTo(0.03, 6);
      expect(h.Omega / w[0]).toBeCloseTo(0.5, 6);
    }
  });

  it("reduces to the static deflection as Ω → 0, and peaks at resonance by 1/2ζ", () => {
    expect(pinnedResponse("uniform", { Omega: 0, a: 0, b: 0 })).toBeCloseTo(5 / 384, 14);
    expect(pinnedResponse("point", { Omega: 0, a: 0, b: 0 })).toBeCloseTo(1 / 48, 9);
    const beam = new Beam({ p: 3, ne: 16, supports: SUPPORTS["pinned–pinned"] });
    const slow = evaluateQoI(beam, "response", { q: () => 1 }, 0.5, harmonicOf(pinned("uniform", 1e-6)));
    expect(slow.Q / evaluateQoI(beam, "deflection", { q: () => 1 }, 0.5).Q).toBeCloseTo(1, 10);
    // At Ω = ω₁ the first mode alone amplifies its static share by 1/(2ζ).
    const zeta = 0.01, at = pinnedResponse("uniform", harmonicOf(pinned("uniform", 1, zeta)));
    const mode1 = 4 / Math.PI / Math.PI ** 4; // φ₁(½) f₁ / λ₁
    expect(at / (mode1 / (2 * zeta))).toBeCloseTo(1, 2);
  });

  it("converges to the modal series at h^(2(p−1)) for a distributed load", () => {
    for (const p of [2, 3]) {
      const H = runHierarchy({ p, k: p - 1, beam: pinned("uniform"), qoi: "response", ne0: 4, levels: 5 });
      expect(H.exact).toBeCloseTo(pinnedResponse("uniform", harmonicOf(pinned("uniform"))), 15);
      expect(H.alphaExact!).toBeCloseTo(2 * (p - 1), 1);
      expect(H.alpha!).toBeCloseTo(2 * (p - 1), 1);
    }
  });
});
