// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { fitRate, runHierarchy, theoryRate, type BeamCase } from "../src/hierarchy";
import {
  STEEL_BAR, deflectionScale, flexuralRigidity, fmt, frequencyScale, massPerLength, si,
} from "../src/ui/units";

const uniform = (supports: BeamCase["supports"], load: BeamCase["load"] = "uniform"): BeamCase =>
  ({ supports, load, section: "uniform" });

describe("discretisation hierarchy", () => {
  it("fits a known slope and ignores points at the round-off floor", () => {
    const h = [1 / 4, 1 / 8, 1 / 16, 1 / 32, 1 / 64], none = h.map(() => 0);
    expect(fitRate(h, h.map((x) => 3 * x ** 4), none)).toBeCloseTo(4, 10);
    // The last two sit on their floor: the fit uses the three before them.
    const floor = [1e-15, 1e-15, 1e-15, 1e-15, 1e-15];
    expect(fitRate(h, [3e-3, 3e-3 / 16, 3e-3 / 256, 1e-15, 2e-15], floor)).toBeCloseTo(4, 10);
    expect(fitRate(h, h.map(() => 1e-16), floor)).toBeNull();
  });

  it("measures a round-off floor that rises as the mesh refines, and fits only above it", () => {
    // The cantilever's λ₁ is small beside its largest discrete eigenvalue, so
    // f64 round-off reaches ~10⁻¹⁰ by a few dozen elements at p = 4.
    const H = runHierarchy({ p: 4, k: 3, beam: uniform("cantilever"), qoi: "omega1", ne0: 4, levels: 6 });
    const noise = H.levels.map((l) => l.noise / H.exact!);
    expect(noise[5]).toBeGreaterThan(10 * noise[0]);
    expect(H.levels[5].err / H.exact!).toBeGreaterThan(1e-12); // at the floor, not at h⁶ below it
    expect(H.alphaExact!).toBeCloseTo(6, 0);
  });

  it("measures α = 2(p − 1) for ω₁, in both the error and the successive differences", () => {
    for (const p of [2, 3, 4])
      for (const k of [1, p - 1]) {
        const H = runHierarchy({ p, k, beam: uniform("clamped–clamped"), qoi: "omega1", ne0: 4, levels: 4 });
        const want = theoryRate("omega1", p, uniform("clamped–clamped"))!;
        expect(H.alphaExact!).toBeCloseTo(want, 0);
        expect(H.alpha!).toBeCloseTo(want, 0);
        // The difference between levels bounds the error of the coarser one, to
        // within the factor 1 − 2^−α the geometric tail predicts.
        const [, l1] = H.levels;
        expect(l1.dQ / H.levels[0].err).toBeCloseTo(1 - 2 ** -want, 1);
        // dofs = m·ne + c: γ → 1, from above while the offset c still shows.
        expect(H.gamma!).toBeGreaterThan(0.95);
        expect(H.gamma!).toBeLessThan(1.25);
      }
  });

  it("measures the field (RMS deflection) rate on a tapered beam with no closed form", () => {
    const tapered: BeamCase = { supports: "cantilever", load: "uniform", section: "tapered" };
    for (const p of [2, 3]) {
      const H = runHierarchy({ p, k: p - 1, beam: tapered, qoi: "field", ne0: 4, levels: 5 });
      expect(H.exact).toBeNull();
      expect(H.levels.every((l) => Number.isNaN(l.err))).toBe(true);
      expect(H.alpha!).toBeCloseTo(theoryRate("field", p, tapered)!, 0);
    }
  });

  it("finds nothing left to converge when the solution is in the space", () => {
    // Uniform load, quartic solution, quartic splines: every level is exact.
    const H = runHierarchy({ p: 4, k: 3, beam: uniform("pinned–pinned"), qoi: "compliance", ne0: 2, levels: 4 });
    for (const l of H.levels) expect(l.err / H.exact!).toBeLessThan(1e-12);
    expect(H.alphaExact).toBeNull();
  });

  it("gives the textbook values on the finest level", () => {
    const H = runHierarchy({ p: 3, k: 2, beam: uniform("cantilever", "point"), qoi: "deflection", ne0: 4, levels: 2 });
    expect(H.exact).toBeCloseTo(1 / 3, 14);
    expect(H.levels[1].Q).toBeCloseTo(1 / 3, 12);
  });
});

describe("dimensional display", () => {
  it("puts the steel bar's first frequency where a hand calculation does", () => {
    const r = STEEL_BAR, EI = flexuralRigidity(r), rA = massPerLength(r);
    expect(EI).toBeCloseTo(2800, 6);
    expect(rA).toBeCloseTo(3.14, 10);
    // Cantilever: f₁ = (1.8751²/2π) √(EI / ρA L⁴).
    const f1 = (1.87510407 ** 2 / (2 * Math.PI)) * Math.sqrt(EI / (rA * r.L ** 4));
    expect((1.87510407 ** 2 * frequencyScale(r)) / (2 * Math.PI)).toBeCloseTo(f1, 10);
    expect(fmt.frequency({ dimensional: true, ref: r }, 1.87510407 ** 2)).toBe(si(f1, "Hz"));
    // Tip deflection under its own reference load: qL⁴/8EI.
    expect(deflectionScale(r, "uniform") / 8).toBeCloseTo((r.q0 * r.L ** 4) / (8 * EI), 14);
  });

  it("formats with engineering prefixes, and leaves nondimensional numbers bare", () => {
    expect(si(0.004464, "m")).toBe("4.464 mm");
    expect(si(16.73, "Hz")).toBe("16.73 Hz");
    expect(si(2800, "N·m²", 3)).toBe("2.80 kN·m²");
    expect(fmt.length({ dimensional: false, ref: STEEL_BAR }, 0.5)).toBe("0.5000");
  });
});
