// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { BEAM_POINTS, MADE_MODES, PLATE_POINTS, design, make } from "../src/manufacture";
import { REFERENCE } from "../src/plate/exact";

describe("the introduction's factory", () => {
  it("makes the design beam the uniform cantilever: ω₁ = 1.8751²", () => {
    const d = design("beam");
    expect(d.omega).toHaveLength(MADE_MODES);
    expect(d.omega[0]).toBeCloseTo(1.875104 ** 2, 4);
    expect(d.depth.every((h) => h === 1)).toBe(true);
    expect(d.shapes[0][BEAM_POINTS - 1]).toBeCloseTo(1, 12); // the tip, at its peak
  });

  it("makes the design plate Chladni's free square, its rigid modes skipped", () => {
    const d = design("plate"), ref = REFERENCE.FFFF!.omega.slice(3);
    for (let n = 0; n < MADE_MODES; n++) expect(Math.abs(d.omega[n] / ref[n] - 1)).toBeLessThan(5e-3);
    expect(d.shapes[0]).toHaveLength(PLATE_POINTS ** 2);
  });

  it("makes part i the same every time, and different parts differently", () => {
    for (const kind of ["beam", "plate"] as const) {
      const a = make(kind, 3), b = make(kind, 3), c = make(kind, 4);
      expect(a.omega).toEqual(b.omega);
      expect(a.depth).toEqual(b.depth);
      expect(a.omega[0]).not.toBe(c.omega[0]);
    }
  });

  it("signs each made mode to agree with the design's", () => {
    for (const kind of ["beam", "plate"] as const) {
      const d = design(kind), p = make(kind, 0);
      for (let n = 0; n < MADE_MODES; n++) {
        let dot = 0;
        for (let i = 0; i < d.shapes[n].length; i++) dot += d.shapes[n][i] * p.shapes[n][i];
        expect(dot).toBeGreaterThanOrEqual(0);
      }
    }
  });
});
