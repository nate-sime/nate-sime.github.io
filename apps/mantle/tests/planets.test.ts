/** The profile registry's Earth extraction: data must reproduce the old app. */

import { describe, expect, it } from "vitest";
import { EARTH, isPlanetProfileModified, PLANETS, planetFor, radiiFor } from "../src/planets";
import { SURFACE_MATERIALS } from "../src/gpu/surfaceAssets";
import { createDimensionalScale, timeUnitFor, velocityUnitFor } from "../src/ui/dimensional";
import { defaultState, geometryFor } from "../src/ui/presets";

describe("planet profiles", () => {
  it("ships one uniquely named, documented Earth profile", () => {
    expect(Object.keys(PLANETS)).toEqual(["earth"]);
    expect(planetFor("earth")).toBe(EARTH);
    expect(EARTH.model.sources.length).toBeGreaterThan(0);
    expect(EARTH.model.caveats.length).toBeGreaterThan(0);
    expect(SURFACE_MATERIALS[EARTH.visual.surface]).toBeDefined();
  });

  it("derives the existing Earth annulus from physical radii", () => {
    const r = radiiFor(EARTH);
    expect(r.depthKm).toBe(2885);
    expect(r.ri).toBeCloseTo(1.208318891, 9);
    expect(r.ro - r.ri).toBeCloseTo(1, 12);
    const g = geometryFor(defaultState());
    expect(g.lo).toBeCloseTo(r.ri, 12);
    expect(g.hi).toBeCloseTo(r.ro, 12);
  });

  it("derives profile modification from planet-owned controls", () => {
    const state = defaultState();
    expect(isPlanetProfileModified(state)).toBe(false);
    state.logRa += 1;
    expect(isPlanetProfileModified(state)).toBe(true);
  });

  it("derives time and velocity display units from the same physical scale", () => {
    const r = radiiFor(EARTH);
    const scale = createDimensionalScale(r.depthKm * 1e3, EARTH.physical.thermalDiffusivityM2s);
    expect(timeUnitFor(scale)).toBeCloseTo(scale.depth ** 2 / scale.kappa, 6);
    expect(velocityUnitFor(scale)).toBeCloseTo(scale.kappa / scale.depth, 20);
  });
});
