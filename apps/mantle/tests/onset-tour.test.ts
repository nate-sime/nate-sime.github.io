// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The "Convection onset" tour's physics, checked against the settings each
 * step actually runs in.
 *
 * `tour.test.ts` checks that the tour is well-formed; this checks that what
 * it *says* is true. Every step of that tour makes a claim about the sign of
 * a growth rate — "here it fades", "above this it grows", "the same Ra that
 * grew now fades" — and each one rests on a threshold that a solver change,
 * a preset change or an edit to the step could move without anything else
 * noticing. So each claim is re-measured here from the step's own state,
 * replayed the way `tour.ts`'s `enter` accumulates it, with `geometryFor`
 * turning that state into the domain exactly as `main.ts` does.
 *
 * Only signs are asserted, and only where the measured rate is several units
 * from zero, on a grid far coarser than the app's: this is a check that the
 * prose points the right way, not a convergence study — `temperature.test.ts`
 * owns the analytic threshold itself.
 */

import { describe, it, expect } from "vitest";
import { Simulation } from "../src/solver/step";
import { TOURS, type TourStep } from "../src/ui/tours";
import { defaultState, geometryFor, type State } from "../src/ui/presets";

const steps: readonly TourStep[] = TOURS["Convection onset"];
const byId = (id: string): TourStep => {
  const s = steps.find((step) => step.id === id);
  if (!s) throw new Error(`no step "${id}" in the onset tour`);
  return s;
};

/** The state step `id` runs in: every patch up to and including its own. */
const stateAt = (id: string): State => {
  const s = defaultState();
  for (const step of steps) {
    if (step.patch) Object.assign(s, step.patch);
    if (step.id === id) return s;
  }
  throw new Error(`no step "${id}" in the onset tour`);
};

/**
 * Growth rate of the seeded pattern at `logRa` in `state`'s domain: the slope
 * of ln v_rms after the seed's other radial components have decayed. v_rms is
 * linear in the disturbance's amplitude while it is small, so this is σ
 * itself, not 2σ as a Nusselt-number measure would be.
 */
const sigma = (state: State, logRa: number): number => {
  const dt = state.dtMax;
  const sim = new Simulation({
    geom: geometryFor(state), nr: 16, na: 32, gnr: 33, gna: 64,
    Ra: 10 ** logRa, dtMax: dt, seed: { amp: 0.01, mode: state.wavenumber },
  });
  for (let n = 0; n < 50; n++) sim.step();
  const a = sim.temp.rmsVelocity(sim.velocity);
  for (let n = 0; n < 100; n++) sim.step();
  return Math.log(sim.temp.rmsVelocity(sim.velocity) / a) / (100 * dt);
};

describe("convection onset tour", () => {
  it("fades at the conduction step's Ra", () => {
    const s = stateAt("conduction");
    expect(sigma(s, s.logRa)).toBeLessThan(-2);
  });

  // The ramp is the step's whole point: it has to start on the stable side
  // and end on the unstable one, in the ring and pattern the tour seeds.
  it("ramps from below the ring's threshold to above it", () => {
    const s = stateAt("onset");
    const ramp = byId("onset").ramp!;
    expect(ramp.from).toBeDefined();
    expect(sigma(s, ramp.from!)).toBeLessThan(-2);
    expect(sigma(s, ramp.to)).toBeGreaterThan(2);
  });

  // "Each pattern has its own threshold … Mode 4 grows fastest … modes 2
  // and 6 fade", at the step's own Ra. Modes 3 and 5 are only claimed to
  // change slowly, which is not a sign to assert.
  it("grows mode 4 fastest at the pattern step's Ra, and fades modes 2 and 6", () => {
    const s = stateAt("patterns");
    expect(s.wavenumber).toBe(4);
    const four = sigma(s, s.logRa);
    expect(four).toBeGreaterThan(2);
    for (const m of [3, 5]) expect(sigma({ ...s, wavenumber: m }, s.logRa)).toBeLessThan(four);
    for (const m of [2, 6]) expect(sigma({ ...s, wavenumber: m }, s.logRa)).toBeLessThan(-1);
  });

  // "Bracket 657.5 with the card: at 600 they fade and at 720 they grow."
  it("starts the free-slip box below 27π⁴/4, and brackets it where the card says", () => {
    const s = stateAt("box-free-slip");
    expect(10 ** s.logRa).toBeLessThan((27 * Math.PI ** 4) / 4);
    expect(sigma(s, s.logRa)).toBeLessThan(-2);
    expect(sigma(s, Math.log10(600))).toBeLessThan(0);
    expect(sigma(s, Math.log10(720))).toBeGreaterThan(0);
  });

  // "Ra is about 800, where the rolls grew a moment ago" — and now fade.
  it("fades under no-slip at an Ra that grew under free slip", () => {
    const s = stateAt("box-no-slip");
    expect(s.radialWalls).toBe("no-slip");
    expect(sigma({ ...s, radialWalls: "free-slip" }, s.logRa)).toBeGreaterThan(2);
    expect(sigma(s, s.logRa)).toBeLessThan(-2);
  });

  // "Raise Ra past about 2,000 … and they grow again" — a quarter past
  // the quoted value, where the measured rate is clear of zero.
  it("grows under no-slip past the quoted threshold", () => {
    const s = stateAt("box-no-slip");
    expect(sigma(s, Math.log10(2000 * 1.25))).toBeGreaterThan(1);
  });
});
