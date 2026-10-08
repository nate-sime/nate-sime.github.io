// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The guided tours' tables. Every invariant `ui/tour.ts` relies on fails at
 * step time, part-way through a tour a reader clicks through once — so they
 * are checked here, without a browser, as the mantle app checks its own.
 */

import { describe, expect, it } from "vitest";
import { OPTIONS, RANGES } from "../src/ui/options";
import { defaultState, type State } from "../src/ui/state";
import { TOURS, TOUR_BASE, TOUR_NAMES, type TourStep } from "../src/ui/tours";
import { CONTROLS, visible } from "../src/ui/visibility";
import { PLATE_MAX_NE, admissibleAt, levelsWithin } from "../src/ui/views/structure";

const tours = Object.entries(TOURS) as [string, readonly TourStep[]][];
const allSteps = tours.flatMap(([tour, steps]) => steps.map((step, i) => ({ tour, steps, step, i, name: `${tour}[${i}] ${step.id}` })));

/** The state a step runs in: every patch up to and including it, over the defaults — what `enter` and `go` apply. */
function stateAt(steps: readonly TourStep[], i: number): State {
  const s = defaultState();
  for (let k = 0; k <= i; k++) Object.assign(s, steps[k].patch ?? {});
  return s;
}

describe("tour tables", () => {
  it("has the plan's four tours, each named by a button and none empty", () => {
    expect(TOUR_NAMES).toHaveLength(4);
    for (const [, steps] of tours) expect(steps.length).toBeGreaterThan(2);
  });

  it("gives every step a unique id", () => {
    const ids = allSteps.map((s) => s.step.id);
    expect(new Set(ids).size).toBe(ids.length);
  });

  it("starts every tour from the same base, so it lands the same whatever the reader had set", () => {
    for (const [tour, steps] of tours) {
      const first = steps[0].patch ?? {};
      for (const key of Object.keys(TOUR_BASE)) expect(first, `${tour}: ${key}`).toHaveProperty(key);
      expect(first.view, tour).toBeDefined();
    }
  });
});

describe("tour patches", () => {
  const defaults = defaultState();
  const options = OPTIONS as unknown as Record<string, Record<string, unknown>>;
  const ranges = RANGES as unknown as Record<string, { min: number; max: number; step: number }>;

  it.each(allSteps)("$name only sets fields State has, to values the pane can show", ({ step }) => {
    for (const [key, value] of Object.entries(step.patch ?? {})) {
      expect(defaults, key).toHaveProperty(key);
      expect(typeof value, key).toBe(typeof defaults[key as keyof State]);
      if (options[key]) expect(Object.values(options[key]), `${key} = ${String(value)}`).toContain(value);
      const r = ranges[key];
      if (r) {
        expect(value as number, key).toBeGreaterThanOrEqual(r.min);
        expect(value as number, key).toBeLessThanOrEqual(r.max);
      }
    }
  });
});

describe("tour steps in the state they leave", () => {
  it.each(allSteps)("$name lights a control that is on screen", ({ steps, step, i }) => {
    if (step.target === null) return;
    expect(CONTROLS).toContain(step.target);
    expect(visible(stateAt(steps, i), step.target), step.target).toBe(true);
  });

  it.each(allSteps)("$name runs on meshes the app will solve", ({ steps, i }) => {
    const st = stateAt(steps, i);
    const multilevel = st.view === "mlmc" || st.view === "live";
    if (st.view === "convergence" || multilevel)
      expect(admissibleAt(st, st.ne0)).toBeNull();
    if (st.view === "montecarlo") {
      const ne = st.ne0 * 2 ** st.mcLevel;
      expect(admissibleAt(st, st.mcLevel > 0 ? ne / 2 : ne)).toBeNull();
      if (st.structure === "plate") expect(ne).toBeLessThanOrEqual(PLATE_MAX_NE);
    }
    if (multilevel && st.structure === "plate")
      expect(levelsWithin(st.ne0, st.levels, PLATE_MAX_NE)).toBeGreaterThanOrEqual(3);
  });

  it.each(allSteps)("$name waits only for samples its run will make", ({ steps, step, i }) => {
    if (!step.dwell || !("samples" in step.dwell)) return;
    const st = stateAt(steps, i);
    expect(["montecarlo", "mlmc", "live"]).toContain(st.view);
    if (st.view === "montecarlo") expect(step.dwell.samples).toBeLessThanOrEqual(st.mcSamples);
    // A multilevel run stops when its sweep converges: the default beam sweep makes 248,480 samples in all.
    else expect(step.dwell.samples).toBeLessThanOrEqual(248_480);
  });
});
