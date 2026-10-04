// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The guided tour's step tables. No GPU and no DOM — `vitest.config.ts` sets
 * no environment, and none of this needs one: these are the invariants
 * `ui/tour.ts` assumes and the pane silently supplies, and every one of them
 * fails part-way through a tour rather than at start-up. A bad step could
 * ship behind nine good ones nobody clicked past.
 *
 * The same trade `presets.test.ts` makes for `QUICK_STARTS`, for the same
 * reason.
 */

import { describe, it, expect } from "vitest";
import {
  DEFAULT_TOUR, SECTION_HELP, TOUR_TARGETS, TOURS, type TourStep,
} from "../src/ui/tours";
import {
  BENCHMARKS, BOX_LENGTH, GEOMETRY, LOG_RA, PARTICLES, PRESETS, QUICK_STARTS, RADIAL_WALLS,
  DEFAULT_PRESET, RADIUS_INNER, SIGMA_Y, SPEEDS, VISCOSITY, defaultState, geometryFor,
  type State,
} from "../src/ui/presets";
import { planetFor } from "../src/planets";

const tours = Object.entries(TOURS) as [string, readonly TourStep[]][];
const allSteps = tours.flatMap(([tour, steps]) =>
  steps.map((step, i) => [`${tour}[${i}] ${step.id}`, step, steps] as const));

const presetTable: Record<string, Partial<State>> = { ...QUICK_STARTS, ...BENCHMARKS };

/**
 * Replay a tour up to and including step `i`, in the order `tour.ts`'s
 * `stage` applies things: a planet not already live (its profile, then the
 * step's patch over it, as `selectPlanet` writes them), then the preset, then
 * the patch. Several checks below are about the state a step actually *runs*
 * in, which is the accumulation of every step before it, not the fields that
 * step happens to name itself.
 */
const stateAt = (steps: readonly TourStep[], i: number): State => {
  const s = defaultState();
  for (let k = 0; k <= i; k++) {
    const planet = steps[k].planet;
    if (planet && s.activePlanet !== planet) {
      const profile = planetFor(planet);
      Object.assign(s, profile.solver.state, {
        isothermal: false, wavenumber: profile.solver.initialWavenumber, activePlanet: planet,
      });
    }
    if (steps[k].preset) Object.assign(s, presetTable[steps[k].preset!]);
    if (steps[k].patch) Object.assign(s, steps[k].patch);
  }
  return s;
};

/** The view a step runs in, carried forward the same way `enter` reconciles it. */
const viewAt = (steps: readonly TourStep[], i: number): "2d" | "3d" => {
  // The app opens on the annulus in 3D — see `build` in main.ts.
  let v: "2d" | "3d" = "3d";
  for (let k = 0; k <= i; k++) if (steps[k].view) v = steps[k].view!;
  return v;
};

describe("tour tables", () => {
  it("opens on a tour that exists", () => {
    expect(TOURS[DEFAULT_TOUR]).toBeDefined();
    expect(TOURS[DEFAULT_TOUR].length).toBeGreaterThan(0);
  });

  it.each(tours)("%s has unique step ids", (_tour, steps) => {
    const ids = steps.map((s) => s.id);
    expect(new Set(ids).size).toBe(ids.length);
  });

  // A step with nothing to say is a step that dims the screen for no reason.
  it.each(allSteps)("%s says something", (_name, step) => {
    expect(step.title.length).toBeGreaterThan(0);
    expect(step.body.length).toBeGreaterThan(0);
    for (const p of step.body) {
      if (typeof p === "string") { expect(p.trim().length).toBeGreaterThan(0); continue; }
      if ("items" in p) {
        expect(p.items.length).toBeGreaterThan(0);
        for (const item of p.items) expect(item.trim().length).toBeGreaterThan(0);
        continue;
      }
      // A table: a header, rows, and no ragged or empty cells.
      expect(p.head.length).toBeGreaterThan(0);
      expect(p.rows.length).toBeGreaterThan(0);
      for (const row of p.rows) expect(row.length).toBe(p.head.length);
      for (const cell of [...p.head, ...p.rows.flat()]) expect(cell.trim().length).toBeGreaterThan(0);
    }
    if (step.watch !== undefined) expect(step.watch.trim().length).toBeGreaterThan(0);
  });
});

describe("tour targets", () => {
  // `main.ts` types its resolver as `Record<TourTargetName, HTMLElement>`, so
  // the compiler already insists every *listed* name is wired to something.
  // This is the other direction: that a step only ever names one on the list.
  it.each(allSteps)("%s points at a target the app exposes", (_name, step) => {
    if (step.target === null) return;
    expect(TOUR_TARGETS).toContain(step.target);
    if (step.highlight !== undefined) expect(TOUR_TARGETS).toContain(step.highlight);
    if (step.companion !== undefined) expect(TOUR_TARGETS).toContain(step.companion);
  });

  // Not a style rule: `tour.ts` lights exactly one element per step, and
  // `holeFor` falls back to dimming the whole screen when it cannot resolve
  // one. A step that meant to point somewhere and silently didn't looks
  // identical to an opening card.
  it("points somewhere on every step but the framing ones", () => {
    const steps = TOURS[DEFAULT_TOUR];
    const blind = steps.filter((s) => s.target === null);
    expect(blind.length).toBeLessThanOrEqual(1);
  });
});

describe("tour actions", () => {
  // The "?" buttons live on the advanced folders, which the simple view
  // hides outright — a step pointing at one in the simple view lights
  // nothing, and falls back to an unanchored card.
  it.each(allSteps)("%s shows the advanced view to point at a section's ?", (_name, step) => {
    if (step.target !== "domainHelp" && step.highlight !== "domainHelp") return;
    expect(step.advanced).toBe(true);
  });

  // The tour's physics was measured on the default grid at Courant 2 (see
  // `STRONG_LID` in tours.ts); the lids' regimes depend on both.
  it("runs every three-planet step on the grid its lids were measured on", () => {
    const steps: readonly TourStep[] = TOURS["Three planet tour"];
    for (let i = 0; i < steps.length; i++) {
      const s = stateAt(steps, i);
      expect(s.resolution).toBe(DEFAULT_PRESET);
      expect(s.courant).toBe(2);
    }
  });

  // `stage` skips a planet that is already live, so stating it is free, and
  // without it "back" from the next planet lands a card on the wrong one.
  it("names the planet on every three-planet step that changes the model", () => {
    const steps: readonly TourStep[] = TOURS["Three planet tour"];
    for (const step of steps)
      if (step.patch) expect(step.planet).toBeDefined();
  });

  // The Venus → Earth comparison is the tour's claim: same rock, one change.
  // Ra moves too, but only by the g·d³ the card's table accounts for.
  it("changes only the yield stress and Ra between Venus and Earth", () => {
    const steps: readonly TourStep[] = TOURS["Three planet tour"];
    const venus = stateAt(steps, steps.findIndex((s) => s.id === "venus"));
    const earth = stateAt(steps, steps.findIndex((s) => s.id === "earth"));
    const differ = (Object.keys(venus) as (keyof State)[]).filter((k) => venus[k] !== earth[k]);
    expect(differ.sort()).toEqual(["activePlanet", "logRa", "sigmaY"]);
    expect(earth.sigmaY).toBeLessThan(venus.sigmaY);
    expect(Math.abs(earth.logRa - venus.logRa)).toBeLessThan(0.05);
  });

  // "finish" puts these back from the snapshot taken at open; anywhere but
  // the last card there is no "finish" to press.
  it.each(allSteps)("%s restores on finish only from its last step", (_name, step, steps) => {
    if (!step.restoreOnFinish) return;
    expect(steps.indexOf(step)).toBe(steps.length - 1);
  });

  it.each(allSteps)("%s names a preset that exists", (_name, step) => {
    if (!step.preset) return;
    expect(presetTable[step.preset]).toBeDefined();
  });

  // `applyPatch` (controls.ts) writes these onto the live `State` and then
  // dispatches by key, so a value outside its own table reaches the solver as
  // an `undefined` lookup rather than being rejected at the pane.
  it.each(allSteps)("%s patches only legal values", (_name, step) => {
    const patch = step.patch;
    if (!patch) return;
    if (patch.viscosity !== undefined) expect(VISCOSITY[patch.viscosity]).toBeDefined();
    if (patch.particles !== undefined) expect(PARTICLES[patch.particles]).toBeDefined();
    if (patch.contours !== undefined) expect(patch.contours).toBeGreaterThanOrEqual(0);
    if (patch.geometry !== undefined) expect(GEOMETRY[patch.geometry]).toBeDefined();
    if (patch.radialWalls !== undefined) expect(RADIAL_WALLS[patch.radialWalls]).toBeDefined();
    if (patch.resolution !== undefined) expect(PRESETS[patch.resolution]).toBeDefined();
    if (patch.speed !== undefined) expect(Object.values(SPEEDS)).toContain(patch.speed);
    if (patch.boxLength !== undefined) {
      expect(patch.boxLength).toBeGreaterThanOrEqual(BOX_LENGTH.min);
      expect(patch.boxLength).toBeLessThanOrEqual(BOX_LENGTH.max);
    }
    // Tweakpane clamps to the binding's bounds on refresh.
    if (patch.sigmaY !== undefined) {
      expect(patch.sigmaY).toBeGreaterThanOrEqual(SIGMA_Y.min);
      expect(patch.sigmaY).toBeLessThanOrEqual(SIGMA_Y.max);
    }
  });

  // The slider's own bounds, `controls.ts`'s `vigour` binding: Tweakpane
  // clamps a value outside them on refresh, so a ramp aiming past the end
  // would animate to a number the pane then quietly changes underneath it.
  it.each(allSteps)("%s ramps within the vigour slider's range", (_name, step) => {
    if (!step.ramp) return;
    for (const v of [step.ramp.to, step.ramp.from ?? step.ramp.to]) {
      expect(v).toBeGreaterThanOrEqual(0);
      expect(v).toBeLessThanOrEqual(7);
    }
    expect(step.ramp.ms).toBeGreaterThan(0);
  });

  // The card's own Ra slider (`raControl`) is a window onto the pane's: it has
  // to sit inside the pane slider's range, and it has to reach the values its
  // step sets, or it would open resting against one end.
  it.each(allSteps)("%s gives its Ra control a range covering the step", (_name, step, steps) => {
    if (!step.raControl) return;
    const { min, max } = step.raControl;
    expect(min).toBeLessThan(max);
    expect(min).toBeGreaterThanOrEqual(LOG_RA.min);
    expect(max).toBeLessThanOrEqual(LOG_RA.max);
    const set = [stateAt(steps, steps.indexOf(step)).logRa, step.ramp?.from, step.ramp?.to];
    for (const v of set) {
      if (v === undefined) continue;
      expect(v).toBeGreaterThanOrEqual(min);
      expect(v).toBeLessThanOrEqual(max);
    }
  });

  // `isothermal` forces Ra to 0 regardless of the slider (see that flag's own
  // header in presets.ts) and hides the vigour blade outright — a ramp under
  // it would drag a control that is neither visible nor being solved with.
  it.each(allSteps)("%s does not ramp under the isothermal override", (_name, step, steps) => {
    if (!step.ramp) return;
    const i = steps.indexOf(step);
    expect(stateAt(steps, i).isothermal).toBe(false);
  });
});

describe("tour camera", () => {
  // `focus` drives the flat view's camera (`animateViewTo` in main.ts). The
  // globe has its own, which that function does not touch, so a focus step
  // still in 3D would narrate a zoom that never happened.
  it.each(allSteps)("%s only focuses in the flat view", (_name, step, steps) => {
    if (!step.focus) return;
    expect(viewAt(steps, steps.indexOf(step))).toBe("2d");
  });

  // A point outside the domain frames background. Checked against the
  // geometry the step is actually running in rather than a literal, so
  // changing `RADIUS_INNER` or moving a step onto a box fails here.
  it.each(allSteps)("%s focuses on a point inside the domain", (_name, step, steps) => {
    if (!step.focus || step.focus === "reset") return;
    const { zoom, x, y } = step.focus;
    expect(zoom).toBeGreaterThanOrEqual(1);
    expect(zoom).toBeLessThanOrEqual(40);     // ZOOM_MIN/ZOOM_MAX, main.ts
    const s = stateAt(steps, steps.indexOf(step));
    const g = geometryFor(s);
    if (g.kind === "annulus") {
      const r = Math.hypot(x, y);
      expect(r).toBeGreaterThanOrEqual(g.lo);
      expect(r).toBeLessThanOrEqual(g.hi);
      if ((s.activePlanet ?? "earth") === "earth") expect(g.lo).toBeCloseTo(RADIUS_INNER, 6);
    } else {
      expect(Math.abs(x)).toBeLessThanOrEqual(g.width);
      expect(y).toBeGreaterThanOrEqual(g.lo);
      expect(y).toBeLessThanOrEqual(g.hi);
    }
  });
});

describe("tour dwells", () => {
  it.each(allSteps)("%s dwells for a positive amount", (_name, step) => {
    if (!step.dwell) return;
    const d = step.dwell;
    expect("steps" in d ? d.steps : d.ms).toBeGreaterThan(0);
  });

  // A dwell counted in solver steps never fills while the solver is paused.
  // It doesn't strand a reader — `tour.ts` gates nothing on it — but the bar
  // would sit at zero under a line promising an effect that cannot arrive.
  it.each(allSteps)("%s is running if it counts solver steps", (_name, step, steps) => {
    if (!step.dwell || !("steps" in step.dwell)) return;
    const s = stateAt(steps, steps.indexOf(step));
    expect(s.paused).toBe(false);
    expect(s.speed).toBeGreaterThan(0);
  });
});

describe("section help", () => {
  const sections = Object.entries(SECTION_HELP) as [string, readonly TourStep[]][];
  const helpSteps = sections.flatMap(([section, steps]) =>
    steps.map((step, i) => [`${section}[${i}] ${step.id}`, step] as const));

  // `tour.ts` looks both tables up by name in one merged record, so a section
  // sharing a tour's name would silently replace one with the other.
  it("shares no name with a guided tour", () => {
    for (const name of Object.keys(SECTION_HELP)) expect(Object.keys(TOURS)).not.toContain(name);
  });

  it.each(sections)("%s has steps with unique ids", (_section, steps) => {
    expect(steps.length).toBeGreaterThan(0);
    const ids = steps.map((s) => s.id);
    expect(new Set(ids).size).toBe(ids.length);
  });

  it.each(helpSteps)("%s says something about a control the app exposes", (_name, step) => {
    expect(step.title.trim().length).toBeGreaterThan(0);
    expect(step.body.length).toBeGreaterThan(0);
    for (const p of step.body) if (typeof p === "string") expect(p.trim().length).toBeGreaterThan(0);
    if (step.watch !== undefined) expect(step.watch.trim().length).toBeGreaterThan(0);
    expect(step.target).not.toBeNull();
    expect(TOUR_TARGETS).toContain(step.target);
    if (step.highlight !== undefined) expect(TOUR_TARGETS).toContain(step.highlight);
  });

  // `SectionHelpStep` already leaves these fields out of the type; this is
  // the same promise checked on the data, so a cast cannot quietly let a
  // help card start changing the run it was opened over.
  it.each(helpSteps)("%s only explains, never changes the run", (_name, step) => {
    for (const key of [
      "patch", "preset", "planet", "ramp", "reseed", "focus", "view",
      "surfaceGuide", "dwell", "advanced", "replay", "restoreOnFinish",
    ] as const) {
      expect(step[key]).toBeUndefined();
    }
  });
});
