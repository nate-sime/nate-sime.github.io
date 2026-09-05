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
import { DEFAULT_TOUR, TOUR_TARGETS, TOURS, type TourStep } from "../src/ui/tours";
import {
  BENCHMARKS, PARTICLES, QUICK_STARTS, RADIUS_INNER, VISCOSITY, defaultState,
  geometryFor, type State,
} from "../src/ui/presets";

const tours = Object.entries(TOURS) as [string, readonly TourStep[]][];
const allSteps = tours.flatMap(([tour, steps]) =>
  steps.map((step, i) => [`${tour}[${i}] ${step.id}`, step, steps] as const));

const presetTable: Record<string, Partial<State>> = { ...QUICK_STARTS, ...BENCHMARKS };

/**
 * Replay a tour up to and including step `i`, in the order `tour.ts`'s
 * `enter` applies things: preset first, then the patch over it. Several
 * checks below are about the state a step actually *runs* in, which is the
 * accumulation of every step before it, not the fields that step happens to
 * name itself.
 */
const stateAt = (steps: readonly TourStep[], i: number): State => {
  const s = defaultState();
  for (let k = 0; k <= i; k++) {
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
    for (const p of step.body) expect(p.trim().length).toBeGreaterThan(0);
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
  });

  // The slider's own bounds, `controls.ts`'s `vigour` binding: Tweakpane
  // clamps a value outside them on refresh, so a ramp aiming past the end
  // would animate to a number the pane then quietly changes underneath it.
  it.each(allSteps)("%s ramps within the vigour slider's range", (_name, step) => {
    if (!step.ramp) return;
    expect(step.ramp.to).toBeGreaterThanOrEqual(0);
    expect(step.ramp.to).toBeLessThanOrEqual(7);
    expect(step.ramp.ms).toBeGreaterThan(0);
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
    const g = geometryFor(stateAt(steps, steps.indexOf(step)));
    if (g.kind === "annulus") {
      const r = Math.hypot(x, y);
      expect(r).toBeGreaterThanOrEqual(g.lo);
      expect(r).toBeLessThanOrEqual(g.hi);
      expect(g.lo).toBeCloseTo(RADIUS_INNER, 6);
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
