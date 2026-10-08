// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Which of the pane's controls each view shows — a pure function of `State`,
 * apart from the pane so it can be tested without a DOM. `controls.ts`
 * applies it after every change; `tests/tour.test.ts` uses it to check that
 * every control a tour step lights is actually on screen in the state that
 * step leaves the app in.
 *
 * A control is a folder or a blade inside one (`FOLDER_OF`); `hidden` is each
 * one's own flag, as Tweakpane holds it, and `visible` is what a reader sees:
 * the blade not hidden, and its folder not hidden either.
 */

import type { State } from "./state";

/** Every folder, and every blade that has a rule or that a tour may point at. */
export const CONTROLS = [
  // top level
  "tours", "view", "structure", "units",
  // folders
  "disc", "beam", "motion", "plate", "hierarchy", "field", "mc", "ml", "reference",
  // discretisation
  "p", "continuity", "ne", "derivative",
  // beam
  "supports", "section", "load", "mode", "polygon",
  // motion
  "motionShow", "motionForce", "motionZeta", "mesh",
  // plate
  "edges", "aspect", "nu", "plateLoad", "plateNe", "plateMode",
  // hierarchy
  "qoi", "ne0", "levels", "force", "zeta",
  // random field
  "kernel", "ell", "sigma", "terms", "massFollows", "loadSigma", "seed",
  // Monte Carlo, multilevel Monte Carlo
  "mcLevel", "mcSamples", "mlSurvey", "mlEps", "liveLevel",
] as const;

export type ControlName = (typeof CONTROLS)[number];

/** The folder each blade lives in; absent for the top level and the folders themselves. */
export const FOLDER_OF: Partial<Record<ControlName, ControlName>> = {
  p: "disc", continuity: "disc", ne: "disc", derivative: "disc",
  supports: "beam", section: "beam", load: "beam", mode: "beam", polygon: "beam",
  motionShow: "motion", motionForce: "motion", motionZeta: "motion", mesh: "motion",
  edges: "plate", aspect: "plate", nu: "plate", plateLoad: "plate", plateNe: "plate", plateMode: "plate",
  qoi: "hierarchy", ne0: "hierarchy", levels: "hierarchy", force: "hierarchy", zeta: "hierarchy",
  kernel: "field", ell: "field", sigma: "field", terms: "field", massFollows: "field", loadSigma: "field", seed: "field",
  mcLevel: "mc", mcSamples: "mc", mlSurvey: "ml", mlEps: "ml", liveLevel: "ml",
};

/** Each control's own hidden flag for this state. */
export function hidden(st: State): Record<ControlName, boolean> {
  const v = st.view;
  const sampling = v === "montecarlo" || v === "mlmc" || v === "live";
  // The views that take either structure, and which one they show.
  const either = v === "convergence" || v === "field" || sampling;
  const onPlate = v === "plate" || (either && st.structure === "plate");
  const beamish = v === "beam" || v === "spectrum" || (either && !onPlate);
  const hier = v !== "convergence" && !sampling;
  const h: Record<ControlName, boolean> = Object.fromEntries(CONTROLS.map((c) => [c, false])) as Record<ControlName, boolean>;
  Object.assign(h, {
    structure: !either,
    units: v === "basis",
    derivative: v !== "basis",
    ne: v === "convergence" || v === "plate" || sampling || (v === "field" && onPlate),
    beam: !beamish,
    supports: v === "field",
    mode: v !== "beam" || st.motion !== "mode",
    polygon: v !== "beam",
    load: v === "spectrum",
    motion: v !== "beam" && v !== "plate",
    motionForce: st.motion !== "response",
    motionZeta: st.motion !== "response",
    plate: !onPlate,
    plateNe: !(v === "plate" || v === "field"),
    plateMode: v !== "plate" || st.motion !== "mode",
    plateLoad: v === "field",
    nu: v === "field",
    hierarchy: hier,
    levels: v === "montecarlo",
    force: hier || st.qoi !== "response",
    zeta: hier || st.qoi !== "response",
    field: v !== "field" && !sampling,
    mc: v !== "montecarlo",
    ml: v !== "mlmc" && v !== "live",
    liveLevel: v !== "live",
    reference: !st.dimensional,
  } satisfies Partial<Record<ControlName, boolean>>);
  return h;
}

/** Whether a reader can see the control: neither it nor its folder hidden. */
export function visible(st: State, name: ControlName): boolean {
  const h = hidden(st), folder = FOLDER_OF[name];
  return !h[name] && !(folder && h[folder]);
}
