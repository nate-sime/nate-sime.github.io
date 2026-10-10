// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Entry point: one figure, one readout, one pane, and a view per stage of
 * MC_PLAN.md that has landed —
 *
 *   1  the spline basis, and its tensor product on the plate's elements
 *                                 (`ui/views/basis.ts`)
 *   2  the beam and the Kirchhoff plate, static and modal, side by side
 *                                 (`ui/views/structures.ts`, drawing `ui/views/beam.ts`
 *                                 and `ui/views/plate.ts`; the plate's solver in `plate/`)
 *   3  the discretisation hierarchy (`ui/views/convergence.ts`)
 *   4  the random input           (`ui/views/field.ts`)
 *   5  plain Monte Carlo          (`ui/views/montecarlo.ts`, sampling in `mc/`)
 *   6  multilevel Monte Carlo     (`ui/views/mlmc.ts`, the algorithm in `mc/mlmc.ts`)
 *   7  MLMC, live                 (`ui/views/live.ts`): a run's coupled samples,
 *                                 solved and moving, beside its statistics as they
 *                                 accumulate; view 2 animates the forced response
 *
 * The hierarchy, field and sampling views take the plate too, when the pane's
 * structure is the plate.
 *
 * A view draws every graph it has at once, each in a panel of the figure
 * (`ui/figure.ts`). Everything is computed on the CPU in f64 and drawn into 2-D
 * canvases; Monte Carlo samples run in a pool of Web Workers, paused while
 * their view is not on screen.
 * The pane's first folder opens four guided tours (`ui/tour.ts`, `ui/tours.ts`)
 * that walk through the stages' lessons, setting the pane as they go.
 * Frames are requested only while something moves (a vibrating mode, a run);
 * otherwise a redraw follows a control change, a resize, or the pointer.
 * `U` toggles dimensional units, as the pane's button does.
 */

import { buildPane } from "./ui/controls";
import { buildTour } from "./ui/tour";
import type { TourName } from "./ui/tours";
import { Figure } from "./ui/figure";
import { defaultState, type View } from "./ui/state";
import { renderBasis } from "./ui/views/basis";
import { renderConvergence } from "./ui/views/convergence";
import { renderField } from "./ui/views/field";
import { renderMlmc } from "./ui/views/mlmc";
import { renderMonteCarlo } from "./ui/views/montecarlo";
import { renderLive } from "./ui/views/live";
import { renderStructures } from "./ui/views/structures";
import type { ViewResult } from "./ui/views/view";
import { workers } from "./ui/views/workers";

const el = (id: string) => document.getElementById(id)!;

el("caption").textContent = `MC-sim ${__APP_VERSION__} · Monte Carlo for vibrating structures`;

const state = defaultState();
const figure = new Figure(el("figure"));
const readout = el("readout");

const RENDER: Record<View, (t: number) => ViewResult> = {
  basis: () => renderBasis(figure, state),
  structures: (t) => renderStructures(figure, state, t),
  convergence: () => renderConvergence(figure, state),
  field: () => renderField(figure, state),
  montecarlo: () => renderMonteCarlo(figure, state),
  mlmc: () => renderMlmc(figure, state),
  live: (t) => renderLive(figure, state, t),
};

/**
 * The readout keeps the tallest height it has had in this view. A rerun cuts
 * its text to one "waiting" line and then grows it back, and the figure takes
 * whatever height the readout leaves: letting it shrink would resize, and
 * blank, every panel twice. A new view, or a resized window, starts it afresh.
 */
let readoutMin = 0;
const resetReadout = () => { readoutMin = 0; readout.style.minHeight = ""; };
function holdReadout(): void {
  const h = readout.offsetHeight;
  if (h > readoutMin) { readoutMin = h; readout.style.minHeight = `${h}px`; }
}

let frame = 0;
let shown: View | null = null;
function render(t = performance.now()): void {
  frame = 0;
  // Leaving a sampling view idles the workers; the next view to sample opens its own run.
  if (state.view !== shown) {
    workers.idle();
    shown = state.view;
    resetReadout();
  }
  let r: ViewResult;
  try {
    r = RENDER[state.view](t);
  } catch (e) {
    r = { readout: `error: ${(e as Error).message}`, animate: false };
  }
  if (readout.textContent !== r.readout) readout.textContent = r.readout;
  holdReadout();
  if (r.animate) frame = requestAnimationFrame(render);
}

const request = () => { if (!frame) frame = requestAnimationFrame(render); };
workers.connect(request);
// The tour needs the pane's controls to point at, and the pane's tour buttons
// need the tour: the pane is built first and handed a callback that closes the loop.
let startTour: (name: TourName) => void = () => {};
const pane = buildPane(state, request, (name) => startTour(name));
startTour = buildTour(el("tour"), el("pane"), {
  element: (name) => pane.element(name),
  applyPatch: (patch) => {
    Object.assign(state, patch);
    pane.refresh();
    request();
  },
  readState: () => state,
  samples: () => workers.samples(),
});

new ResizeObserver(request).observe(el("figure"));
window.addEventListener("resize", () => { resetReadout(); request(); });
window.addEventListener("keydown", (e) => {
  if (e.key.toLowerCase() !== "u" || e.target instanceof HTMLInputElement) return;
  state.dimensional = !state.dimensional;
  pane.refresh();
  request();
});
request();
