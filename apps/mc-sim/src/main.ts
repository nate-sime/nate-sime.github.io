// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Entry point: one figure, one readout, one pane, and a view per stage of
 * MC_PLAN.md that has landed —
 *
 *   1  the spline basis           (`ui/views/basis.ts`)
 *   2  the beam, static and modal (`ui/views/beam.ts`)
 *   3  the discretisation hierarchy, and the spectrum
 *                                 (`ui/views/convergence.ts`, `ui/views/spectrum.ts`)
 *   4  the random input           (`ui/views/field.ts`)
 *   5  plain Monte Carlo          (`ui/views/montecarlo.ts`, sampling in `mc/`)
 *   6  multilevel Monte Carlo     (`ui/views/mlmc.ts`, the algorithm in `mc/mlmc.ts`)
 *   7  the Kirchhoff plate        (`ui/views/plate.ts`, the solver in `plate/`); the
 *                                 hierarchy, field and sampling views take it too,
 *                                 when the pane's structure is the plate
 *   8  MLMC, live                 (`ui/views/live.ts`): a run's coupled samples,
 *                                 solved and moving, beside its statistics; the
 *                                 beam and plate views animate the forced response
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
import { renderBeam } from "./ui/views/beam";
import { renderConvergence } from "./ui/views/convergence";
import { renderField } from "./ui/views/field";
import { renderMlmc } from "./ui/views/mlmc";
import { renderMonteCarlo } from "./ui/views/montecarlo";
import { renderLive } from "./ui/views/live";
import { renderPlate } from "./ui/views/plate";
import { renderSpectrum } from "./ui/views/spectrum";
import type { ViewResult } from "./ui/views/view";
import { workers } from "./ui/views/workers";

const el = (id: string) => document.getElementById(id)!;

el("caption").textContent = `MC-sim ${__APP_VERSION__} · Monte Carlo for vibrating structures`;

const state = defaultState();
const figure = new Figure(el("figure"));
const readout = el("readout");
const one = () => figure.panels(1)[0];

const RENDER: Record<View, (t: number) => ViewResult> = {
  basis: () => renderBasis(one(), state),
  beam: (t) => renderBeam(figure, state, t),
  plate: (t) => renderPlate(figure, state, t),
  convergence: () => renderConvergence(one(), state),
  spectrum: () => renderSpectrum(one(), state),
  field: () => renderField(figure, state),
  montecarlo: () => renderMonteCarlo(figure, state),
  mlmc: () => renderMlmc(figure, state),
  live: (t) => renderLive(figure, state, t),
};

let frame = 0;
let shown: View | null = null;
function render(t = performance.now()): void {
  frame = 0;
  // Leaving a sampling view idles the workers; the next view to sample opens its own run.
  if (state.view !== shown) {
    workers.idle();
    shown = state.view;
  }
  let r: ViewResult;
  try {
    r = RENDER[state.view](t);
  } catch (e) {
    r = { readout: `error: ${(e as Error).message}`, animate: false };
  }
  if (readout.textContent !== r.readout) readout.textContent = r.readout;
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
window.addEventListener("keydown", (e) => {
  if (e.key.toLowerCase() !== "u" || e.target instanceof HTMLInputElement) return;
  state.dimensional = !state.dimensional;
  pane.refresh();
  request();
});
request();
