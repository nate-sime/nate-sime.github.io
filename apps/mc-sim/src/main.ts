// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Entry point: one canvas, one readout, one pane, and a view per stage of
 * MC_PLAN.md that has landed —
 *
 *   1  the spline basis           (`ui/views/basis.ts`)
 *   2  the beam, static and modal (`ui/views/beam.ts`)
 *   3  the discretisation hierarchy, and the spectrum
 *                                 (`ui/views/convergence.ts`, `ui/views/spectrum.ts`)
 *   4  the random input           (`ui/views/field.ts`)
 *   5  plain Monte Carlo          (`ui/views/montecarlo.ts`, sampling in `mc/`)
 *
 * Everything is computed on the CPU in f64 and drawn into a 2-D canvas; Monte
 * Carlo samples run in a pool of Web Workers, paused while their view is not
 * on screen.
 * Frames are requested only while something moves (a vibrating mode);
 * otherwise a redraw follows a control change, a resize, or the pointer.
 * `U` toggles dimensional units, as the pane's button does.
 */

import { buildPane } from "./ui/controls";
import { Plot } from "./ui/plot";
import { defaultState, type View } from "./ui/state";
import { renderBasis } from "./ui/views/basis";
import { renderBeam } from "./ui/views/beam";
import { renderConvergence } from "./ui/views/convergence";
import { renderField } from "./ui/views/field";
import { monteCarlo, renderMonteCarlo } from "./ui/views/montecarlo";
import { renderSpectrum } from "./ui/views/spectrum";
import type { ViewResult } from "./ui/views/view";

const el = (id: string) => document.getElementById(id)!;

el("caption").textContent = `MC-sim ${__APP_VERSION__} · Monte Carlo for vibrating structures`;

const state = defaultState();
const canvas = el("view") as HTMLCanvasElement;
const plot = new Plot(canvas);
const readout = el("readout");

const RENDER: Record<View, (t: number) => ViewResult> = {
  basis: () => renderBasis(plot, state),
  beam: (t) => renderBeam(plot, state, t),
  convergence: () => renderConvergence(plot, state),
  spectrum: () => renderSpectrum(plot, state),
  field: () => renderField(plot, state),
  montecarlo: () => renderMonteCarlo(plot, state),
};

let frame = 0;
let shown: View | null = null;
function render(t = performance.now()): void {
  frame = 0;
  plot.resize();
  if (state.view !== shown) {
    if (shown === "montecarlo") monteCarlo.pause();
    else if (state.view === "montecarlo") monteCarlo.resume();
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
monteCarlo.connect(request);
const pane = buildPane(state, request);

new ResizeObserver(request).observe(canvas);
window.addEventListener("keydown", (e) => {
  if (e.key.toLowerCase() !== "u" || e.target instanceof HTMLInputElement) return;
  state.dimensional = !state.dimensional;
  pane.refresh();
  request();
});
request();
