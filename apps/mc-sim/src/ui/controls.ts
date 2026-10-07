// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The Tweakpane pane. It owns no solver state: every control writes into the
 * one `State` object and calls `onChange`, and each view shows only the
 * controls it reads.
 */

import type { BindingApi } from "@tweakpane/core";
import { Pane, type ButtonApi, type FolderApi } from "tweakpane";
import { SUPPORTS } from "../beam/beam";
import { QOI_NAME } from "./views/convergence";
import type { State, View } from "./state";

const VIEWS: Record<string, View> = {
  "1 · spline basis": "basis",
  "2 · beam": "beam",
  "3 · convergence": "convergence",
  "3 · spectrum": "spectrum",
};

export interface PaneHandle {
  /** Re-sync labels and visibility after `State` changed from outside the pane. */
  refresh(): void;
}

export function buildPane(st: State, onChange: () => void): PaneHandle {
  const pane = new Pane({ title: "MC-sim", container: document.getElementById("pane") ?? undefined });
  const changed = () => { sync(); onChange(); };

  pane.addBinding(st, "view", { label: "view", options: VIEWS }).on("change", changed);
  const units: ButtonApi = pane.addButton({ title: "" });
  units.on("click", () => { st.dimensional = !st.dimensional; changed(); });

  // ---- discretisation ----
  const disc = pane.addFolder({ title: "discretisation" });
  const p = disc.addBinding(st, "p", { label: "degree p", min: 1, max: 6, step: 1 });
  const cont = disc.addBinding(st, "continuity", {
    label: "continuity",
    options: { "maximal Cᵖ⁻¹": "max", "C¹": "c1", "C⁰": "c0" },
  });
  const ne = disc.addBinding(st, "ne", { label: "elements", min: 1, max: 64, step: 1 });
  const der = disc.addBinding(st, "derivative", { label: "derivative", options: { "B": 0, "B′": 1, "B″": 2 } });
  for (const b of [p, cont, ne, der] as BindingApi[]) b.on("change", changed);

  // ---- beam ----
  const beam = pane.addFolder({ title: "beam" });
  const sup = beam.addBinding(st, "supports", {
    label: "supports", options: Object.fromEntries(Object.keys(SUPPORTS).map((k) => [k, k])),
  });
  const sec = beam.addBinding(st, "section", { label: "section", options: { uniform: "uniform", "tapered (depth → ½)": "tapered" } });
  const show = beam.addBinding(st, "show", { label: "show", options: { "vibration modes": "modes", "static deflection": "deflection" } });
  const load = beam.addBinding(st, "load", { label: "load", options: { uniform: "uniform", "point (tip / midspan)": "point" } });
  const mode = beam.addBinding(st, "mode", { label: "mode", min: 1, max: 8, step: 1 });
  const poly = beam.addBinding(st, "polygon", { label: "control polygon" });
  for (const b of [sup, sec, show, load, mode, poly] as BindingApi[]) b.on("change", changed);

  // ---- hierarchy ----
  const hier = pane.addFolder({ title: "hierarchy" });
  const qoi = hier.addBinding(st, "qoi", {
    label: "quantity", options: Object.fromEntries(Object.entries(QOI_NAME).map(([k, v]) => [v, k])),
  });
  const ne0 = hier.addBinding(st, "ne0", { label: "coarsest ne₀", options: { "2": 2, "4": 4, "8": 8 } });
  const levels = hier.addBinding(st, "levels", { label: "levels", min: 2, max: 9, step: 1 });
  for (const b of [qoi, ne0, levels] as BindingApi[]) b.on("change", changed);

  // ---- reference beam (dimensional display only) ----
  const ref: FolderApi = pane.addFolder({ title: "reference beam (units)", expanded: false });
  const refs = [
    ref.addBinding(st, "L", { label: "span L [m]", min: 0.05, max: 50 }),
    ref.addBinding(st, "E_GPa", { label: "E [GPa]", min: 0.01, max: 1000 }),
    ref.addBinding(st, "b_mm", { label: "width b [mm]", min: 0.1, max: 5000 }),
    ref.addBinding(st, "h_mm", { label: "depth h [mm]", min: 0.1, max: 5000 }),
    ref.addBinding(st, "rho", { label: "ρ [kg/m³]", min: 1, max: 25000 }),
    ref.addBinding(st, "q0", { label: "q₀ [N/m]", min: 0, max: 1e6 }),
    ref.addBinding(st, "P0", { label: "P₀ [N]", min: 0, max: 1e6 }),
  ];
  for (const b of refs) b.on("change", onChange);

  function sync(): void {
    const v = st.view;
    units.title = st.dimensional ? "units: dimensional  ⇄" : "units: nondimensional  ⇄";
    const beamish = v === "beam" || v === "convergence" || v === "spectrum";
    der.hidden = v !== "basis";
    ne.hidden = v === "convergence";
    beam.hidden = !beamish;
    show.hidden = mode.hidden = poly.hidden = v !== "beam";
    mode.hidden ||= st.show !== "modes";
    load.hidden = v === "spectrum" || (v === "beam" && st.show === "modes");
    hier.hidden = v !== "convergence";
    ref.hidden = !st.dimensional;
    units.hidden = v === "basis";
  }
  sync();
  return { refresh: () => { pane.refresh(); sync(); } };
}
