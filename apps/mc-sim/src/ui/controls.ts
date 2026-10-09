// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The Tweakpane pane. It owns no solver state: every control writes into the
 * one `State` object and calls `onChange`, and each view shows only the
 * controls it reads — which ones is `visibility.ts`'s pure rule, applied here.
 *
 * The first folder starts the guided tours (`tour.ts`), which point at the
 * controls below by name: `PaneHandle.element` is how they find them, since
 * Tweakpane renders no IDs of its own.
 */

import type { BindingApi } from "@tweakpane/core";
import { Pane, type ButtonApi, type FolderApi } from "tweakpane";
import { OPTIONS, RANGES } from "./options";
import { workers } from "./views/workers";
import type { State } from "./state";
import { TOUR_NAMES, type TourName } from "./tours";
import { CONTROLS, hidden, type ControlName } from "./visibility";

export interface PaneHandle {
  /** Re-sync labels and visibility after `State` changed from outside the pane. */
  refresh(): void;
  /** The element of a control, for the tour to light. */
  element(name: ControlName): HTMLElement;
}

export function buildPane(st: State, onChange: () => void, onTour: (name: TourName) => void): PaneHandle {
  const pane = new Pane({ title: "MC-sim", container: document.getElementById("pane") ?? undefined });
  const changed = () => { sync(); onChange(); };

  // ---- guided tours ----
  const tours = pane.addFolder({ title: "guided tours" });
  for (const name of TOUR_NAMES) tours.addButton({ title: name }).on("click", () => onTour(name));

  const view = pane.addBinding(st, "view", { label: "view", options: OPTIONS.view });
  view.on("change", changed);
  const structure = pane.addBinding(st, "structure", { label: "structure", options: OPTIONS.structure });
  structure.on("change", changed);
  const units: ButtonApi = pane.addButton({ title: "" });
  units.on("click", () => { st.dimensional = !st.dimensional; changed(); });

  // ---- discretisation ----
  const disc = pane.addFolder({ title: "discretisation" });
  const p = disc.addBinding(st, "p", { label: "degree p", ...RANGES.p });
  const cont = disc.addBinding(st, "continuity", { label: "continuity", options: OPTIONS.continuity });
  const ne = disc.addBinding(st, "ne", { label: "elements", ...RANGES.ne });
  const der = disc.addBinding(st, "derivative", { label: "derivative", options: OPTIONS.derivative });
  for (const b of [p, cont, ne, der] as BindingApi[]) b.on("change", changed);

  // ---- beam ----
  const beam = pane.addFolder({ title: "beam" });
  const sup = beam.addBinding(st, "supports", { label: "supports", options: OPTIONS.supports });
  const sec = beam.addBinding(st, "section", { label: "section", options: OPTIONS.section });
  const load = beam.addBinding(st, "load", { label: "load", options: OPTIONS.load });
  const mode = beam.addBinding(st, "mode", { label: "mode", ...RANGES.mode });
  const poly = beam.addBinding(st, "polygon", { label: "control polygon" });
  for (const b of [sup, sec, load, mode, poly] as BindingApi[]) b.on("change", changed);

  // ---- what the beam and plate view animates ----
  const motion = pane.addFolder({ title: "motion" });
  const mo = {
    motion: motion.addBinding(st, "motion", { label: "show", options: OPTIONS.motion }),
    force: motion.addBinding(st, "forceRatio", { label: "forcing Ω/ω₁", ...RANGES.forceRatio }),
    zeta: motion.addBinding(st, "zeta", { label: "damping ζ", ...RANGES.zeta }),
    mesh: motion.addBinding(st, "mesh", { label: "mesh" }),
  };
  // Forcing and damping are shared with the hierarchy folder: a change here re-reads both.
  for (const b of Object.values(mo) as BindingApi[]) b.on("change", () => { pane.refresh(); changed(); });

  // ---- plate ----
  const plate = pane.addFolder({ title: "plate" });
  const pb = {
    edges: plate.addBinding(st, "edges", { label: "edges (x=0, y=0, x=a, y=1)", options: OPTIONS.edges }),
    aspect: plate.addBinding(st, "aspect", { label: "aspect a", options: OPTIONS.aspect }),
    nu: plate.addBinding(st, "nu", { label: "Poisson ν", ...RANGES.nu }),
    load: plate.addBinding(st, "load", { label: "load", options: OPTIONS.plateLoad }),
    ne: plate.addBinding(st, "plateNe", { label: "elements / side", ...RANGES.plateNe }),
  };
  // Load is shared with the beam folder: a change here re-reads both.
  for (const b of Object.values(pb) as BindingApi[]) b.on("change", () => { pane.refresh(); changed(); });

  // ---- hierarchy ----
  const hier = pane.addFolder({ title: "hierarchy" });
  const qoi = hier.addBinding(st, "qoi", { label: "quantity", options: OPTIONS.qoi });
  const ne0 = hier.addBinding(st, "ne0", { label: "coarsest ne₀", options: OPTIONS.ne0 });
  const levels = hier.addBinding(st, "levels", { label: "levels", ...RANGES.levels });
  const force = hier.addBinding(st, "forceRatio", { label: "forcing Ω/ω₁", ...RANGES.forceRatio });
  const zeta = hier.addBinding(st, "zeta", { label: "damping ζ", ...RANGES.zeta });
  for (const b of [qoi, ne0, levels, force, zeta] as BindingApi[]) b.on("change", changed);

  // ---- random field ----
  const rf = pane.addFolder({ title: "random field" });
  const rfb: BindingApi[] = [
    rf.addBinding(st, "kernel", { label: "kernel", options: OPTIONS.kernel }),
    rf.addBinding(st, "ell", { label: "corr. length ℓ/L", ...RANGES.ell }),
    rf.addBinding(st, "sigma", { label: "σ of log EI", ...RANGES.sigma }),
    rf.addBinding(st, "terms", { label: "KL terms M", ...RANGES.terms }),
    rf.addBinding(st, "massFollows", { label: "mass follows depth" }),
    rf.addBinding(st, "loadSigma", { label: "load σ_q", ...RANGES.loadSigma }),
    rf.addBinding(st, "seed", { label: "seed", ...RANGES.seed }),
  ];
  for (const b of rfb as BindingApi[]) b.on("change", changed);

  // ---- Monte Carlo ----
  const mc = pane.addFolder({ title: "Monte Carlo" });
  const mcb: BindingApi[] = [
    mc.addBinding(st, "mcLevel", { label: "level ℓ", ...RANGES.mcLevel }),
    mc.addBinding(st, "mcSamples", { label: "samples N", options: OPTIONS.mcSamples }),
  ];
  for (const b of mcb as BindingApi[]) b.on("change", changed);
  mc.addButton({ title: "pause / resume" }).on("click", () => workers.toggle());
  mc.addButton({ title: "next seed" }).on("click", () => { st.seed++; pane.refresh(); changed(); });

  // ---- multilevel Monte Carlo ----
  const ml = pane.addFolder({ title: "multilevel Monte Carlo" });
  const mlb: BindingApi[] = [
    ml.addBinding(st, "mlSurvey", { label: "survey N / level", options: OPTIONS.mlSurvey }),
    ml.addBinding(st, "mlEps", { label: "finest ε (rel.)", options: OPTIONS.mlEps }),
  ];
  for (const b of mlb as BindingApi[]) b.on("change", changed);
  const liveLevel = ml.addBinding(st, "liveLevel", { label: "level shown ℓ", ...RANGES.liveLevel });
  liveLevel.on("change", changed);
  const cmpTol = ml.addBinding(st, "cmpTol", { label: "tolerance shown", options: OPTIONS.cmpTol });
  cmpTol.on("change", changed);
  // The same samples again — a run is a function of its seed — but arriving, to be watched.
  const restart = ml.addButton({ title: "run again from zero" });
  restart.on("click", () => workers.restart());
  ml.addButton({ title: "pause / resume" }).on("click", () => workers.toggle());
  ml.addButton({ title: "next seed" }).on("click", () => { st.seed++; pane.refresh(); changed(); });

  // ---- reference beam (dimensional display only) ----
  const ref: FolderApi = pane.addFolder({ title: "reference structure (units)", expanded: false });
  const refs = [
    ref.addBinding(st, "L", { label: "span L [m]", min: 0.05, max: 50 }),
    ref.addBinding(st, "E_GPa", { label: "E [GPa]", min: 0.01, max: 1000 }),
    ref.addBinding(st, "b_mm", { label: "width b [mm]", min: 0.1, max: 5000 }),
    ref.addBinding(st, "h_mm", { label: "depth h [mm]", min: 0.1, max: 5000 }),
    ref.addBinding(st, "rho", { label: "ρ [kg/m³]", min: 1, max: 25000 }),
    ref.addBinding(st, "q0", { label: "q₀ [N/m]", min: 0, max: 1e6 }),
    ref.addBinding(st, "P0", { label: "P₀ [N]", min: 0, max: 1e6 }),
  ];
  const refPlate = [
    ref.addBinding(st, "t_mm", { label: "plate t [mm]", min: 0.1, max: 500 }),
    ref.addBinding(st, "qPa", { label: "plate q₀ [Pa]", min: 0, max: 1e7 }),
  ];
  for (const b of refPlate) b.on("change", onChange);
  for (const b of refs) b.on("change", onChange);

  // Every control by name: what `visibility.ts` rules on and the tour points at.
  const controls: Record<ControlName, { hidden: boolean; element: HTMLElement }> = {
    tours, view, structure, units,
    disc, beam, motion, plate, hierarchy: hier, field: rf, mc, ml, reference: ref,
    p, continuity: cont, ne, derivative: der,
    supports: sup, section: sec, load, mode, polygon: poly,
    motionShow: mo.motion, motionForce: mo.force, motionZeta: mo.zeta, mesh: mo.mesh,
    edges: pb.edges, aspect: pb.aspect, nu: pb.nu, plateLoad: pb.load, plateNe: pb.ne,
    qoi, ne0, levels, force, zeta,
    kernel: rfb[0], ell: rfb[1], sigma: rfb[2], terms: rfb[3], massFollows: rfb[4], loadSigma: rfb[5], seed: rfb[6],
    mcLevel: mcb[0], mcSamples: mcb[1], mlSurvey: mlb[0], mlEps: mlb[1], liveLevel, cmpTol, restart,
  };

  function sync(): void {
    units.title = st.dimensional ? "units: dimensional  ⇄" : "units: nondimensional  ⇄";
    const h = hidden(st);
    for (const name of CONTROLS) controls[name].hidden = h[name];
    const onPlate = !h.plate;
    for (const b of refPlate) b.hidden = !onPlate;
    for (const i of [2, 3, 5]) refs[i].hidden = onPlate; // section b × h and q₀ per length: the beam's only
  }
  sync();
  return {
    refresh: () => { pane.refresh(); sync(); },
    element: (name) => controls[name].element,
  };
}
