// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * What each pane control offers: the choices of a list, the range of a
 * slider. Kept out of `controls.ts`, which needs Tweakpane and so a DOM, so
 * that `tests/tour.test.ts` can check a tour only ever sets a control to a
 * value the pane could show.
 */

import { SUPPORTS } from "../beam/beam";
import { EDGES } from "../plate/plate";
import { MAX_TERMS } from "../random/kl";
import type { View } from "./state";
import { QOI_NAME } from "./views/convergence";
import { KERNEL_NAME } from "./views/field";
import { PLATE_MODES } from "./views/plate";

const named = <T extends string>(keys: readonly T[]) => Object.fromEntries(keys.map((k) => [k, k])) as Record<T, T>;

/** Lists: label → value. */
export const OPTIONS = {
  view: {
    "1 · spline basis": "basis",
    "2 · beam": "beam",
    "3 · convergence": "convergence",
    "3 · spectrum": "spectrum",
    "4 · random field": "field",
    "5 · Monte Carlo": "montecarlo",
    "6 · multilevel MC": "mlmc",
    "7 · plate": "plate",
    "8 · MLMC, live": "live",
  } satisfies Record<string, View>,
  structure: { beam: "beam", "plate (Kirchhoff)": "plate" },
  continuity: { "maximal Cᵖ⁻¹": "max", "C¹": "c1", "C⁰": "c0" },
  derivative: { "B": 0, "B′": 1, "B″": 2 },
  supports: named(Object.keys(SUPPORTS) as (keyof typeof SUPPORTS)[]),
  section: { uniform: "uniform", "tapered (depth → ½)": "tapered" },
  load: { uniform: "uniform", "point (tip / midspan)": "point" },
  plateLoad: { "uniform pressure": "uniform", "point (centre / free edge)": "point" },
  motion: { "a natural mode": "mode", "forced response at Ω": "response" },
  edges: named(Object.keys(EDGES) as (keyof typeof EDGES)[]),
  aspect: { "1": 1, "1.5": 1.5, "2": 2 },
  qoi: Object.fromEntries(Object.entries(QOI_NAME).map(([k, v]) => [v, k])),
  ne0: { "2": 2, "4": 4, "8": 8 },
  kernel: Object.fromEntries(Object.entries(KERNEL_NAME).map(([k, v]) => [v, k])),
  mcSamples: { "10²": 100, "10³": 1000, "10⁴": 10000, "10⁵": 100000 },
  mlSurvey: { "500": 500, "10³": 1000, "2·10³": 2000, "10⁴": 10000 },
  mlEps: { "10⁻²": 1e-2, "3·10⁻³": 3e-3, "10⁻³": 1e-3, "3·10⁻⁴": 3e-4, "10⁻⁴": 1e-4 },
} as const;

/** Sliders: inclusive range and step. */
export const RANGES = {
  p: { min: 1, max: 6, step: 1 },
  ne: { min: 1, max: 64, step: 1 },
  mode: { min: 1, max: 8, step: 1 },
  forceRatio: { min: 0.05, max: 3, step: 0.01 },
  zeta: { min: 0.002, max: 0.3, step: 0.001 },
  nu: { min: 0, max: 0.49, step: 0.01 },
  plateNe: { min: 1, max: 32, step: 1 },
  plateMode: { min: 1, max: PLATE_MODES, step: 1 },
  levels: { min: 2, max: 9, step: 1 },
  ell: { min: 0.02, max: 2, step: 0.01 },
  sigma: { min: 0, max: 1.5, step: 0.01 },
  terms: { min: 1, max: MAX_TERMS, step: 1 },
  loadSigma: { min: 0, max: 1, step: 0.01 },
  seed: { min: 1, max: 99999, step: 1 },
  mcLevel: { min: 0, max: 6, step: 1 },
  liveLevel: { min: 1, max: 8, step: 1 },
} as const;
