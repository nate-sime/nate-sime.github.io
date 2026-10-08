// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Everything the pane can set, in one plain object. Views read it; nothing
 * else owns it. The reference beam is held in the units a person types (GPa,
 * mm) and converted once, in `referenceOf`.
 */

import type { SupportName } from "../beam/beam";
import type { BeamCase, LoadCase, QoI, Section } from "../hierarchy";
import type { FieldSpec } from "../random/field";
import type { Kernel } from "../random/kl";
import type { Display, ReferenceBeam } from "./units";

export type View = "basis" | "beam" | "convergence" | "spectrum" | "field" | "montecarlo" | "mlmc";
export type Continuity = "max" | "c1" | "c0";

export interface State {
  view: View;
  dimensional: boolean;
  // discretisation
  p: number;
  continuity: Continuity;
  ne: number;
  derivative: number;
  // beam
  supports: SupportName;
  section: Section;
  load: LoadCase;
  mode: number;
  polygon: boolean;
  // hierarchy
  qoi: QoI;
  ne0: number;
  levels: number;
  /** Forced response: Ω / ω₁ of the uniform beam, and the Rayleigh damping ratio. */
  forceRatio: number;
  zeta: number;
  // random field
  kernel: Kernel;
  ell: number;
  sigma: number;
  terms: number;
  massFollows: boolean;
  loadSigma: number;
  seed: number;
  // Monte Carlo
  mcLevel: number;
  mcSamples: number;
  // multilevel Monte Carlo
  /** Samples per level in the survey. */
  mlSurvey: number;
  /** The finest tolerance, relative to |Q| of the mean beam; the sweep runs 16, 8, 4, 2 and 1 times it. */
  mlEps: number;
  // reference beam, in pane units
  L: number;
  E_GPa: number;
  b_mm: number;
  h_mm: number;
  rho: number;
  q0: number;
  P0: number;
}

export const defaultState = (): State => ({
  view: "basis",
  dimensional: false,
  p: 3,
  continuity: "max",
  ne: 6,
  derivative: 0,
  supports: "cantilever",
  section: "uniform",
  load: "uniform",
  mode: 1,
  polygon: false,
  qoi: "omega1",
  ne0: 4,
  levels: 6,
  forceRatio: 0.8,
  zeta: 0.02,
  kernel: "matern32",
  ell: 0.2,
  sigma: 0.3,
  terms: 24,
  massFollows: false,
  loadSigma: 0,
  seed: 1,
  mcLevel: 2,
  mcSamples: 10000,
  mlSurvey: 2000,
  mlEps: 3e-4,
  L: 1,
  E_GPa: 210,
  b_mm: 20,
  h_mm: 20,
  rho: 7850,
  q0: 100,
  P0: 100,
});

/** Continuity index k for a degree, clamped to what the degree can carry. */
export function continuityK(c: Continuity, p: number): number {
  const k = c === "max" ? p - 1 : c === "c1" ? 1 : 0;
  return Math.max(0, Math.min(p - 1, k));
}

export const continuityName = (k: number) => `C${"⁰¹²³⁴⁵⁶⁷⁸⁹"[k]}`;

export function referenceOf(s: State): ReferenceBeam {
  return { L: s.L, E: s.E_GPa * 1e9, b: s.b_mm / 1e3, h: s.h_mm / 1e3, rho: s.rho, q0: s.q0, P0: s.P0 };
}

export const displayOf = (s: State): Display => ({ dimensional: s.dimensional, ref: referenceOf(s) });

export const beamCaseOf = (s: State): BeamCase => ({
  supports: s.supports, load: s.load, section: s.section, forcing: { ratio: s.forceRatio, zeta: s.zeta },
});

export const fieldSpecOf = (s: State): FieldSpec => ({
  kernel: s.kernel, ell: s.ell, sigma: s.sigma, terms: s.terms, massFollows: s.massFollows, loadSigma: s.loadSigma,
});
