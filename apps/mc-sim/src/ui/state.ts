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
import type { EdgeName } from "../plate/plate";
import type { PlateCase } from "../plate/qoi";
import type { FieldSpec } from "../random/field";
import type { Kernel } from "../random/kl";
import type { Display, ReferenceBeam, ReferencePlate } from "./units";

export type View = "intro" | "basis" | "structures" | "convergence" | "field" | "montecarlo" | "mlmc" | "live";
/** What the beam and plate view animates: a natural mode, or the steady forced response at Ω. */
export type Motion = "mode" | "response";
export type Structure = "beam" | "plate";
export type Continuity = "max" | "c1" | "c0";

export interface State {
  view: View;
  dimensional: boolean;
  /** What the hierarchy, field and sampling views work on; the beam and plate view shows both. */
  structure: Structure;
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
  motion: Motion;
  /** Draw the element boundaries over the solution. */
  mesh: boolean;
  // plate
  edges: EdgeName;
  /** Side along x over side along y. */
  aspect: number;
  nu: number;
  /** Plate meshes: elements per side. */
  plateNe: number;
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
  /** The live view: the level whose coupled pairs it shows. */
  liveLevel: number;
  /** The live view's tolerance: which of the sweep's, 0 the coarsest (16 ε_min) to 4 the finest. */
  cmpTol: number;
  // reference beam, in pane units
  L: number;
  E_GPa: number;
  b_mm: number;
  h_mm: number;
  rho: number;
  q0: number;
  P0: number;
  // reference plate, in pane units (E, ρ, P₀ and the side L shared with the beam)
  t_mm: number;
  qPa: number;
}

export const defaultState = (): State => ({
  view: "intro",
  dimensional: false,
  structure: "beam",
  p: 3,
  continuity: "max",
  ne: 6,
  derivative: 0,
  supports: "cantilever",
  section: "uniform",
  load: "uniform",
  mode: 1,
  polygon: false,
  motion: "mode",
  mesh: false,
  edges: "SSSS",
  aspect: 1,
  nu: 0.3,
  plateNe: 8,
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
  mlSurvey: 200,
  mlEps: 3e-4,
  liveLevel: 2,
  cmpTol: 4,
  L: 1,
  E_GPa: 210,
  b_mm: 20,
  h_mm: 20,
  rho: 7850,
  q0: 100,
  P0: 100,
  t_mm: 5,
  qPa: 1000,
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

export function plateReferenceOf(s: State): ReferencePlate {
  return { L: s.L, E: s.E_GPa * 1e9, nu: s.nu, t: s.t_mm / 1e3, rho: s.rho, q0: s.qPa, P0: s.P0 };
}

/** How numbers are shown: for the pane's structure, or the one a view is about. */
export const displayOf = (s: State, structure: Structure = s.structure): Display =>
  ({ dimensional: s.dimensional, ref: structure === "plate" ? plateReferenceOf(s) : referenceOf(s) });

export const plateCaseOf = (s: State): PlateCase => ({
  edges: s.edges, load: s.load, aspect: s.aspect, nu: s.nu, forcing: { ratio: s.forceRatio, zeta: s.zeta },
});

export const beamCaseOf = (s: State): BeamCase => ({
  supports: s.supports, load: s.load, section: s.section, forcing: { ratio: s.forceRatio, zeta: s.zeta },
});

export const fieldSpecOf = (s: State): FieldSpec => ({
  kernel: s.kernel, ell: s.ell, sigma: s.sigma, terms: s.terms, massFollows: s.massFollows, loadSigma: s.loadSigma,
});
