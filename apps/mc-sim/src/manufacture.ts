// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The introduction's factory: "manufacture" one beam or plate — draw its
 * random stiffness, solve its lowest modes — and hand back what a picture of
 * it needs. The same field, solver and Philox draws as every other view, at
 * deliberately exaggerated settings: σ of log EI is 0.5 (the tours use 0.3),
 * and the section's depth follows the stiffness (EI ∝ h³, so h ∝ e^{1/3}, and
 * the mass with it), so what varies is visibly the part's thickness.
 *
 * Part i of a kind is a pure function of i, as a Monte Carlo sample is: the
 * beams read Philox stream 0 and the plates stream 1, and part −1 is the
 * design — the uniform structure every made one is compared with.
 *
 * The beam is a cantilever. The plate is free on all four edges — Chladni's
 * plate — because its second and third modes (the "ring" and the "X") lie
 * close together, and a random stiffness mixes them: the nodal lines of a made
 * plate wander far from the design's, which is the point of the picture. Its
 * three rigid modes, at ω = 0, are skipped.
 *
 * Pure: no DOM; `ui/views/intro.ts` draws it.
 */

import { Beam, SUPPORTS } from "./beam/beam";
import { gaussLegendre } from "./quad";
import { Plate, plateTables, type PlateTables } from "./plate/plate";
import { FieldAt, draw, realise, type FieldSpec, type PointField } from "./random/field";
import { FieldOnGrid } from "./random/field2d";
import { MAX_TERMS, karhunenLoeve, termsOf, truncate, type KLData } from "./random/kl";
import { splineSpace, tabulate, type Tabulation } from "./spline";

export const MADE_FIELD: FieldSpec = {
  kernel: "matern32", ell: 0.3, sigma: 0.5, terms: 24, massFollows: true, loadSigma: 0,
};
export const MADE_SEED = 2026;
/** Modes offered for each part. */
export const MADE_MODES = 3;
/** Points along the beam, and per side of the plate, that a part is drawn at. */
export const BEAM_POINTS = 81;
export const PLATE_POINTS = 31;

const P = 3, BEAM_NE = 24, PLATE_NE = 14, RIGID = 3;

export type Kind = "beam" | "plate";

export interface Part {
  readonly kind: Kind;
  /** −1 for the design. */
  readonly index: number;
  /** Depth relative to the design's, (e/e₀)^{1/3}, at the drawing points (flat [iy·n + ix] on the plate). */
  readonly depth: Float64Array;
  /** ω of the lowest elastic modes, nondimensional. */
  readonly omega: Float64Array;
  /** Each mode at the drawing points: peak |φ| = 1, signed to agree with the design's. */
  readonly shapes: readonly Float64Array[];
  /** Wall time of the draw, assembly and eigensolve. */
  readonly ms: number;
}

const linspace = (n: number) => Float64Array.from({ length: n }, (_, i) => i / (n - 1));
export const beamX = linspace(BEAM_POINTS);
export const plateX = linspace(PLATE_POINTS);

interface Shop {
  readonly M: number;
  /** The field at the quadrature points (the solve), and at the drawing points (the picture). */
  readonly solveAt: PointField;
  readonly drawAt: PointField;
  readonly tab: Tabulation | PlateTables;
  readonly ones: Float64Array;
}

let kl: KLData | null = null;
const shops = new Map<Kind, Shop>();
const designs = new Map<Kind, Part>();

function shop(kind: Kind): Shop {
  let s = shops.get(kind);
  if (s) return s;
  kl ??= truncate(karhunenLoeve(MADE_FIELD.kernel, MADE_FIELD.ell, MAX_TERMS), MADE_FIELD.terms);
  if (kind === "beam") {
    const tab = tabulate(splineSpace(P, BEAM_NE), gaussLegendre(P + 1), 2);
    s = { M: termsOf(kl), solveAt: new FieldAt(kl, tab.x), drawAt: new FieldAt(kl, beamX), tab, ones: new Float64Array(tab.x.length).fill(1) };
  } else {
    const tab = plateTables(P, PLATE_NE);
    const solveAt = new FieldOnGrid(kl, kl, MADE_FIELD.terms, tab.tx.x, tab.ty.x);
    s = { M: solveAt.M, solveAt, drawAt: new FieldOnGrid(kl, kl, MADE_FIELD.terms, plateX, plateX), tab, ones: new Float64Array(solveAt.n).fill(1) };
  }
  shops.set(kind, s);
  return s;
}

/** The design: the uniform part, every made one's reference. */
export function design(kind: Kind): Part {
  let d = designs.get(kind);
  if (!d) designs.set(kind, (d = make(kind, -1)));
  return d;
}

/** Part `index` of a kind, or the design for −1. */
export function make(kind: Kind, index: number): Part {
  const t0 = performance.now(), s = shop(kind);
  let e = s.ones, mu = s.ones, depth = new Float64Array(s.drawAt.n).fill(1);
  if (index >= 0) {
    const d = draw(MADE_FIELD, s.M, MADE_SEED, index, kind === "beam" ? 0 : 1);
    ({ e, mu } = realise(s.solveAt, MADE_FIELD, d, s.ones, s.ones));
    depth = s.drawAt.lognormal(MADE_FIELD.sigma, d.xi).map(Math.cbrt);
  }
  let values: Float64Array, shapes: Float64Array[];
  if (kind === "beam") {
    const beam = new Beam({ p: P, ne: BEAM_NE, supports: SUPPORTS.cantilever, stiffness: e, mass: mu }, s.tab as Tabulation);
    const m = beam.modes(MADE_MODES);
    values = m.values;
    shapes = m.vectors.map((c) => beamX.map((x) => beam.evaluate(c, x)[0]));
  } else {
    const plate = new Plate({ p: P, ne: PLATE_NE, edges: "FFFF", stiffness: e, mass: mu }, s.tab as PlateTables);
    const m = plate.modes(MADE_MODES + RIGID);
    values = m.values.slice(RIGID);
    shapes = m.vectors.slice(RIGID).map((c) => plate.grid(c, plateX, plateX));
  }
  const ref = index >= 0 ? design(kind).shapes : null;
  shapes.forEach((f, n) => normalise(f, ref?.[n]));
  return { kind, index, depth, omega: values.map((v) => Math.sqrt(Math.max(v, 0))), shapes, ms: performance.now() - t0 };
}

/** Scale to peak |f| = 1, signed to agree with `ref` (or, for the design, with its largest value positive). */
function normalise(f: Float64Array, ref?: Float64Array): void {
  let peak = 0;
  for (const v of f) if (Math.abs(v) > Math.abs(peak)) peak = v;
  let sign = Math.sign(peak) || 1;
  if (ref) {
    let dot = 0;
    for (let i = 0; i < f.length; i++) dot += f[i] * ref[i];
    sign = dot < 0 ? -1 : 1;
  }
  const k = sign / Math.abs(peak || 1);
  for (let i = 0; i < f.length; i++) f[i] *= k;
}
