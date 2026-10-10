// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * What the hierarchy, field and sampling views need to know about the
 * structure the pane has chosen — beam or plate — in one place: whether a
 * mesh can carry it, what a solve costs, the random field it reads, and the
 * spec a Monte Carlo run is built from.
 */

import { SUPPORTS, admissible as beamAdmissible, freeDofs } from "../../beam/beam";
import type { McSpec } from "../../mc/sampler";
import { admissible as plateAdmissible, solveWork } from "../../plate/plate";
import { plateQoiPoint } from "../../plate/qoi";
import { productTerms } from "../../random/field2d";
import { MAX_TERMS, karhunenLoeve, termsOf, truncate, type KLData } from "../../random/kl";
import { beamCaseOf, continuityK, fieldSpecOf, plateCaseOf, type State } from "../state";
import { qoiPoint } from "../../hierarchy";

/**
 * The finest plate mesh the hierarchy and sampling views solve, per side. A
 * cubic plate on 32 × 32 is a 1,000-unknown banded solve, and a sample there
 * with its parent takes ~70 ms; on 64 × 64 it is ~0.5 s, and a 2,000-sample
 * MLMC survey of that level alone would take minutes on eight workers.
 */
export const PLATE_MAX_NE = 32;

export const isPlate = (st: State) => st.structure === "plate";

/** Why the pane's structure cannot be solved on `ne` elements (per side), or null. */
export function admissibleAt(st: State, ne: number, p = st.p): string | null {
  const k = continuityK(st.continuity, p);
  return isPlate(st)
    ? plateAdmissible({ p, k, ne, edges: st.edges, aspect: st.aspect, nu: st.nu })
    : beamAdmissible({ p, k, ne, supports: SUPPORTS[st.supports] });
}

/**
 * Cost of one solve on `ne` elements in the unit MLMC allocates by: free
 * coefficients for a beam (its banded solve is O(dofs·p²), p fixed), and
 * dofs × bandwidth² — the banded factorisation's work — for a plate.
 */
export function workAt(st: State, ne: number): number {
  const k = continuityK(st.continuity, st.p);
  return isPlate(st)
    ? solveWork({ p: st.p, k, ne, edges: st.edges, aspect: st.aspect, nu: st.nu })
    : freeDofs({ p: st.p, k, ne, supports: SUPPORTS[st.supports] });
}

export const workUnit = (st: State) => (isPlate(st) ? "dofs × bandwidth²" : "free coefficients");

/** Levels a plate run may use from ne₀ without passing `max` elements per side. */
export function levelsWithin(ne0: number, levels: number, max: number): number {
  let L = levels;
  while (L > 1 && ne0 * 2 ** (L - 1) > max) L--;
  return L;
}

const cache = new Map<string, KLData>();

/** The full expansion for a kernel and length on [0, 1], solved once and kept. */
function fullKL(kernel: State["kernel"], ell: number): KLData {
  const key = `${kernel}/${ell}`;
  let kl = cache.get(key);
  if (!kl) {
    if (cache.size > 8) cache.delete(cache.keys().next().value!);
    cache.set(key, (kl = karhunenLoeve(kernel, ell, MAX_TERMS)));
  }
  return kl;
}

/**
 * The expansion(s) a run reads: along the beam; or along the plate's x (for
 * length ℓ/a, read at x/a) and y. A plate keeps M products of the two, and
 * none needs more than M terms of either factor.
 */
export function klsOf(st: State): { kl: KLData; klY?: KLData; M: number } {
  if (!isPlate(st)) {
    const kl = truncate(fullKL(st.kernel, st.ell), st.terms);
    return { kl, M: termsOf(kl) };
  }
  const kl = truncate(fullKL(st.kernel, st.ell / st.aspect), st.terms);
  const klY = st.aspect === 1 ? kl : truncate(fullKL(st.kernel, st.ell), st.terms);
  return { kl, klY: st.aspect === 1 ? undefined : klY, M: productTerms(kl, klY, st.terms).values.length };
}

export function mcSpecOf(st: State): McSpec {
  const k = continuityK(st.continuity, st.p);
  return {
    p: st.p, k, beam: beamCaseOf(st), plate: isPlate(st) ? plateCaseOf(st) : undefined,
    qoi: st.qoi, field: fieldSpecOf(st), seed: st.seed,
  };
}

/** Where the point load acts and the deflection is read, as text. */
export function pointText(st: State): string {
  if (!isPlate(st)) return `x = ${qoiPoint(st.supports)}`;
  const at = plateQoiPoint(st);
  return `(x, y) = (${+at.x.toFixed(3)}, ${at.y})`;
}

/** "beam, cantilever" or "plate SSSS, 1 × 1": the structure in a few words. */
export function structureText(st: State): string {
  return isPlate(st) ? `${st.edges} plate, ${+st.aspect.toFixed(3)} × 1, ν = ${st.nu}` : `${st.supports} beam`;
}
