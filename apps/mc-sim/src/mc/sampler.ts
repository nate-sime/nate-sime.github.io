// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * One Monte Carlo sample: draw ω, build the beam it describes on a level of the
 * hierarchy, solve, and read off the quantity of interest.
 *
 * Every sample on level ℓ is also solved on level ℓ − 1 with the *same* ω —
 * the same ξ, through the same KL expansion, at the coarse mesh's quadrature
 * points. The pair gives the correction Y = Q_ℓ − Q_{ℓ−1}, whose mean is how
 * far the level is from its parent (a measure of bias plain Monte Carlo cannot
 * otherwise get) and whose variance is much smaller than Q's: the two facts
 * multilevel Monte Carlo is built from.
 *
 * Sample i depends only on (spec, i) — never on which worker runs it, nor in
 * what order — so an estimate after N samples is the same number on one core
 * or eight, and any sample can be re-run alone.
 *
 * Pure: no DOM, no Worker; `worker.ts` wraps it, and the tests call it directly.
 */

import { Beam, SUPPORTS } from "../beam/beam";
import { evaluateQoI, qoiPoint, sectionOf, type BeamCase, type QoI } from "../hierarchy";
import { gaussLegendre } from "../quad";
import { FieldAt, draw, realise, type FieldSpec } from "../random/field";
import { termsOf, type KLData } from "../random/kl";
import { sampleAt, splineSpace, tabulate, type Tabulation } from "../spline";

export interface McSpec {
  readonly p: number;
  readonly k: number;
  readonly beam: BeamCase;
  readonly qoi: QoI;
  readonly field: FieldSpec;
  readonly seed: number;
}

/** Points the deflection field is recorded at, for its mean and bands. */
export const FIELD_POINTS = 65;
export const fieldX = Float64Array.from({ length: FIELD_POINTS }, (_, i) => i / (FIELD_POINTS - 1));

export interface Sample {
  readonly Q: number;
  /** Q on the parent level, same ω; NaN on the coarsest. */
  readonly Qc: number;
  /** Deflection at `fieldX`, if asked for. */
  readonly w: Float64Array | null;
}

interface LevelData {
  readonly tab: Tabulation;
  readonly field: FieldAt;
  readonly e0: Float64Array;
  readonly mu0: Float64Array;
}

export class Sampler {
  readonly M: number;
  private readonly levels = new Map<number, LevelData>();

  constructor(readonly spec: McSpec, readonly kl: KLData) {
    this.M = termsOf(kl);
  }

  /** Tables for a level: built on first use, then shared by every sample on it. */
  private level(ne: number): LevelData {
    let L = this.levels.get(ne);
    if (!L) {
      const { p, k, beam } = this.spec;
      const tab = tabulate(splineSpace(p, ne, k), gaussLegendre(p + 1), 2);
      const sec = sectionOf(beam.section);
      L = {
        tab,
        field: new FieldAt(this.kl, tab.x),
        e0: sampleAt(tab, sec.stiffness ?? (() => 1)),
        mu0: sampleAt(tab, sec.mass ?? (() => 1)),
      };
      this.levels.set(ne, L);
    }
    return L;
  }

  /** Q for sample `index` on the mesh of `ne` elements; the static solution too if `wantW`. */
  solve(index: number, ne: number, wantW = false): { Q: number; w: Float64Array | null } {
    const { spec } = this, L = this.level(ne);
    const r = realise(L.field, spec.field, draw(spec.field, this.M, spec.seed, index), L.e0, L.mu0);
    const beam = new Beam({ p: spec.p, k: spec.k, ne, supports: SUPPORTS[spec.beam.supports], stiffness: r.e, mass: r.mu }, L.tab);
    const xq = qoiPoint(spec.beam.supports);
    const load = spec.beam.load === "uniform" ? { q: r.q ?? (() => 1) } : { forces: [{ x: xq, P: r.P }] };
    const { Q, c } = evaluateQoI(beam, spec.qoi, load, xq);
    let w: Float64Array | null = null;
    if (wantW) {
      const cs = c.length ? c : beam.solve(load);
      w = fieldX.map((x) => beam.evaluate(cs, x)[0]);
    }
    return { Q, w };
  }

  /** Sample `index` on `ne` elements and, if `coarse`, on ne/2 with the same ω. */
  sample(index: number, ne: number, coarse: boolean, wantW = false): Sample {
    const fine = this.solve(index, ne, wantW);
    const Qc = coarse ? this.solve(index, ne / 2).Q : NaN;
    return { Q: fine.Q, Qc, w: fine.w };
  }
}
