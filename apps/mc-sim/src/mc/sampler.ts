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
 * Sample i depends only on (spec, stream, i) — never on which worker runs it,
 * nor in what order — so an estimate after N samples is the same number on one
 * core or eight, and any sample can be re-run alone. Plain Monte Carlo reads
 * stream 0; multilevel Monte Carlo gives each level a stream of its own.
 *
 * The structure is the beam, or — when the spec names one — a plate, whose
 * samples read a separable two-dimensional field (`random/field2d.ts`) at the
 * tensor grid of its quadrature points. Everything else is shared.
 *
 * Pure: no DOM, no Worker; `worker.ts` wraps it, and the tests call it directly.
 */

import { Beam, SUPPORTS } from "../beam/beam";
import { evaluateQoI, harmonicOf, qoiPoint, sectionOf, type BeamCase, type Harmonic, type QoI } from "../hierarchy";
import { Plate, plateTables, type PlateTables } from "../plate/plate";
import { evaluatePlateQoI, plateHarmonicOf, plateLoadOf, plateQoiPoint, type PlateCase } from "../plate/qoi";
import { gaussLegendre } from "../quad";
import { FieldAt, draw, realise, type FieldSpec, type PointField } from "../random/field";
import { FieldOnGrid } from "../random/field2d";
import { termsOf, type KLData } from "../random/kl";
import { sampleAt, splineSpace, tabulate, type Tabulation } from "../spline";

export interface McSpec {
  readonly p: number;
  readonly k: number;
  readonly beam: BeamCase;
  /** If given, the structure is this plate, not the beam: ne is then elements per side. */
  readonly plate?: PlateCase;
  readonly qoi: QoI;
  readonly field: FieldSpec;
  readonly seed: number;
}

/**
 * Points the deflection field is recorded at, for its mean and bands: along
 * the beam, or along the plate's midline y = ½ (scaled by its aspect ratio).
 */
export const FIELD_POINTS = 65;
export const fieldX = Float64Array.from({ length: FIELD_POINTS }, (_, i) => i / (FIELD_POINTS - 1));

export interface Sample {
  readonly Q: number;
  /** Q on the parent level, same ω; NaN on the coarsest. */
  readonly Qc: number;
  /** Deflection at `fieldX` (the amplitude |u|, for the forced response), if asked for. */
  readonly w: Float64Array | null;
}

/**
 * One sample on one level, solved and kept for drawing: the structure, Q, and
 * the shapes a picture of it needs. Built by the same code the workers run, so
 * its Q is theirs to the last bit.
 */
export interface Inspected {
  readonly Q: number;
  readonly beam?: Beam;
  readonly plate?: Plate;
  /** Static deflection, full coefficients. */
  readonly w: Float64Array;
  /** The forced response, real and imaginary parts, when the quantity is the response. */
  readonly u?: { readonly re: Float64Array; readonly im: Float64Array };
  /** The first mode, M-normalised and signed so its largest coefficient is positive, when the quantity is ω₁. */
  readonly mode?: Float64Array;
}

interface Built {
  readonly beam?: Beam;
  readonly plate?: Plate;
  readonly Q: number;
  readonly c: Float64Array;
  readonly ci?: Float64Array;
  readonly solveStatic: () => Float64Array;
}

interface LevelData {
  /** The beam's tabulation, or the plate's pair. */
  readonly tab: Tabulation | PlateTables;
  readonly field: PointField;
  readonly e0: Float64Array;
  readonly mu0: Float64Array;
}

export class Sampler {
  readonly M: number;
  private readonly levels = new Map<number, LevelData>();
  private readonly harmonic?: Harmonic;

  /**
   * `kl` is the expansion along the beam — or along the plate's x, for length
   * ℓ/a — and `klY` along the plate's y (the same object for a square plate).
   * A plate keeps the M largest products of the two (`random/field2d.ts`).
   */
  constructor(readonly spec: McSpec, readonly kl: KLData, readonly klY: KLData = kl) {
    this.M = spec.plate ? Math.min(spec.field.terms, termsOf(kl) * termsOf(klY)) : termsOf(kl);
    if (spec.qoi === "response") this.harmonic = spec.plate ? plateHarmonicOf(spec.plate) : harmonicOf(spec.beam);
  }

  /** Tables for a level: built on first use, then shared by every sample on it. */
  private level(ne: number): LevelData {
    let L = this.levels.get(ne);
    if (!L) {
      const { p, k, beam, plate } = this.spec;
      if (plate) {
        const tab = plateTables(p, ne, k, plate.aspect), n = tab.tx.x.length * tab.ty.x.length;
        const ones = new Float64Array(n).fill(1);
        L = { tab, field: new FieldOnGrid(this.kl, this.klY, this.M, tab.tx.x, tab.ty.x, plate.aspect), e0: ones, mu0: ones };
      } else {
        const tab = tabulate(splineSpace(p, ne, k), gaussLegendre(p + 1), 2);
        const sec = sectionOf(beam.section);
        L = {
          tab,
          field: new FieldAt(this.kl, tab.x),
          e0: sampleAt(tab, sec.stiffness ?? (() => 1)),
          mu0: sampleAt(tab, sec.mass ?? (() => 1)),
        };
      }
      this.levels.set(ne, L);
    }
    return L;
  }

  /** Draw sample `index` of `stream`, build its structure on `ne` elements, and evaluate Q. */
  private build(index: number, ne: number, stream: number): Built {
    const { spec } = this, L = this.level(ne);
    const r = realise(L.field, spec.field, draw(spec.field, this.M, spec.seed, index, stream), L.e0, L.mu0);
    if (spec.plate) {
      const pc = spec.plate;
      const plate = new Plate({ p: spec.p, k: spec.k, ne, edges: pc.edges, aspect: pc.aspect, nu: pc.nu, stiffness: r.e, mass: r.mu }, L.tab as PlateTables);
      const load = plateLoadOf(pc, r.q, r.P);
      const { Q, c, ci } = evaluatePlateQoI(plate, spec.qoi, load, plateQoiPoint(pc), this.harmonic);
      return { plate, Q, c, ci, solveStatic: () => plate.solve(load) };
    }
    const beam = new Beam({ p: spec.p, k: spec.k, ne, supports: SUPPORTS[spec.beam.supports], stiffness: r.e, mass: r.mu }, L.tab as Tabulation);
    const xq = qoiPoint(spec.beam.supports);
    const load = spec.beam.load === "uniform" ? { q: r.q ?? (() => 1) } : { forces: [{ x: xq, P: r.P }] };
    const { Q, c, ci } = evaluateQoI(beam, spec.qoi, load, xq, this.harmonic);
    return { beam, Q, c, ci, solveStatic: () => beam.solve(load) };
  }

  /** Q for sample `index` of `stream` on the mesh of `ne` elements; the deflection too if `wantW`. */
  solve(index: number, ne: number, wantW = false, stream = 0): { Q: number; w: Float64Array | null } {
    const b = this.build(index, ne, stream);
    if (!wantW) return { Q: b.Q, w: null };
    const cs = b.c.length ? b.c : b.solveStatic();
    let w: Float64Array, wi: Float64Array | null = null;
    if (b.plate) {
      const xs = fieldX.map((x) => x * b.plate!.aspect), mid = [0.5];
      w = b.plate.grid(cs, xs, mid);
      if (b.ci) wi = b.plate.grid(b.ci, xs, mid);
    } else {
      w = fieldX.map((x) => b.beam!.evaluate(cs, x)[0]);
      if (b.ci) wi = fieldX.map((x) => b.beam!.evaluate(b.ci!, x)[0]);
    }
    return { Q: b.Q, w: wi ? w.map((v, i) => Math.hypot(v, wi![i])) : w };
  }

  /** Sample `index` of `stream` on `ne` elements, kept whole for drawing. */
  inspect(index: number, ne: number, stream = 0): Inspected {
    const b = this.build(index, ne, stream);
    const w = this.spec.qoi === "omega1" || this.spec.qoi === "response" ? b.solveStatic() : b.c;
    let mode: Float64Array | undefined;
    if (this.spec.qoi === "omega1") {
      mode = (b.beam ?? b.plate)!.modes(1).vectors[0];
      let peak = 0;
      for (const v of mode) if (Math.abs(v) > Math.abs(peak)) peak = v;
      if (peak < 0) mode = mode.map((v) => -v);
    }
    const u = b.ci ? { re: b.c, im: b.ci } : undefined;
    return { Q: b.Q, beam: b.beam, plate: b.plate, w, u, mode };
  }

  /** Sample `index` of `stream` on `ne` elements and, if `coarse`, on ne/2 with the same ω. */
  sample(index: number, ne: number, coarse: boolean, wantW = false, stream = 0): Sample {
    const fine = this.solve(index, ne, wantW, stream);
    const Qc = coarse ? this.solve(index, ne / 2, false, stream).Q : NaN;
    return { Q: fine.Q, Qc, w: fine.w };
  }
}
