// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Running statistics for a Monte Carlo run, merged in sample order.
 *
 * Batches come back from the workers in whatever order they finish. They are
 * held until every earlier sample has arrived and only then folded in, one
 * sample at a time, by Welford's update — so the statistics after N samples
 * are those of samples 0 … N − 1, bit for bit, however many workers produced
 * them and however the batches were cut. That is what makes a run with a given
 * seed a reproducible experiment rather than a race.
 *
 * What is kept: every Q and Q_c (for the histogram, and for the moments of
 * any prefix of the run), the moments of Q, of its parent-level partner Q_c and
 * of the correction Y = Q − Q_c, the pointwise moments of the deflection field,
 * the first few sample paths, and a trajectory of the running mean at
 * geometrically spaced N (and every power of ten) for the convergence plots.
 * `view(n)` gives all of it over a prefix — the field's moments only at the
 * powers of ten, where they are kept.
 */

export class Moments {
  n = 0;
  mean = 0;
  private m2 = 0;

  /** A Welford state taken from `state()`. */
  static of([n, mean, m2]: readonly number[]): Moments {
    const m = new Moments();
    m.n = n; m.mean = mean; m.m2 = m2;
    return m;
  }

  state(): [number, number, number] {
    return [this.n, this.mean, this.m2];
  }

  add(x: number): void {
    this.n++;
    const d = x - this.mean;
    this.mean += d / this.n;
    this.m2 += d * (x - this.mean);
  }

  /** Unbiased sample variance; NaN below two samples. */
  get variance(): number {
    return this.n > 1 ? this.m2 / (this.n - 1) : NaN;
  }

  get sd(): number {
    return Math.sqrt(this.variance);
  }

  /** Standard error of the mean, σ̂/√N. */
  get se(): number {
    return this.sd / Math.sqrt(this.n);
  }
}

/** Pointwise Welford over a vector-valued sample. */
export class FieldMoments {
  n = 0;
  readonly mean: Float64Array;
  private readonly m2: Float64Array;

  constructor(readonly size: number) {
    this.mean = new Float64Array(size);
    this.m2 = new Float64Array(size);
  }

  add(w: ArrayLike<number>, offset = 0): void {
    this.n++;
    for (let i = 0; i < this.size; i++) {
      const x = w[offset + i], d = x - this.mean[i];
      this.mean[i] += d / this.n;
      this.m2[i] += d * (x - this.mean[i]);
    }
  }

  sd(): Float64Array {
    return this.m2.map((v) => (this.n > 1 ? Math.sqrt(v / (this.n - 1)) : NaN));
  }

  clone(): FieldMoments {
    const f = new FieldMoments(this.size);
    f.n = this.n;
    f.mean.set(this.mean);
    f.m2.set(this.m2);
    return f;
  }
}

/** What a worker returns for samples [from, from + count). */
export interface Batch {
  readonly from: number;
  readonly count: number;
  readonly Q: Float64Array;
  /** Parent-level Q, same ω; NaN on the coarsest level. */
  readonly Qc: Float64Array;
  /** Deflection at the field points, flat [sample·size + i]; null if not recorded. */
  readonly W: Float64Array | null;
  /** Worker time spent on the batch, ms. */
  readonly ms: number;
}

export interface Checkpoint {
  readonly n: number;
  readonly mean: number;
  readonly sd: number;
  /** Mean and sd of the correction Y = Q − Q_c (NaN on the coarsest level). */
  readonly dMean: number;
  readonly dSd: number;
}

/** What plain Monte Carlo draws: the statistics of samples 0 … n − 1. */
export interface Stats {
  readonly n: number;
  readonly q: Moments;
  /** Moments of Q − Q_c; empty on a stream with no parent level. */
  readonly dq: Moments;
  readonly min: number;
  readonly max: number;
  readonly values: Float64Array;
  readonly trajectory: readonly Checkpoint[];
  /** Pointwise moments of the deflection, over its own field.n ≤ n samples. */
  readonly field: FieldMoments;
}

/** Sample paths kept for drawing. */
export const PATHS = 12;

/** Every this many samples the running moments are kept, so a prefix's moments cost at most this many updates. */
const SNAP = 1024;

export class Accumulator {
  /** Samples folded in: always a prefix 0 … n − 1. */
  n = 0;
  readonly q = new Moments();
  readonly qc = new Moments();
  readonly dq = new Moments();
  /**
   * The multilevel correction: Y = Q − Q_c, or Q itself on a stream with no
   * parent level (Y₀ = Q₀, the first term of the telescoping sum).
   */
  readonly y = new Moments();
  readonly field: FieldMoments;
  readonly paths: Float64Array[] = [];
  readonly trajectory: Checkpoint[] = [];
  min = Infinity;
  max = -Infinity;
  /** Worker time summed over every batch folded in, ms. */
  cpuMs = 0;
  private Qs = new Float64Array(1024);
  private Qcs = new Float64Array(1024);
  /** Welford states of Q and Y at n = 0, SNAP, 2·SNAP, … — six numbers each. */
  private snaps: number[] = [0, 0, 0, 0, 0, 0];
  private readonly pending = new Map<number, Batch>();
  private nextCheckpoint = 2;
  /** The field's moments at n = 1, 10, 100, …: a prefix's field, which Welford cannot undo. */
  private readonly fieldSnaps: FieldMoments[] = [];
  private nextDecade = 1;

  constructor(readonly fieldSize: number) {
    this.field = new FieldMoments(fieldSize);
  }

  /** Every Q folded in so far, in sample order. */
  get values(): Float64Array {
    return this.Qs.subarray(0, this.n);
  }

  /** Every Q_c folded in so far (NaN on a stream with no parent level). */
  get coarseValues(): Float64Array {
    return this.Qcs.subarray(0, this.n);
  }

  /**
   * The moments of Q and of Y over samples 0 … n − 1 alone, n ≤ this.n — what
   * the run had said when it had n samples. Multilevel Monte Carlo asks this of
   * one store for many sample counts: every tolerance reads its own prefix.
   */
  prefix(n: number): { q: Moments; y: Moments } {
    if (n > this.n) throw new Error(`prefix of ${n} samples asked of ${this.n}`);
    const k = Math.floor(n / SNAP), o = 6 * k;
    const q = Moments.of(this.snaps.slice(o, o + 3)), y = Moments.of(this.snaps.slice(o + 3, o + 6));
    for (let i = k * SNAP; i < n; i++) {
      const Q = this.Qs[i], Qc = this.Qcs[i];
      q.add(Q);
      y.add(Number.isNaN(Qc) ? Q : Q - Qc);
    }
    return { q, y };
  }

  /**
   * The statistics of samples 0 … n − 1 alone: what the run had shown at n,
   * so lowering N shows fewer samples without throwing the rest away. The
   * field is exact at a power of ten (every N the view offers) and otherwise
   * the last power of ten below n.
   */
  view(n: number): Stats {
    if (n >= this.n) return this;
    const { q, y } = this.prefix(n);
    const values = this.values.subarray(0, n);
    let min = Infinity, max = -Infinity;
    for (const v of values) { if (v < min) min = v; if (v > max) max = v; }
    // A stream has a parent level for every sample or for none, so Y is Q − Q_c throughout or Q throughout.
    const dq = this.dq.n ? y : new Moments();
    const field = this.fieldSnaps.filter((f) => f.n <= n).pop() ?? new FieldMoments(this.fieldSize);
    return { n, q, dq, min, max, values, trajectory: this.trajectory.filter((c) => c.n <= n), field };
  }

  push(b: Batch): void {
    if (b.from < this.n) return; // already folded (a duplicate)
    this.pending.set(b.from, b);
    for (let next = this.pending.get(this.n); next; next = this.pending.get(this.n)) {
      this.pending.delete(this.n);
      this.fold(next);
    }
  }

  private fold(b: Batch): void {
    if (this.n + b.count > this.Qs.length) {
      const size = Math.max(2 * this.Qs.length, this.n + b.count);
      const grow = (a: Float64Array) => { const g = new Float64Array(size); g.set(a.subarray(0, this.n)); return g; };
      this.Qs = grow(this.Qs);
      this.Qcs = grow(this.Qcs);
    }
    const size = this.fieldSize;
    for (let s = 0; s < b.count; s++) {
      const Q = b.Q[s], Qc = b.Qc[s];
      this.Qs[this.n] = Q;
      this.Qcs[this.n] = Qc;
      this.q.add(Q);
      if (Q < this.min) this.min = Q;
      if (Q > this.max) this.max = Q;
      if (!Number.isNaN(Qc)) { this.qc.add(Qc); this.dq.add(Q - Qc); }
      this.y.add(Number.isNaN(Qc) ? Q : Q - Qc);
      if (b.W) {
        this.field.add(b.W, s * size);
        if (this.paths.length < PATHS) this.paths.push(b.W.slice(s * size, (s + 1) * size));
      }
      this.n++;
      if (this.n % SNAP === 0) this.snaps.push(...this.q.state(), ...this.y.state());
      if (this.n === this.nextDecade) {
        if (size) this.fieldSnaps.push(this.field.clone());
        this.nextDecade *= 10;
      }
      if (this.n >= this.nextCheckpoint) {
        this.trajectory.push({
          n: this.n, mean: this.q.mean, sd: this.q.sd,
          dMean: this.dq.n ? this.dq.mean : NaN, dSd: this.dq.n > 1 ? this.dq.sd : NaN,
        });
        // Every power of ten is a checkpoint too, so a view of the first 10ᵏ ends on one.
        this.nextCheckpoint = Math.min(Math.max(this.n + 1, Math.ceil(this.n * 1.06)), this.nextDecade);
      }
    }
    this.cpuMs += b.ms;
  }

  /** Samples received but waiting on an earlier batch. */
  get waiting(): number {
    let s = 0;
    for (const b of this.pending.values()) s += b.count;
    return s;
  }
}

/** Two-sided 95% normal quantile. */
export const Z95 = 1.959963984540054;

export interface Histogram {
  readonly edges: Float64Array;
  /** Density: counts / (N · width), so it overlays a pdf. */
  readonly density: Float64Array;
}

/**
 * A density histogram on [lo, hi] with Scott's bin width 3.49 σ̂ N^{−1/3}, held
 * between 8 and 80 bins — enough to show a skew, not so many that it is noise.
 */
export function histogram(values: ArrayLike<number>, lo: number, hi: number, sd: number): Histogram {
  const n = values.length;
  const width = Number.isFinite(sd) && sd > 0 ? 3.49 * sd * n ** (-1 / 3) : (hi - lo) / 8;
  const bins = Math.max(8, Math.min(80, Math.ceil((hi - lo) / width) || 8));
  const span = hi > lo ? hi - lo : Math.abs(lo) * 1e-6 || 1e-12;
  const edges = Float64Array.from({ length: bins + 1 }, (_, i) => lo + (span * i) / bins);
  const counts = new Float64Array(bins);
  for (let i = 0; i < n; i++) counts[Math.min(bins - 1, Math.floor(((values[i] - lo) / span) * bins))]++;
  const h = span / bins;
  return { edges, density: counts.map((c) => c / (n * h)) };
}
