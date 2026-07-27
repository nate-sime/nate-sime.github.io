/**
 * The two Nusselt numbers against time — the data behind the plot in the
 * bottom-left corner, and the scales it is drawn on.
 *
 * Nu is the diagnostic this solver is *verified* by, and the readout's two
 * numbers show only half of what it says. The inner and outer fluxes are
 * independent reductions over independent boundary rows, so their agreement at
 * steady state is a genuine global heat balance rather than a fitted quantity —
 * but "at steady state" is the load-bearing clause, and a pair of instantaneous
 * numbers cannot tell you whether the run has settled, is still oscillating, or
 * is mid-way through absorbing a change to Ra. That is a time series, so it is
 * shown as one.
 *
 * **Both series share one axis, and must.** They are the same quantity at two
 * radii; the whole reading is whether the two curves have met. A second y scale
 * would let them meet at an arbitrary offset, which is the one thing this plot
 * exists to rule out.
 *
 * DOM-free on purpose, for the same reason `equation.ts` is: the ring buffer's
 * restart rule and the tick arithmetic are the parts that can go silently wrong,
 * and neither needs a browser to check. `nuplot.ts` renders it.
 */

/** A boundary, in the order the readout names them. */
export type NuSeries = "inner" | "outer";

/**
 * Series colours: categorical slots 1 and 2 of the site's reference palette,
 * stepped for a dark surface, and validated as a pair against every surface the
 * panel can sit on — worst-case CVD ΔE 26.8, normal-vision 31.8, both above
 * 3:1 contrast even where the brightest part of the field shows through the
 * panel's 8% transparency.
 *
 * The *assignment* is semantic rather than slot order. This page already spends
 * a colour scale on temperature: warm is hot, dark and cool is cold, everywhere
 * on screen. The inner boundary is the hot one and the outer is the cold one, so
 * the warm slot goes to `inner`. Handing `inner` slot 1 because it is named
 * first in the readout would put the cold boundary in orange in a window where
 * orange means hot — a legend can label that, but it cannot stop it misreading
 * at a glance.
 */
export const NU_COLOUR: Record<NuSeries, string> = {
  inner: "#d95926",
  outer: "#3987e5",
};

export interface NuSample {
  /** Simulation time the reduction was taken at. */
  t: number;
  /** Step count at that time — what the display window is measured in. */
  step: number;
  inner: number;
  outer: number;
}

export interface NuExtent {
  t0: number;
  t1: number;
  /** Over *both* series: they share the axis. */
  lo: number;
  hi: number;
}

/**
 * Samples retained. Diagnostics arrive at the poll rate — once per 15 frames, so
 * roughly 4–8 a second — making this half an hour or so of history, and the
 * ceiling on what the "all" display window can show. Deliberately generous: it
 * is the *display* window that decides how much is on screen, and the buffer's
 * only job is to not be the thing that limits it. Three typed arrays plus the
 * step axis is about 400 kB, against tens of megabytes of quadrature buffers in
 * the Krylov tier.
 *
 * Retention is a sample count, not an interval of simulation time: `speed` spans
 * a factor of 256, so a fixed time budget would hold minutes at one end of that
 * list and milliseconds at the other.
 */
export const NU_CAPACITY = 16_384;

/** A rolling buffer of the last `capacity` polls. */
export class NuTrace {
  private readonly ts: Float64Array;
  /** f64, not i32: at 16 steps a frame a run passes 2³¹ steps inside a fortnight. */
  private readonly st: Float64Array;
  private readonly vi: Float32Array;
  private readonly vo: Float32Array;
  /** Next slot to write; the oldest sample held is `head − count`. */
  private head = 0;
  private count = 0;

  constructor(readonly capacity = NU_CAPACITY) {
    if (!Number.isInteger(capacity) || capacity < 2)
      throw new Error(`capacity ${capacity} must be an integer ≥ 2`);
    this.ts = new Float64Array(capacity);
    this.st = new Float64Array(capacity);
    this.vi = new Float32Array(capacity);
    this.vo = new Float32Array(capacity);
  }

  get length(): number { return this.count; }

  clear(): void { this.head = 0; this.count = 0; }

  /**
   * Record one poll, and say whether anything was stored — the plot redraws on a
   * new sample rather than on every frame.
   *
   * Two things are dropped rather than plotted. A non-finite value: the
   * reductions read NaN until the first readback lands, and the readout already
   * says so with a dash. And a repeat of the last `t`: the poll is asynchronous
   * and deliberately allowed to be one behind, so the frame loop sees the same
   * stats several times over.
   *
   * A clock *before* the last one restarts the buffer instead of extending it.
   * Reseeding sets both the time and the step count back to zero, and so does a
   * rebuild; either way the samples on the far side of that are from a different
   * run and joining them with a line would draw a trajectory nothing followed.
   * Making that the rule here rather than a `clear()` at each call site means a
   * poll still in flight across the reset cannot slip a stale sample in behind it.
   */
  push({ t, step, inner, outer }: NuSample): boolean {
    if (!(Number.isFinite(t) && Number.isFinite(step)
          && Number.isFinite(inner) && Number.isFinite(outer)))
      return false;
    if (this.count > 0) {
      const k = (this.head + this.capacity - 1) % this.capacity;
      if (t < this.ts[k] || step < this.st[k]) this.clear();
      else if (t === this.ts[k]) return false;
    }
    this.ts[this.head] = t;
    this.st[this.head] = step;
    this.vi[this.head] = inner;
    this.vo[this.head] = outer;
    this.head = (this.head + 1) % this.capacity;
    if (this.count < this.capacity) this.count++;
    return true;
  }

  /**
   * Index of the oldest sample lying within `steps` of the newest — the start of
   * the slice the plot draws. `Infinity` gives 0, the whole buffer.
   *
   * Scanned backwards from the newest rather than bisected: it stops after the
   * samples it is going to return, so a narrow window costs a walk over that
   * window and not over the buffer. The step axis is monotone by construction
   * (`push` restarts on a clock that goes backwards), so the first sample past
   * the cutoff ends it.
   */
  first(steps: number): number {
    if (this.count === 0 || !(steps > 0)) return 0;
    if (!Number.isFinite(steps)) return 0;
    const cutoff = this.st[this.index(this.count - 1)] - steps;
    let i = this.count - 1;
    while (i > 0 && this.st[this.index(i - 1)] >= cutoff) i--;
    return i;
  }

  /** Sample `i`, counting from the oldest held. */
  at(i: number): NuSample {
    const k = this.index(i);
    return { t: this.ts[k], step: this.st[k], inner: this.vi[k], outer: this.vo[k] };
  }

  get last(): NuSample | null {
    return this.count === 0 ? null : this.at(this.count - 1);
  }

  /**
   * Bounds of the samples from `from` onwards, `null` if there are none. Taken
   * over the *displayed* slice and not the whole buffer, so narrowing the window
   * rescales the y axis onto what is on screen — which is the point of narrowing
   * it: a settled band a thousandth of the transient's height is invisible on an
   * axis the transient still sets.
   */
  extent(from = 0): NuExtent | null {
    if (from >= this.count || from < 0) return null;
    let lo = Infinity, hi = -Infinity;
    for (let i = from; i < this.count; i++) {
      const k = this.index(i);
      lo = Math.min(lo, this.vi[k], this.vo[k]);
      hi = Math.max(hi, this.vi[k], this.vo[k]);
    }
    return {
      t0: this.ts[this.index(from)],
      t1: this.ts[this.index(this.count - 1)],
      lo, hi,
    };
  }

  private index(i: number): number {
    if (!Number.isInteger(i) || i < 0 || i >= this.count)
      throw new RangeError(`sample ${i} of ${this.count}`);
    return (this.head - this.count + i + this.capacity) % this.capacity;
  }
}

/**
 * The 1–2–5 step nearest `span / target`, nearest on a log scale — so `target`
 * is what it says, roughly, and the count lands between about 2 and 4 ticks
 * either side of it. On a panel this short that is the whole useful range: a
 * fourth gridline across 54 px is texture, not a scale.
 */
export function niceStep(span: number, target: number): number {
  if (!(span > 0) || !(target > 0)) return 1;
  const raw = span / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const n = raw / mag;
  return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
}

export interface NuAxis {
  lo: number;
  hi: number;
  /** Rounded tick values inside [lo, hi]; possibly empty. */
  ticks: number[];
  /** Decimals to print every tick at, so the column reads as a column. */
  decimals: number;
}

/**
 * An axis over a data extent: padded so the curves never run into the frame,
 * with rounded ticks at whatever multiples of the step fall *inside* that range.
 *
 * The extent is padded rather than snapped out to the bracketing ticks. Snapping
 * gives a tick on each edge, which looks tidier, but it can inflate the span by
 * two-thirds — and on this plot the vertical resolution *is* the reading. Once
 * the two curves are within a per cent of each other, the question is whether
 * they are still approaching, and rounding the axis outwards is exactly what
 * would hide it.
 *
 * A flat extent — one sample, or a settled trace whose f32 values stopped
 * moving — gets a band around the value instead of a zero-height axis.
 */
export function niceAxis(lo: number, hi: number, target = 3, pad = 0.08): NuAxis {
  if (!(Number.isFinite(lo) && Number.isFinite(hi))) return { lo: 0, hi: 1, ticks: [], decimals: 0 };
  if (hi < lo) [lo, hi] = [hi, lo];
  const half = hi > lo ? (hi - lo) * pad : Math.max(Math.abs(hi), 1) * 0.05;
  const a = lo - half, b = hi + half;
  const step = niceStep(b - a, target);
  const ticks: number[] = [];
  // Bounded independently of `step`: this runs every time a sample lands, and a
  // step that came out subnormal must not spin the frame loop.
  for (let k = Math.ceil(a / step); ticks.length < 64 && k * step <= b; k++)
    ticks.push(k * step);
  return { lo: a, hi: b, ticks, decimals: decimalsFor(step) };
}

/** Decimals a value stepping by `step` needs to be printed at without repeats. */
export const decimalsFor = (step: number): number =>
  step > 0 ? Math.max(0, Math.min(6, -Math.floor(Math.log10(step)))) : 0;

/**
 * Decimals for the two ends of the time axis — a *shared* count, since they are
 * two coordinates on one axis and a pair printed at different precisions reads as
 * a pair of unrelated numbers.
 *
 * A tenth of the span resolves the window, which is what the labels are for. The
 * exception is a window whose start is small but not zero: the first poll of a run
 * lands a few steps in, so t₀ is something like 0.003 against a span of order one,
 * and a tenth of the span rounds that to "0.0". On its own that would be a fair
 * label for a value negligible on its own axis — but the dimensional row beneath
 * it resolves the same instant to three significant figures whatever the span, and
 * prints "3.86 Gyr". Both are correct, and together they read as the two rows
 * contradicting each other. So the count rises to whatever shows t₀'s leading
 * digit.
 */
export function axisDecimals(t0: number, span: number): number {
  const base = decimalsFor(span / 10);
  if (!Number.isFinite(t0) || t0 === 0) return base;
  return Number(t0.toFixed(base)) === 0 ? decimalsFor(Math.abs(t0)) : base;
}
