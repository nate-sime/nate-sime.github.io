/**
 * The Nusselt trace's data layer. No DOM: what is checked here is the ring
 * buffer's restart rule and the tick arithmetic, and neither is a property of
 * the drawing.
 *
 * Two failures are worth guarding against because both are silent. A stale
 * sample surviving a reseed joins two runs with a straight line, producing a
 * trajectory nothing followed — and it looks exactly like a transient. And a
 * degenerate axis (one sample, or a settled trace whose f32 values stopped
 * moving) divides by a zero span: NaN coordinates draw nothing at all, so the
 * panel would simply go blank at the moment the run reached the state the plot
 * exists to show.
 */

import { describe, it, expect } from "vitest";
import {
  NU_COLOUR, NuTrace, axisDecimals, decimalsFor, niceAxis, niceStep, type NuSeries,
} from "../src/ui/nusselt";
import {
  REFERENCE, TIME_UNIT, TIME_UNIT_YEARS, YEAR, dimensionalTime, referenceNote,
} from "../src/ui/dimensional";

const series = Object.keys(NU_COLOUR) as NuSeries[];

/** A sample whose step count tracks `t`, for the cases that only care about t. */
const sample = (t: number, inner: number, outer: number) =>
  ({ t, step: Math.round(t * 1e4), inner, outer });

describe("series colours", () => {
  it("names both boundaries, with distinct colours", () => {
    expect(series).toEqual(["inner", "outer"]);
    expect(NU_COLOUR.inner).not.toBe(NU_COLOUR.outer);
  });

  it.each(series)("%s is an opaque hex the canvas and the legend can share", (s) => {
    expect(NU_COLOUR[s]).toMatch(/^#[0-9a-f]{6}$/);
  });

  /**
   * The one property of the *assignment* rather than the palette: this page
   * spends a colour scale on temperature, and the inner boundary is the hot one.
   * A cool inner curve would read as the cold boundary in a window where warm
   * means hot everywhere else. Compared in red-minus-blue, which is all the
   * claim needs.
   */
  it("gives the hot boundary the warm colour", () => {
    const warmth = (hex: string) =>
      parseInt(hex.slice(1, 3), 16) - parseInt(hex.slice(5, 7), 16);
    expect(warmth(NU_COLOUR.inner)).toBeGreaterThan(warmth(NU_COLOUR.outer));
  });
});

describe("NuTrace", () => {
  it("starts empty and reports nothing to scale to", () => {
    const tr = new NuTrace(8);
    expect(tr.length).toBe(0);
    expect(tr.last).toBeNull();
    expect(tr.extent()).toBeNull();
  });

  it("holds samples oldest-first", () => {
    const tr = new NuTrace(8);
    expect(tr.push({ t: 0.1, step: 10, inner: 2, outer: 3 })).toBe(true);
    expect(tr.push({ t: 0.2, step: 20, inner: 4, outer: 5 })).toBe(true);
    expect(tr.length).toBe(2);
    expect(tr.at(0)).toEqual({ t: 0.1, step: 10, inner: 2, outer: 3 });
    expect(tr.at(1)).toEqual({ t: 0.2, step: 20, inner: 4, outer: 5 });
    expect(tr.last).toEqual({ t: 0.2, step: 20, inner: 4, outer: 5 });
  });

  it("rolls the oldest out once full, and stays ordered across the wrap", () => {
    const tr = new NuTrace(4);
    for (let i = 0; i < 10; i++) tr.push(sample(i, i, -i));
    expect(tr.length).toBe(4);
    expect(Array.from({ length: 4 }, (_, i) => tr.at(i).t)).toEqual([6, 7, 8, 9]);
    expect(tr.extent()).toMatchObject({ t0: 6, t1: 9, lo: -9, hi: 9 });
  });

  // The readout prints a dash until the first readback lands; the plot must not
  // plot the NaN behind it.
  it("drops non-finite samples", () => {
    const tr = new NuTrace(8);
    expect(tr.push({ t: NaN, step: NaN, inner: NaN, outer: NaN })).toBe(false);
    expect(tr.push({ t: 0.1, step: 1, inner: NaN, outer: 3 })).toBe(false);
    expect(tr.push({ t: 0.1, step: 1, inner: 2, outer: Infinity })).toBe(false);
    // The step axis is what the display window is measured against, so a sample
    // without one is no more plottable than a sample without a value.
    expect(tr.push({ t: 0.1, step: NaN, inner: 2, outer: 3 })).toBe(false);
    expect(tr.length).toBe(0);
  });

  // The poll is asynchronous and allowed to be one behind, so the frame loop
  // offers the same reading many times over.
  it("drops a repeat of the sample it already holds", () => {
    const tr = new NuTrace(8);
    expect(tr.push(sample(0.1, 2, 3))).toBe(true);
    expect(tr.push(sample(0.1, 2, 3))).toBe(false);
    expect(tr.push(sample(0.1, 9, 9))).toBe(false);
    expect(tr.length).toBe(1);
  });

  /**
   * The rule the whole class exists for. Reseed sets the clock to zero, and a
   * poll already in flight can land after it — so the guard has to be here,
   * where every sample passes, and not at the call sites that know a reset
   * happened. Either counter going backwards is enough.
   */
  it("restarts the buffer when the clock goes backwards", () => {
    const tr = new NuTrace(8);
    for (const t of [0.1, 0.2, 0.3]) tr.push(sample(t, 5, 5));
    expect(tr.push(sample(0, 1, 1))).toBe(true);
    expect(tr.length).toBe(1);
    expect(tr.at(0)).toMatchObject({ t: 0, inner: 1, outer: 1 });
  });

  it("restarts the buffer when the step count goes backwards", () => {
    const tr = new NuTrace(8);
    for (const s of [100, 200, 300]) tr.push({ t: s / 1000, step: s, inner: 5, outer: 5 });
    // A rebuild at a smaller dt: the step count restarts while t happens to land
    // beyond where the previous run had reached.
    expect(tr.push({ t: 0.4, step: 5, inner: 1, outer: 1 })).toBe(true);
    expect(tr.length).toBe(1);
  });

  it("refuses to index outside the buffer", () => {
    const tr = new NuTrace(4);
    tr.push(sample(1, 1, 1));
    expect(() => tr.at(1)).toThrow(RangeError);
    expect(() => tr.at(-1)).toThrow(RangeError);
  });

  it("rejects a capacity that cannot hold a line", () => {
    expect(() => new NuTrace(1)).toThrow();
    expect(() => new NuTrace(2.5)).toThrow();
  });
});

/**
 * The display window. Measured in solver steps, which is the point: the poll rate
 * is the frame loop's, so a window counted in samples would cover a different
 * stretch of the simulation at every playback speed.
 */
describe("NuTrace.first", () => {
  /** 20 samples, 50 steps apart: steps 0, 50, … 950. */
  const filled = (): NuTrace => {
    const tr = new NuTrace(64);
    for (let i = 0; i < 20; i++)
      tr.push({ t: i * 0.05, step: i * 50, inner: 2, outer: 2 });
    return tr;
  };

  it("shows everything for an unbounded window", () => {
    expect(filled().first(Infinity)).toBe(0);
    expect(new NuTrace(8).first(Infinity)).toBe(0);
  });

  it("keeps exactly the samples within the window of the newest", () => {
    const tr = filled();
    // Newest is step 950; a 200-step window reaches back to 750, which is
    // samples at 750, 800, 850, 900, 950 — five of them, indices 15…19.
    expect(tr.first(200)).toBe(15);
    expect(tr.at(tr.first(200)).step).toBe(750);
    expect(tr.first(100)).toBe(17);
    // Wider than the run: everything, not an out-of-range index.
    expect(tr.first(1e6)).toBe(0);
  });

  it("always leaves at least the newest sample to draw", () => {
    const tr = filled();
    for (const w of [0, 1, 49, -5, NaN]) {
      const from = tr.first(w);
      expect(from).toBeGreaterThanOrEqual(0);
      expect(from).toBeLessThan(tr.length);
    }
    expect(tr.first(1)).toBe(19);
  });

  it("is empty-safe", () => {
    const tr = new NuTrace(8);
    expect(tr.first(500)).toBe(0);
    expect(tr.extent(tr.first(500))).toBeNull();
  });

  /**
   * The reason the window exists. A settled band is a thousandth of the initial
   * transient's height, so on an axis the transient still sets it is a flat line;
   * scaling to the visible slice is what resolves it.
   */
  it("scales the extent to the window, not to the whole run", () => {
    const tr = new NuTrace(64);
    for (let i = 0; i < 20; i++)
      tr.push({
        t: i * 0.05, step: i * 50,
        inner: i < 4 ? 1 : 2.48, outer: i < 4 ? 1 : 2.4801,
      });

    // Over the whole run the transient sets the axis, and the settled band it
    // ends in is 3e-4 of a span of 1.5 — a flat line at the top of the plot.
    const all = tr.extent()!;
    expect(all.lo).toBe(1);
    expect(niceAxis(all.lo, all.hi).hi - niceAxis(all.lo, all.hi).lo)
      .toBeGreaterThan(1.5);

    // Over the last 200 steps the transient is gone and the axis is the band.
    const win = tr.extent(tr.first(200))!;
    expect(win.lo).toBeCloseTo(2.48, 6);
    expect(niceAxis(win.lo, win.hi).hi - niceAxis(win.lo, win.hi).lo)
      .toBeLessThan(0.01);
  });
});

describe("niceStep", () => {
  it.each([
    [1, 3, 0.5],
    [6, 3, 2],
    [10, 3, 5],
    [0.03, 3, 0.01],
    [300, 3, 100],
    [4, 4, 1],
  ])("span %f into ~%i steps → %f", (span, target, want) => {
    expect(niceStep(span, target)).toBeCloseTo(want, 12);
  });

  // `target` is honoured to within the coarseness of a 1–2–5 ladder, which is
  // what the axis is allowed to be: 2 to 5 gridlines, never one and never ten.
  it("keeps the tick count near the target across every decade", () => {
    for (let e = -6; e <= 6; e++)
      for (const m of [1, 1.3, 1.7, 2.4, 3.9, 6.2, 8.8]) {
        const span = m * 10 ** e;
        const count = span / niceStep(span, 3);
        expect(count).toBeGreaterThanOrEqual(2);
        expect(count).toBeLessThanOrEqual(5);
      }
  });

  it("only ever returns a 1, 2 or 5 decade", () => {
    for (let e = -6; e <= 6; e++)
      for (const m of [1, 1.3, 1.7, 2.4, 3.9, 6.2, 8.8]) {
        const step = niceStep(m * 10 ** e, 3);
        const mant = step / 10 ** Math.floor(Math.log10(step) + 1e-9);
        expect([1, 2, 5]).toContain(Math.round(mant * 10) / 10);
      }
  });

  it("does not divide by a degenerate span", () => {
    expect(niceStep(0, 3)).toBe(1);
    expect(niceStep(-1, 3)).toBe(1);
    expect(niceStep(NaN, 3)).toBe(1);
  });
});

describe("niceAxis", () => {
  it("pads the extent so a curve never runs into the frame", () => {
    const ax = niceAxis(2, 6);
    expect(ax.lo).toBeLessThan(2);
    expect(ax.hi).toBeGreaterThan(6);
  });

  /**
   * Not snapped out to the bracketing ticks. Once the two curves are within a
   * per cent of each other the question is whether they are still approaching,
   * and rounding the axis outwards is what would hide it — so the padded span
   * must stay close to the data's.
   */
  it("keeps the span close to the data's rather than rounding it out", () => {
    const ax = niceAxis(4.80, 4.90);
    expect(ax.hi - ax.lo).toBeLessThan(0.13);
  });

  it("puts every tick inside the axis, at a multiple of the step", () => {
    for (const [lo, hi] of [[1, 2], [0.004, 0.009], [12, 480], [-3, 7]]) {
      const ax = niceAxis(lo, hi);
      expect(ax.ticks.length).toBeGreaterThan(0);
      for (const t of ax.ticks) {
        expect(t).toBeGreaterThanOrEqual(ax.lo);
        expect(t).toBeLessThanOrEqual(ax.hi);
      }
      const step = ax.ticks.length > 1 ? ax.ticks[1] - ax.ticks[0] : 0;
      if (step > 0) expect(Math.abs(ax.ticks[0] / step % 1)).toBeLessThan(1e-6);
    }
  });

  /**
   * Once the initial transient rolls out of the window the axis is left with a
   * band a few parts in a thousand wide, which is where distinct labels stop
   * being free: at two decimals every tick in `2.4787 … 2.4817` prints as "2.48".
   * The gutter is measured from these strings, so they also have to stay short
   * enough to leave a plot behind them.
   */
  it("labels every tick distinctly at its own decimals", () => {
    for (const [lo, hi] of [[1, 2], [0.004, 0.009], [4.8, 4.9], [12, 480],
                            [2.4787, 2.4817], [112.4, 112.6]]) {
      const ax = niceAxis(lo, hi);
      const labels = ax.ticks.map((v) => v.toFixed(ax.decimals));
      expect(new Set(labels).size).toBe(labels.length);
      for (const s of labels) expect(s.length).toBeLessThanOrEqual(8);
    }
  });

  /**
   * A single sample, or a run settled to the point where f32 stops resolving the
   * change. A zero-height axis maps every value to NaN and the panel goes blank
   * exactly when the trace has reached the state it is there to show.
   */
  it("gives a flat extent a band rather than a zero span", () => {
    for (const v of [4.83, 0, -2]) {
      const ax = niceAxis(v, v);
      expect(ax.hi).toBeGreaterThan(ax.lo);
      expect(v).toBeGreaterThanOrEqual(ax.lo);
      expect(v).toBeLessThanOrEqual(ax.hi);
      expect(Number.isFinite(ax.hi - ax.lo)).toBe(true);
    }
  });

  it("survives a reversed or non-finite extent", () => {
    expect(niceAxis(6, 2)).toMatchObject({ lo: expect.any(Number) });
    expect(niceAxis(6, 2).hi).toBeGreaterThan(niceAxis(6, 2).lo);
    for (const bad of [[NaN, 1], [1, Infinity]] as const) {
      const ax = niceAxis(bad[0], bad[1]);
      expect(ax.hi).toBeGreaterThan(ax.lo);
      expect(ax.ticks).toEqual([]);
    }
  });

  it("bounds the tick count whatever the extent", () => {
    for (const [lo, hi] of [[0, 1e12], [1e-9, 2e-9], [-1e6, 1e6]])
      expect(niceAxis(lo, hi).ticks.length).toBeLessThanOrEqual(64);
  });
});

/**
 * The dimensional clock. The scaling is one multiplication, so what is worth
 * pinning is not the arithmetic but the two things that make the readout
 * trustworthy: that the *stated* reference is the one actually used — a note
 * naming R_o and κ while the conversion used something else is worse than no note
 * — and that the figures stay readable across a range that runs from Myr to Tyr.
 */
describe("dimensional time", () => {
  it("converts through the reference it prints", () => {
    // One diffusion time across the outer radius, from the reference alone.
    const expected = REFERENCE.Ro ** 2 / REFERENCE.kappa / YEAR;
    expect(TIME_UNIT_YEARS).toBeCloseTo(expected, -6);
    const note = referenceNote();
    expect(note).toContain(REFERENCE.kappa.toExponential(0));
    expect(note).toContain(`${(REFERENCE.Ro / 1e3).toFixed(0)} km`);
    expect(note).toContain(TIME_UNIT_YEARS.toExponential(2));
  });

  /**
   * The number the write-up's physical argument rests on: diffusion across the
   * mantle is orders of magnitude slower than the age of the Earth, which is why
   * it convects. If this fell to the tens of Gyr the reference would have
   * silently become something other than whole-mantle.
   */
  it("puts one diffusion time in the trillions of years", () => {
    expect(TIME_UNIT_YEARS).toBeGreaterThan(1e12);
    expect(TIME_UNIT_YEARS).toBeLessThan(2e12);
  });

  it("picks a unit that keeps three significant figures", () => {
    for (const t of [1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1, 5]) {
      const s = dimensionalTime(t);
      expect(s).toMatch(/^\d+(\.\d+)? (yr|kyr|Myr|Gyr|Tyr)$/);
      // Never four digits before the point — that is what the ladder is for.
      expect(s.split(".")[0].replace(/\D/g, "").length).toBeLessThanOrEqual(3);
    }
    expect(dimensionalTime(1)).toMatch(/Tyr$/);
    expect(dimensionalTime(1e-4)).toMatch(/Myr$/);
  });

  it("says nothing rather than something wrong for a clock that has not started", () => {
    expect(dimensionalTime(NaN)).toBe("—");
    expect(dimensionalTime(Infinity)).toBe("—");
    expect(dimensionalTime(0)).toBe("0 yr");
  });

  it("is monotone in the nondimensional time", () => {
    const ts = [1e-5, 1e-4, 1e-3, 0.01, 0.1, 1, 10];
    const yrs = ts.map((t) => (t * TIME_UNIT) / YEAR);
    for (let i = 1; i < ts.length; i++) expect(yrs[i]).toBeGreaterThan(yrs[i - 1]);
  });
});

/**
 * The x axis's two endpoint labels, and the one case the span alone gets wrong.
 * Both rows of the axis show the same two instants, so a precision that hides a
 * nonzero window start in the top row while the years row resolves it makes the
 * two rows appear to contradict each other.
 */
describe("axisDecimals", () => {
  it("resolves the window from its span", () => {
    expect(axisDecimals(9.548, 0.048)).toBe(3);   // 9.548 … 9.596
    expect(axisDecimals(0.7, 4.09)).toBe(1);      // 0.7 … 4.8, not 1 … 5
    expect(axisDecimals(4.6, 5.0)).toBe(1);
  });

  /** The real case: the first poll of a run lands a few steps in. */
  it("raises the count rather than printing a nonzero start as zero", () => {
    const td = axisDecimals(0.003, 1.235);
    expect(td).toBe(3);
    expect((0.003).toFixed(td)).toBe("0.003");
    expect(Number((0.003).toFixed(td))).not.toBe(0);
  });

  it("shares one count between the two ends", () => {
    // Both labels are formatted at the returned count — the assertion is simply
    // that one number comes back, and that it serves the smaller end.
    for (const [t0, span] of [[0.003, 1.235], [1e-5, 2], [9.5, 0.05], [0, 1]]) {
      const td = axisDecimals(t0, span);
      expect(Number.isInteger(td)).toBe(true);
      if (t0 !== 0) expect(Number(t0.toFixed(td))).not.toBe(0);
    }
  });

  it("leaves a window that genuinely starts at zero alone", () => {
    expect(axisDecimals(0, 1.235)).toBe(decimalsFor(0.1235));
    expect((0).toFixed(axisDecimals(0, 1.235))).toBe("0.0");
  });

  it("stays inside the formatter's range", () => {
    for (const [t0, span] of [[1e-30, 1], [0, 0], [NaN, 1], [1, NaN]]) {
      const td = axisDecimals(t0, span);
      expect(td).toBeGreaterThanOrEqual(0);
      expect(td).toBeLessThanOrEqual(6);
    }
  });
});

describe("decimalsFor", () => {
  it.each([[1, 0], [2, 0], [10, 0], [0.5, 1], [0.2, 1], [0.05, 2], [0.001, 3]])(
    "a step of %f prints at %i decimals", (step, want) => {
      expect(decimalsFor(step)).toBe(want);
    });

  it("clamps rather than asking for hundreds of digits", () => {
    expect(decimalsFor(1e-30)).toBe(6);
    expect(decimalsFor(0)).toBe(0);
  });
});
