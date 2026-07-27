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
  NU_COLOUR, NuTrace, decimalsFor, niceAxis, niceStep, type NuSeries,
} from "../src/ui/nusselt";

const series = Object.keys(NU_COLOUR) as NuSeries[];

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
    expect(tr.push(0.1, 2, 3)).toBe(true);
    expect(tr.push(0.2, 4, 5)).toBe(true);
    expect(tr.length).toBe(2);
    expect(tr.at(0)).toEqual({ t: 0.1, inner: 2, outer: 3 });
    expect(tr.at(1)).toEqual({ t: 0.2, inner: 4, outer: 5 });
    expect(tr.last).toEqual({ t: 0.2, inner: 4, outer: 5 });
  });

  it("rolls the oldest out once full, and stays ordered across the wrap", () => {
    const tr = new NuTrace(4);
    for (let i = 0; i < 10; i++) tr.push(i, i, -i);
    expect(tr.length).toBe(4);
    expect(Array.from({ length: 4 }, (_, i) => tr.at(i).t)).toEqual([6, 7, 8, 9]);
    expect(tr.extent()).toMatchObject({ t0: 6, t1: 9, lo: -9, hi: 9 });
  });

  // The readout prints a dash until the first readback lands; the plot must not
  // plot the NaN behind it.
  it("drops non-finite samples", () => {
    const tr = new NuTrace(8);
    expect(tr.push(NaN, NaN, NaN)).toBe(false);
    expect(tr.push(0.1, NaN, 3)).toBe(false);
    expect(tr.push(0.1, 2, Infinity)).toBe(false);
    expect(tr.length).toBe(0);
  });

  // The poll is asynchronous and allowed to be one behind, so the frame loop
  // offers the same reading many times over.
  it("drops a repeat of the sample it already holds", () => {
    const tr = new NuTrace(8);
    expect(tr.push(0.1, 2, 3)).toBe(true);
    expect(tr.push(0.1, 2, 3)).toBe(false);
    expect(tr.push(0.1, 9, 9)).toBe(false);
    expect(tr.length).toBe(1);
  });

  /**
   * The rule the whole class exists for. Reseed sets the clock to zero, and a
   * poll already in flight can land after it — so the guard has to be here,
   * where every sample passes, and not at the call sites that know a reset
   * happened.
   */
  it("restarts the window when the clock goes backwards", () => {
    const tr = new NuTrace(8);
    for (const t of [0.1, 0.2, 0.3]) tr.push(t, 5, 5);
    expect(tr.push(0.0, 1, 1)).toBe(true);
    expect(tr.length).toBe(1);
    expect(tr.at(0)).toEqual({ t: 0, inner: 1, outer: 1 });
  });

  it("refuses to index outside the window", () => {
    const tr = new NuTrace(4);
    tr.push(1, 1, 1);
    expect(() => tr.at(1)).toThrow(RangeError);
    expect(() => tr.at(-1)).toThrow(RangeError);
  });

  it("rejects a capacity that cannot hold a line", () => {
    expect(() => new NuTrace(1)).toThrow();
    expect(() => new NuTrace(2.5)).toThrow();
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
