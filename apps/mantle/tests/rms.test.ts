/**
 * The RMS velocity trace's data layer — `RmsTrace` is `NuTrace` with one
 * series, so what is worth re-checking here is only what a single series
 * changes: the sample shape, and that the colour is distinct from both
 * Nusselt series it shares a page with. The ring buffer's restart rule and the
 * axis arithmetic are already covered in `nusselt.test.ts`, against the same
 * functions this module re-exports rather than re-derives.
 */

import { describe, it, expect } from "vitest";
import { RMS_COLOUR, RmsTrace } from "../src/ui/rms";
import { NU_COLOUR } from "../src/ui/nusselt";

describe("RMS_COLOUR", () => {
  it("is an opaque hex the canvas and the legend can share", () => {
    expect(RMS_COLOUR).toMatch(/^#[0-9a-f]{6}$/);
  });

  it("is distinct from both Nusselt series it shares a page with", () => {
    expect(RMS_COLOUR).not.toBe(NU_COLOUR.inner);
    expect(RMS_COLOUR).not.toBe(NU_COLOUR.outer);
  });
});

describe("RmsTrace", () => {
  it("starts empty and reports nothing to scale to", () => {
    const tr = new RmsTrace(8);
    expect(tr.length).toBe(0);
    expect(tr.last).toBeNull();
    expect(tr.extent()).toBeNull();
  });

  it("holds samples oldest-first", () => {
    const tr = new RmsTrace(8);
    expect(tr.push({ t: 0.1, step: 10, v: 2 })).toBe(true);
    expect(tr.push({ t: 0.2, step: 20, v: 4 })).toBe(true);
    expect(tr.length).toBe(2);
    expect(tr.at(0)).toEqual({ t: 0.1, step: 10, v: 2 });
    expect(tr.last).toEqual({ t: 0.2, step: 20, v: 4 });
  });

  it("rolls the oldest out once full, and stays ordered across the wrap", () => {
    const tr = new RmsTrace(4);
    for (let i = 0; i < 10; i++) tr.push({ t: i, step: i, v: i });
    expect(tr.length).toBe(4);
    expect(Array.from({ length: 4 }, (_, i) => tr.at(i).t)).toEqual([6, 7, 8, 9]);
    expect(tr.extent()).toMatchObject({ t0: 6, t1: 9, lo: 6, hi: 9 });
  });

  it("drops non-finite samples", () => {
    const tr = new RmsTrace(8);
    expect(tr.push({ t: NaN, step: NaN, v: NaN })).toBe(false);
    expect(tr.push({ t: 0.1, step: 1, v: Infinity })).toBe(false);
    expect(tr.push({ t: 0.1, step: NaN, v: 2 })).toBe(false);
    expect(tr.length).toBe(0);
  });

  it("drops a repeat of the sample it already holds", () => {
    const tr = new RmsTrace(8);
    expect(tr.push({ t: 0.1, step: 1, v: 2 })).toBe(true);
    expect(tr.push({ t: 0.1, step: 1, v: 2 })).toBe(false);
    expect(tr.length).toBe(1);
  });

  it("restarts the buffer when the clock goes backwards", () => {
    const tr = new RmsTrace(8);
    for (const t of [0.1, 0.2, 0.3]) tr.push({ t, step: t * 1e4, v: 5 });
    expect(tr.push({ t: 0, step: 0, v: 1 })).toBe(true);
    expect(tr.length).toBe(1);
    expect(tr.at(0)).toMatchObject({ t: 0, v: 1 });
  });

  it("refuses to index outside the buffer", () => {
    const tr = new RmsTrace(4);
    tr.push({ t: 1, step: 1, v: 1 });
    expect(() => tr.at(1)).toThrow(RangeError);
    expect(() => tr.at(-1)).toThrow(RangeError);
  });

  it("rejects a capacity that cannot hold a line", () => {
    expect(() => new RmsTrace(1)).toThrow();
    expect(() => new RmsTrace(2.5)).toThrow();
  });

  it("keeps exactly the samples within the window of the newest", () => {
    const tr = new RmsTrace(64);
    for (let i = 0; i < 20; i++) tr.push({ t: i * 0.05, step: i * 50, v: 2 });
    expect(tr.first(200)).toBe(15);
    expect(tr.at(tr.first(200)).step).toBe(750);
    expect(tr.first(Infinity)).toBe(0);
  });
});
