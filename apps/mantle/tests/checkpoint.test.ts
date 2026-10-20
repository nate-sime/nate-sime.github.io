// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Checkpoints, without a GPU: the file format, the names, which tags fit which
 * run, what `restore` would put back, and the resampling across grids. The
 * round trip through a real solver is in `gpu.test.ts` ("checkpoints").
 */

import { describe, it, expect } from "vitest";
import {
  CHECKPOINT_FORMAT, CheckpointError, RESTORE_KEYS, changedKeys, decodeCheckpoints,
  defaultName, describe as summarise, domainMismatch, domainOf, encodeCheckpoints, resample,
  restorePatch, thumbnail, uniqueName, type Checkpoint, type TaggedState,
} from "../src/checkpoint";
import { COLORMAPS } from "../src/colormaps";
import { annulus, box } from "../src/geometry";
import { defaultState, geometryFor, type State } from "../src/ui/presets";

const GRID = { nr: 4, na: 8, gnr: 5, gna: 8 };

/** A small tag on the default annulus, with fields that exercise f32 rounding and sign. */
function tagged(patch: Partial<State> = {}, over: Partial<Checkpoint> = {}): Checkpoint {
  const state = { ...defaultState(), ...patch };
  return {
    name: "a",
    created: "2026-10-19T12:00:00.000Z",
    appVersion: "test",
    grid: GRID,
    domain: domainOf(geometryFor(state)),
    state: structuredClone(state) as unknown as TaggedState,
    clock: { time: 0.153, steps: 1234, dt: 1e-4 },
    stats: { nuInner: 3.5, nuOuter: 23.6, vrms: 41.2, tMax: 1.8, heatBalance: -0.0021 },
    T: Float32Array.from({ length: GRID.gnr * GRID.gna }, (_, k) => Math.sin(k) / 3 + 0.5),
    psi: Float32Array.from({ length: GRID.nr * GRID.na }, (_, k) => -1e3 * Math.cos(k)),
    ...over,
  };
}

const f32 = (v: number): number => Math.fround(v);

describe("checkpoint file", () => {
  it("round-trips every field bit for bit, and the metadata as data", () => {
    const a = tagged({ H: 20, logRa: 6 });
    const b = tagged({}, {
      name: "b", psi: null,
      stats: { nuInner: NaN, nuOuter: NaN, vrms: NaN, tMax: NaN, heatBalance: NaN },
    });
    const [ra, rb] = decodeCheckpoints(encodeCheckpoints([a, b]));
    expect(ra.name).toBe("a");
    expect(Array.from(ra.T)).toEqual(Array.from(a.T));
    expect(Array.from(ra.psi!)).toEqual(Array.from(a.psi!));
    expect(ra.grid).toEqual(a.grid);
    expect(ra.domain).toEqual(a.domain);
    expect(ra.clock).toEqual(a.clock);
    expect(ra.stats).toEqual(a.stats);
    expect(ra.state.H).toBe(20);
    expect(ra.state.logRa).toBe(6);
    expect(rb.psi).toBeNull();
    // JSON has no NaN: a missing reading goes out as null and comes back NaN.
    expect(Number.isNaN(rb.stats.heatBalance)).toBe(true);
  });

  it("keeps a rename, and an empty library is a valid file", () => {
    const c = tagged();
    c.name = "settled, H = 20";
    expect(decodeCheckpoints(encodeCheckpoints([c]))[0].name).toBe("settled, H = 20");
    expect(decodeCheckpoints(encodeCheckpoints([]))).toEqual([]);
  });

  it("states its format and stores little-endian f32", () => {
    const bytes = encodeCheckpoints([tagged()]);
    expect(new TextDecoder().decode(bytes.subarray(0, 8))).toBe("MANTLECK");
    const view = new DataView(bytes.buffer);
    expect(view.getUint32(8, true)).toBe(CHECKPOINT_FORMAT);
    const header = JSON.parse(new TextDecoder().decode(bytes.subarray(16, 16 + view.getUint32(12, true))));
    const start = Math.ceil((16 + view.getUint32(12, true)) / 4) * 4;
    expect(start % 4).toBe(0);
    expect(view.getFloat32(start + header.entries[0].fields.T.offset, true)).toBe(f32(0.5));
  });

  it("refuses what it cannot read, with a reason", () => {
    const good = encodeCheckpoints([tagged()]);
    const bad = (mutate: (b: Uint8Array) => Uint8Array, why: RegExp): void => {
      const b = mutate(good.slice());
      expect(() => decodeCheckpoints(b)).toThrow(CheckpointError);
      expect(() => decodeCheckpoints(b)).toThrow(why);
    };
    bad(() => new Uint8Array(64), /not a mantle checkpoint/);
    bad((b) => { new DataView(b.buffer).setUint32(8, CHECKPOINT_FORMAT + 1, true); return b; }, /newer version/);
    bad((b) => b.slice(0, b.length - 4), /truncated/);
    bad((b) => { new DataView(b.buffer).setFloat32(b.length - 4, NaN, true); return b; }, /not a number/);
    bad((b) => { b[17] = 0x7b; return b; }, /header/);   // a stray brace in the JSON
    const short = tagged({}, { T: new Float32Array(GRID.gnr * GRID.gna - 1) });
    expect(() => decodeCheckpoints(encodeCheckpoints([short]))).toThrow(/does not match its grid/);
  });
});

describe("checkpoint names", () => {
  const nameOf = (patch: Partial<State>, time = 0.153, taken: string[] = []): string => {
    const c = tagged(patch);
    return defaultName(c.state, c.domain, time, taken);
  };

  it("says where, how vigorous, how heated, and when", () => {
    expect(nameOf({ activePlanet: "earth", logRa: 6, H: 20 })).toBe("Earth · Ra 1.0e6 · H 20.0 · t 0.153");
  });

  it("leaves H out at zero, and says Ra 0 when the run is isothermal", () => {
    expect(nameOf({ activePlanet: "earth", logRa: 4, H: 0 })).toBe("Earth · Ra 1.0e4 · t 0.153");
    expect(nameOf({ activePlanet: "earth", isothermal: true, H: 0 }, 0)).toBe("Earth · Ra 0 · t 0");
  });

  it("names a box by the width it is drawn at, walled or not", () => {
    const walled = { geometry: "Cartesian box", boxLength: 2.5, walls: "free-slip walls", logRa: 4, H: 0 } as const;
    expect(nameOf(walled)).toBe("box ×2.5 · Ra 1.0e4 · t 0.153");
    expect(nameOf({ ...walled, walls: "periodic" })).toBe("box ×2.5 · Ra 1.0e4 · t 0.153");
  });

  it("names a custom planet by its own name", () => {
    const customPlanet = {
      name: "Vulcan", innerRadiusKm: 1000, outerRadiusKm: 3000, surfaceSource: "procedural",
      surface: {
        kind: "procedural", seed: 1, rockiness: 0.5, terrainScale: 0.5, oceanCoverage: 0.5,
        plantLife: 0.5, iceCaps: 0.5, cloudCover: 0.5, atmosphereHue: "#73b8ff", atmosphereDensity: 0.3,
      },
    } as const;
    expect(nameOf({ activePlanet: null, customPlanet, logRa: 5, H: 0 })).toBe("Vulcan · Ra 1.0e5 · t 0.153");
  });

  it("never repeats a name", () => {
    const base = "Earth · Ra 1.0e6 · t 0.153";
    expect(nameOf({ activePlanet: "earth", logRa: 6, H: 0 }, 0.153, [base])).toBe(`${base} (2)`);
    expect(uniqueName("a", ["a", "a (2)"])).toBe("a (3)");
    expect(uniqueName("b", ["a"])).toBe("b");
  });

  it("summarises what the name leaves out, the heat balance signed", () => {
    expect(summarise(tagged({ viscosity: "constant" })))
      .toBe("ψ 4×8 · constant · step 1,234 · heat balance −0.0021");
  });
});

describe("which tags fit which run", () => {
  const ann = domainOf(annulus());

  it("takes the same domain at any resolution, and any radial walls", () => {
    expect(domainMismatch(ann, domainOf(annulus()))).toBeNull();
    expect(domainMismatch(ann, domainOf(annulus(undefined, undefined, "no-slip")))).toBeNull();
    expect(domainMismatch(domainOf(box(2)), domainOf(box(2, "periodic", "no-slip")))).toBeNull();
  });

  it("refuses another domain, saying which way it differs", () => {
    expect(domainMismatch(ann, domainOf(box(2)))).toBe("on an annulus, not a box");
    expect(domainMismatch(ann, domainOf(annulus(0.55, 1.55)))).toMatch(/^on an annulus r 1\.21 … 2\.21, not r 0\.55/);
    expect(domainMismatch(domainOf(box(2, "free-slip")), domainOf(box(2)))).toBe("in a walled box, not a periodic one");
    expect(domainMismatch(domainOf(box(2)), domainOf(box(3)))).toBe("in a box 2 wide, not 3");
    // A walled box of width L is solved on 2L — the same period as a periodic
    // box twice as wide — and is still a different domain.
    expect(domainMismatch(domainOf(box(1, "free-slip")), domainOf(box(2)))).not.toBeNull();
  });
});

describe("restoring a tag's settings", () => {
  it("brings back every physics setting, and none of the view", () => {
    const c = tagged({ logRa: 7.5, H: 20, viscosity: "μ(T, d)", colormap: "viridis" });
    const current = { ...defaultState(), geometry: "Cartesian box" as const, H: 3, colormap: "magma" as const };
    const r = restorePatch(c, current);
    if ("problem" in r) throw new Error(r.problem);
    expect(Object.keys(r.patch).sort()).toEqual([...RESTORE_KEYS].sort());
    expect(r.patch.H).toBe(20);
    expect(r.patch.viscosity).toBe("μ(T, d)");
    expect(r.patch.geometry).toBe("spherical annulus");
    expect(changedKeys(r.patch, current)).toEqual(expect.arrayContaining(["geometry", "H", "logRa", "viscosity"]));
    expect(changedKeys(r.patch, { ...current, ...r.patch })).toEqual([]);
  });

  it("gives a setting the tag predates its default, not the pane's", () => {
    const c = tagged();
    const { H: _, ...older } = c.state;
    const r = restorePatch({ ...c, state: older }, { ...defaultState(), H: 7 });
    if ("problem" in r) throw new Error(r.problem);
    expect(r.patch.H).toBe(defaultState().H);
  });

  it("copies a custom planet rather than sharing it", () => {
    const customPlanet = {
      name: "Vulcan", innerRadiusKm: 1000, outerRadiusKm: 3000, surfaceSource: "earth-daymap",
      surface: {
        kind: "procedural", seed: 1, rockiness: 0.5, terrainScale: 0.5, oceanCoverage: 0.5,
        plantLife: 0.5, iceCaps: 0.5, cloudCover: 0.5, atmosphereHue: "#73b8ff", atmosphereDensity: 0.3,
      },
    } as const;
    const c = tagged({ activePlanet: null, customPlanet });
    const r = restorePatch(c, defaultState());
    if ("problem" in r) throw new Error(r.problem);
    expect(r.patch.customPlanet).toEqual(customPlanet);
    expect(r.patch.customPlanet).not.toBe(c.state.customPlanet);
  });

  it("refuses a setting this version does not offer", () => {
    const c = tagged();
    const r = restorePatch({ ...c, state: { ...c.state, resolution: "enormous" } }, defaultState());
    expect("problem" in r && r.problem).toMatch(/resolution/);
  });

  it("refuses settings that no longer reproduce the field's domain", () => {
    const c = tagged({ geometry: "Cartesian box", boxLength: 3, walls: "periodic" },
      { domain: domainOf(box(2)) });
    const r = restorePatch(c, defaultState());
    expect("problem" in r && r.problem).toMatch(/no longer give the domain.*2 wide, not 3/);
  });
});

describe("resampling across grids", () => {
  const field = (nr: number, na: number, f: (s: number, phi: number) => number): Float32Array =>
    Float32Array.from({ length: nr * na }, (_, k) =>
      f(Math.floor(k / na) / (nr - 1), (2 * Math.PI * (k % na)) / na));

  it("returns a same-sized grid unchanged, as a copy", () => {
    const src = field(5, 8, (s, p) => s + Math.cos(p));
    const out = resample(src, 5, 8, 5, 8);
    expect(Array.from(out)).toEqual(Array.from(src));
    expect(out).not.toBe(src);
  });

  it("lands exactly on the source wherever the nodes coincide, boundary rows included", () => {
    // Exactly the Dirichlet data on the two boundary rows, as a run's T is.
    const src = field(5, 8, (s, p) => (1 - s) + 0.2 * s * (1 - s) * Math.cos(2 * p));
    const out = resample(src, 5, 8, 9, 16);
    for (let i = 0; i < 5; i++)
      for (let j = 0; j < 8; j++) expect(out[2 * i * 16 + 2 * j]).toBe(src[i * 8 + j]);
    for (let j = 0; j < 16; j++) {
      expect(out[j]).toBe(src[0]);
      expect(out[8 * 16 + j]).toBe(src[4 * 8]);
    }
  });

  it("creates no new extrema", () => {
    let seed = 7;
    const rand = (): number => ((seed = (seed * 16807) % 2147483647) / 2147483647);
    const src = Float32Array.from({ length: 9 * 16 }, rand);
    const lo = Math.min(...src), hi = Math.max(...src);
    for (const [nr, na] of [[17, 32], [5, 8], [13, 24]]) {
      const out = resample(src, 9, 16, nr, na);
      expect(Math.min(...out)).toBeGreaterThanOrEqual(lo);
      expect(Math.max(...out)).toBeLessThanOrEqual(hi);
    }
  });

  // Second order, boundary rows included: the linear ghost row is what keeps
  // the first cell in from the ends from dropping to first order.
  it("converges at second order on a smooth field as the source is refined", () => {
    const f = (s: number, p: number): number => s * s + 0.3 * Math.cos(p) * Math.sin(Math.PI * s);
    const exact = field(97, 192, f);
    const err = (nr: number, na: number): number => {
      const out = resample(field(nr, na, f), nr, na, 97, 192);
      return Math.max(...out.map((v, k) => Math.abs(v - exact[k])));
    };
    const coarse = err(13, 24), fine = err(25, 48);
    expect(fine).toBeLessThan(1.5e-3);
    expect(coarse / fine).toBeGreaterThan(3.5);
  });
});

describe("thumbnails", () => {
  const alpha = (px: Uint8ClampedArray, size: number, x: number, y: number): number => px[4 * (y * size + x) + 3];

  it("draws the annulus as a ring, transparent inside and out", () => {
    const c = tagged();
    const N = 32;
    const px = thumbnail(c, N, COLORMAPS.inferno);
    expect(alpha(px, N, N / 2, N / 2)).toBe(0);   // the core
    expect(alpha(px, N, 0, 0)).toBe(0);           // outside the planet
    expect(alpha(px, N, N - 2, N / 2)).toBe(255); // in the mantle, near the surface
  });

  it("puts a box's hot boundary at the bottom", () => {
    const g = box(1, "periodic");
    const gnr = 5, gna = 8;
    const T = Float32Array.from({ length: gnr * gna }, (_, k) => 1 - Math.floor(k / gna) / (gnr - 1));
    const N = 16;
    const stops = COLORMAPS.inferno;
    const px = thumbnail({ T, grid: { nr: 4, na: 8, gnr, gna }, domain: domainOf(g) }, N, stops);
    const rgb = (y: number): number[] => Array.from(px.subarray(4 * y * N, 4 * y * N + 3));
    expect(rgb(N - 1)).toEqual(stops[stops.length - 1].map((v) => Math.round(255 * v)));
    expect(rgb(0)).toEqual(stops[0].map((v) => Math.round(255 * v)));
  });
});
