// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { describe, expect, it } from "vitest";
import { bracket, decode, type Manifest, prepareSlices, tracerAt, trail, undelta } from "../src/data";

/** The exporter's encoding, for round trips: frame k minus frame k-1, wrapping. */
function delta<T extends Uint8Array | Uint16Array>(a: T, frame: number): T {
  const d = a.slice() as T;
  for (let i = a.length - 1; i >= frame; i--) d[i] = a[i] - a[i - frame];
  return d;
}

describe("frame differencing", () => {
  it("round-trips uint8, uint16 and int8 (through its bytes), wrapping", () => {
    const u8 = Uint8Array.from([250, 3, 7, 4, 255, 0, 1, 254, 128]);
    const d8 = delta(u8, 3);
    undelta(d8, 3);
    expect([...d8]).toEqual([...u8]);

    const u16 = Uint16Array.from([65535, 0, 12345, 1, 65534, 40000]);
    const d16 = delta(u16, 2);
    undelta(d16, 2);
    expect([...d16]).toEqual([...u16]);

    const i8 = Int8Array.from([-127, 5, 127, -3, 0, -128]);
    const bytes = delta(new Uint8Array(i8.buffer.slice(0)), 2);
    const back = decode(bytes.buffer, { file: "", dtype: "int8", shape: [3, 2], bytes: 0, delta: true });
    expect([...back]).toEqual([...i8]);
  });
});

describe("bracket", () => {
  const keys = [0, 3, 6, 9, 10];
  it("finds the keyframe interval and the weight of its later end", () => {
    expect(bracket(keys, 0)).toEqual({ k: 0, s: 0 });
    expect(bracket(keys, 1.5)).toEqual({ k: 0, s: 0.5 });
    expect(bracket(keys, 3)).toEqual({ k: 1, s: 0 });
    expect(bracket(keys, 9.5)).toEqual({ k: 3, s: 0.5 });
    expect(bracket(keys, 10)).toEqual({ k: 3, s: 1 });
    expect(bracket(keys, 12)).toEqual({ k: 3, s: 1 });
    expect(bracket(keys, -1)).toEqual({ k: 0, s: 0 });
  });
});

describe("tracers", () => {
  // One tracer over 6 steps, seeded at step 1 and gone after step 4,
  // moving +1 in x per step.
  const steps = 6;
  const P = new Float32Array(3 * steps).fill(NaN);
  for (let s = 1; s <= 4; s++) P.set([s, 0, 0], 3 * s);
  const out = new Float32Array(64);
  const xs = (n: number) => Array.from({ length: n }, (_, q) => out[3 * q]);

  it("is linear between steps and absent outside its life", () => {
    expect(tracerAt(P, steps, 0, 2.25, out, 0)).toBe(true);
    expect(out[0]).toBeCloseTo(2.25);
    expect(tracerAt(P, steps, 0, 0.5, out, 0)).toBe(false); // not yet seeded
    expect(tracerAt(P, steps, 0, 4, out, 0)).toBe(true);
    expect(tracerAt(P, steps, 0, 4.5, out, 0)).toBe(false); // left the domain
  });

  it("trails back `length` steps, newest first, cut at its seeding", () => {
    expect(xs(trail(P, steps, 0, 3.5, 1.5, out, 0))).toEqual([3.5, 3, 2]);
    expect(xs(trail(P, steps, 0, 3.5, 1.25, out, 0))).toEqual([3.5, 3, 2.25]);
    expect(xs(trail(P, steps, 0, 3, 1, out, 0))).toEqual([3, 2]);
    expect(xs(trail(P, steps, 0, 3.5, 10, out, 0))).toEqual([3.5, 3, 2, 1]);
    expect(xs(trail(P, steps, 0, 2.5, 0, out, 0))).toEqual([2.5]);
    expect(trail(P, steps, 0, 5, 2, out, 0)).toBe(0);
  });
});

describe("prepareSlices", () => {
  it("signs the edge distance and carries the temperature outward", () => {
    // 6 × 1 strip: inside (codes 10, 20, 30) then outside.
    const q = Uint8Array.from([10, 20, 30, 0, 0, 0]);
    const rg = prepareSlices(q, 6, 1, 1);
    const R = [0, 1, 2, 3, 4, 5].map((i) => rg[2 * i]);
    const G = [0, 1, 2, 3, 4, 5].map((i) => rg[2 * i + 1]);
    expect(R).toEqual([10, 20, 30, 30, 30, 30]);
    expect(G[2]).toBe(128 + 4); // half a texel inside
    expect(G[3]).toBe(128 - 4); // half a texel outside
    expect(G[2] > G[1] === false).toBe(true); // deeper inside is no closer to the edge
    expect(G[5]).toBeLessThan(G[4]);
  });
});

describe("exported data", () => {
  const dir = new URL("../public/data/", import.meta.url);
  const m = JSON.parse(readFileSync(new URL("manifest.json", dir), "utf-8")) as Manifest;
  const read = (f: { file: string }) => {
    const b = gunzipSync(readFileSync(new URL(f.file, dir)));
    return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) as ArrayBuffer;
  };

  it("matches its manifest", () => {
    const K = m.keyframes.length;
    expect(m.keyframes[0]).toBe(0);
    expect(m.keyframes[K - 1]).toBe(m.steps - 1);
    const s = m.surface;
    for (const f of [s.radius, s.temperature, s.edgeDistance]) expect(f.shape).toEqual([K, s.ny, s.nx]);
    expect(m.slices.shape).toEqual([K, m.slices.y.length, m.slices.nz, m.slices.nx]);
    expect(m.tracers.shape.slice(1)).toEqual([m.steps, 3]);
    expect(m.outline.counts.length).toBe(m.steps);
    expect(m.outline.shape).toEqual([m.outline.counts.reduce((a, b) => a + b, 0), 2, 3]);
  });

  it("decodes to a slab that grows, and is smooth wherever it is drawn", () => {
    const s = m.surface;
    const { nx, ny } = s;
    const N = nx * ny;
    expect(s.edgeDistance.dtype).toBe("int8"); // signed: negative is off the slab
    const R = decode(read(s.radius), s.radius) as Uint16Array;
    const E = decode(read(s.edgeDistance), s.edgeDistance);
    expect(E).toBeInstanceOf(Int8Array);
    const [lo, hi] = s.radius.range;
    const r = (k: number, p: number) => lo + (R[k * N + p] * (hi - lo)) / 65535;
    const counts: number[] = [];
    for (let k = 0; k < m.keyframes.length; k++) {
      let on = 0;
      let jump = 0;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const p = j * nx + i;
          if (E[k * N + p] >= 0) on++;
          // Any edge with an end on the slab is (partly) drawn: neighbours
          // there must not be far apart in radius, or the mesh shows a wall.
          for (const q of [i + 1 < nx ? p + 1 : -1, j + 1 < ny ? p + nx : -1]) {
            if (q < 0 || (E[k * N + p] < 0 && E[k * N + q] < 0)) continue;
            jump = Math.max(jump, Math.abs(r(k, p) - r(k, q)));
          }
        }
      }
      expect(jump).toBeLessThan(8); // a 4 km cell, dipping at most ~52°
      if (k > 0) expect(on).toBeGreaterThanOrEqual(counts[k - 1]); // the slab top only lengthens
      counts.push(on);
    }
    // From well under half the grid to most of it (16k -> 35k of 42k cells).
    expect(counts[0]).toBeLessThan(0.5 * N);
    expect(counts[counts.length - 1]).toBeGreaterThan(0.75 * N);
  });

  it("decodes tracers seeded on the published plane at x = 0, every 25 steps", () => {
    const T = decode(read(m.tracers), m.tracers) as Float32Array;
    const n = m.tracers.shape[0];
    let atZero = 0;
    for (let i = 0; i < n; i++) {
      // The first step each tracer exists is a seed point, on x = 0.
      for (let s = 0; s < m.steps; s++) {
        const a = (i * m.steps + s) * 3;
        if (Number.isNaN(T[a])) continue;
        expect(Math.abs(T[a])).toBeLessThan(1e-3);
        expect(s % 25).toBe(0);
        atZero++;
        break;
      }
    }
    expect(atZero).toBe(n);
    expect(n).toBe(4 * 77);
  });
});
