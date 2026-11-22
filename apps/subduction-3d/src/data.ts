// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// Loading what tools/export_data.py writes into public/data/, and the pure
// bookkeeping on it (frame differencing, keyframe brackets, tracer trails) —
// kept free of three.js so the tests can exercise it directly.

export type Vec3 = [number, number, number];

export interface ArrayFile {
  file: string;
  dtype: "uint8" | "int8" | "uint16" | "float32";
  shape: number[];
  bytes: number;
  /** Frame k is stored as frame k minus frame k-1, wrapping. */
  delta?: boolean;
}

export interface Manifest {
  source: string;
  steps: number;
  myrPerStep: number;
  keyframes: number[];
  temperatureRange: [number, number];
  surface: {
    nx: number;
    ny: number;
    origin: [number, number];
    spacing: number;
    R0: number;
    radius: ArrayFile & { range: [number, number] };
    temperature: ArrayFile;
    edgeDistance: ArrayFile & { kmPerUnit: number };
  };
  slices: ArrayFile & { y: number[]; nx: number; nz: number; x: [number, number]; z: [number, number] };
  outline: ArrayFile & { counts: number[] };
  tracers: ArrayFile;
  camera: { position: Vec3; target: Vec3; up: Vec3; fov: number };
}

export interface Dataset {
  manifest: Manifest;
  /** [keyframe][ny][nx], quantised over manifest.surface.radius.range. */
  radius: Uint16Array;
  /** [keyframe][ny][nx], quantised over manifest.temperatureRange. */
  slabTemperature: Uint8Array;
  /** [keyframe][ny][nx], signed distance to the slab's edge in kmPerUnit. */
  edge: Int8Array;
  /** [keyframe][slice][nz][nx]: 0 outside the domain, 1..255 the temperature. */
  slices: Uint8Array;
  /** Segments of the domain outline, every step back to back: (x, y, z) × 2 each. */
  outline: Float32Array;
  /** outlineStart[step] is the first segment of that step; [steps] is the total. */
  outlineStart: Int32Array;
  /** [tracer][step][3], NaN when that tracer does not exist at that step. */
  tracers: Float32Array;
  tracerCount: number;
}

/** Inflate a gzip body, unless something on the way (a dev server) already did. */
export async function gunzip(buf: ArrayBuffer): Promise<ArrayBuffer> {
  const b = new Uint8Array(buf);
  if (b.length < 2 || b[0] !== 0x1f || b[1] !== 0x8b) return buf;
  const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Response(stream).arrayBuffer();
}

/** Undo the frame differencing in place: frame k += frame k-1. The typed-array
 *  store wraps modulo 2^8 or 2^16 just as the exporter's subtraction did. */
export function undelta(a: Uint8Array | Uint16Array, frameSize: number): void {
  for (let i = frameSize; i < a.length; i++) a[i] = a[i] + a[i - frameSize];
}

export function decode(buf: ArrayBuffer, f: ArrayFile): Uint8Array | Int8Array | Uint16Array | Float32Array {
  const n = f.shape.reduce((a, b) => a * b, 1);
  const frame = n / f.shape[0];
  switch (f.dtype) {
    case "uint16": {
      const a = new Uint16Array(buf, 0, n);
      if (f.delta) undelta(a, frame);
      return a;
    }
    case "uint8":
    case "int8": {
      const a = new Uint8Array(buf, 0, n);
      if (f.delta) undelta(a, frame);
      return f.dtype === "int8" ? new Int8Array(buf, 0, n) : a;
    }
    case "float32":
      return new Float32Array(buf, 0, n);
  }
}

export async function load(base: string, onProgress: (fraction: number) => void): Promise<Dataset> {
  const res = await fetch(base + "manifest.json");
  if (!res.ok) throw new Error(`manifest.json: HTTP ${res.status}`);
  const m = (await res.json()) as Manifest;
  const files = [m.surface.radius, m.surface.temperature, m.surface.edgeDistance, m.slices, m.outline, m.tracers];
  const total = files.reduce((s, f) => s + f.bytes, 0);
  let done = 0;
  onProgress(0);
  const arrays = await Promise.all(
    files.map(async (f) => {
      const r = await fetch(base + f.file);
      if (!r.ok) throw new Error(`${f.file}: HTTP ${r.status}`);
      const raw = await r.arrayBuffer();
      done += f.bytes;
      onProgress(done / total);
      return decode(await gunzip(raw), f);
    }),
  );
  const outlineStart = new Int32Array(m.steps + 1);
  for (let s = 0; s < m.steps; s++) outlineStart[s + 1] = outlineStart[s] + m.outline.counts[s];
  return {
    manifest: m,
    radius: arrays[0] as Uint16Array,
    slabTemperature: arrays[1] as Uint8Array,
    edge: arrays[2] as Int8Array,
    slices: arrays[3] as Uint8Array,
    outline: arrays[4] as Float32Array,
    outlineStart,
    tracers: arrays[5] as Float32Array,
    tracerCount: m.tracers.shape[0],
  };
}

/** The keyframe interval holding step t: keyframes[k] <= t <= keyframes[k + 1],
 *  with s the weight of the later one. */
export function bracket(keys: readonly number[], t: number): { k: number; s: number } {
  const last = keys.length - 1;
  if (t <= keys[0]) return { k: 0, s: 0 };
  if (t >= keys[last]) return { k: last - 1, s: 1 };
  let lo = 0;
  let hi = last;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (keys[mid] <= t) lo = mid;
    else hi = mid;
  }
  return { k: lo, s: (t - keys[lo]) / (keys[lo + 1] - keys[lo]) };
}

/** Position of tracer `i` at fractional step t, linear between steps, written
 *  to out[o..o+2]; false if it does not exist then (not yet seeded, or gone). */
export function tracerAt(P: Float32Array, steps: number, i: number, t: number, out: Float32Array, o: number): boolean {
  if (t < 0 || t > steps - 1) return false;
  const s0 = Math.floor(t);
  const f = t - s0;
  const a = (i * steps + s0) * 3;
  if (Number.isNaN(P[a])) return false;
  if (f === 0) {
    out[o] = P[a];
    out[o + 1] = P[a + 1];
    out[o + 2] = P[a + 2];
    return true;
  }
  const b = a + 3;
  if (Number.isNaN(P[b])) return false;
  out[o] = P[a] + (P[b] - P[a]) * f;
  out[o + 1] = P[a + 1] + (P[b + 1] - P[a + 1]) * f;
  out[o + 2] = P[a + 2] + (P[b + 2] - P[a + 2]) * f;
  return true;
}

/** The pathline of tracer `i` over the `length` steps up to t, newest point
 *  first, cut short where the tracer was seeded. Writes points to out from
 *  index o (3 floats each) and returns how many; 0 if the tracer is not alive
 *  at t. */
export function trail(
  P: Float32Array, steps: number, i: number, t: number, length: number, out: Float32Array, o: number,
): number {
  if (!tracerAt(P, steps, i, t, out, o)) return 0;
  let n = 1;
  const end = t - length;
  for (let s = Math.ceil(t) - 1; s >= 0 && s > end; s--) {
    const a = (i * steps + s) * 3;
    if (Number.isNaN(P[a])) return n;
    out[o + 3 * n] = P[a];
    out[o + 3 * n + 1] = P[a + 1];
    out[o + 3 * n + 2] = P[a + 2];
    n++;
  }
  if (end > 0 && length > 0 && tracerAt(P, steps, i, end, out, o + 3 * n)) n++;
  return n;
}

/**
 * Slice texels as the renderer wants them: two channels per texel.
 *
 *  - R: the temperature, carried outward past the domain's edge from the
 *    nearest texel inside, so blending two frames whose domains differ, or
 *    filtering across the edge, never mixes in the "outside" code.
 *  - G: a signed distance to the domain's edge, 128 + 8 × texels (positive
 *    inside), clamped. Blended between frames, its zero contour moves
 *    smoothly as the domain grows instead of jumping a keyframe at a time.
 *
 * Distances are a two-pass 3×3 chamfer (1, √2) — exact enough for an edge.
 * This runs over every texel of every keyframe (~7M) as the page opens, so
 * the loops are kept plain: one distance field per texel (an inside texel's
 * distance to the outside; an outside texel's to the inside — the other is
 * 0 and cannot change), neighbours unrolled, no allocation.
 */
export function prepareSlices(q: Uint8Array, nx: number, nz: number, layers: number): Uint8Array {
  const out = new Uint8Array(2 * nx * nz * layers);
  const N = nx * nz;
  const dist = new Float32Array(N);
  const val = new Uint8Array(N); // the temperature code, carried outward
  const D = Math.SQRT2;
  let L = q;
  // Relax texel p from neighbour r at distance w. A neighbour on the other
  // side of the edge is a seed, at distance 0; outside texels also take the
  // temperature of whichever texel gave them their distance.
  const relax = (p: number, r: number, w: number): void => {
    const d = (L[p] === 0) !== (L[r] === 0) ? w : dist[r] + w;
    if (d < dist[p]) {
      dist[p] = d;
      if (L[p] === 0) val[p] = val[r];
    }
  };
  for (let l = 0; l < layers; l++) {
    L = q.subarray(l * N, (l + 1) * N);
    for (let p = 0; p < N; p++) {
      dist[p] = 1e9;
      val[p] = L[p];
    }
    for (let j = 0; j < nz; j++) {
      for (let i = 0, p = j * nx; i < nx; i++, p++) {
        if (i > 0) relax(p, p - 1, 1);
        if (j > 0) {
          relax(p, p - nx, 1);
          if (i > 0) relax(p, p - nx - 1, D);
          if (i < nx - 1) relax(p, p - nx + 1, D);
        }
      }
    }
    for (let j = nz - 1; j >= 0; j--) {
      for (let i = nx - 1, p = j * nx + nx - 1; i >= 0; i--, p--) {
        if (i < nx - 1) relax(p, p + 1, 1);
        if (j < nz - 1) {
          relax(p, p + nx, 1);
          if (i < nx - 1) relax(p, p + nx + 1, D);
          if (i > 0) relax(p, p + nx - 1, D);
        }
      }
    }
    const o = l * N * 2;
    for (let p = 0; p < N; p++) {
      const sd = L[p] !== 0 ? dist[p] - 0.5 : 0.5 - dist[p];
      out[o + 2 * p] = Math.max(1, val[p]);
      out[o + 2 * p + 1] = Math.max(0, Math.min(255, Math.round(128 + 8 * sd)));
    }
  }
  return out;
}
