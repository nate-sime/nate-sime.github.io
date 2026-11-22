// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// The model's coordinates are kilometres with the Earth's centre at the origin
// and the domain on the +z axis, so "up" is +z. The scene keeps kilometres but
// is centred on the domain and turned so up is +y, three.js's convention:
//
//     (x, y, z)  ->  (x - cx, z - cz, -(y - cy))
//
// a proper rotation (determinant +1), so handedness and normals survive.

import { Vector3 } from "three";

export type Centre = [number, number, number];

export function toScene(c: Centre, x: number, y: number, z: number, out: Float32Array, o: number): void {
  out[o] = x - c[0];
  out[o + 1] = z - c[2];
  out[o + 2] = -(y - c[1]);
}

export function vec(c: Centre, p: readonly number[]): Vector3 {
  return new Vector3(p[0] - c[0], p[2] - c[2], -(p[1] - c[1]));
}

/** A direction, which does not move with the centre. */
export function dir(p: readonly number[]): Vector3 {
  return new Vector3(p[0], p[2], -p[1]);
}

/** Convert packed (x, y, z) triples in place; NaN stays NaN. */
export function toSceneInPlace(c: Centre, a: Float32Array): void {
  for (let o = 0; o < a.length; o += 3) toScene(c, a[o], a[o + 1], a[o + 2], a, o);
}

/** The middle of the bounding box of packed (x, y, z) triples, ignoring NaN. */
export function centreOf(a: Float32Array): Centre {
  const lo = [Infinity, Infinity, Infinity];
  const hi = [-Infinity, -Infinity, -Infinity];
  for (let o = 0; o < a.length; o += 3) {
    for (let d = 0; d < 3; d++) {
      const v = a[o + d];
      if (v < lo[d]) lo[d] = v;
      if (v > hi[d]) hi[d] = v;
    }
  }
  return [(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2];
}
