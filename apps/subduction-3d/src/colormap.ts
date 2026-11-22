// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { ClampToEdgeWrapping, DataTexture, LinearFilter, RGBAFormat, SRGBColorSpace } from "three";

/** [position 0..1, r, g, b] stops, sRGB. */
export type Stops = readonly (readonly [number, number, number, number])[];

// ParaView's "Reds" as the published figure used it, over 272.87..1575.93 K
// (statev2temp.pvsm). Evenly spaced, so the positions are just i / 16.
const REDS_RGB: readonly (readonly [number, number, number])[] = [
  [1.0, 0.960784, 0.941176], [0.998155, 0.92203, 0.885813], [0.996186, 0.880692, 0.826759],
  [0.992357, 0.809581, 0.732349], [0.98842, 0.736747, 0.635894], [0.988235, 0.656409, 0.543191],
  [0.988235, 0.575702, 0.450673], [0.986344, 0.496886, 0.371396], [0.984375, 0.418147, 0.292657],
  [0.96143, 0.326059, 0.232987], [0.937809, 0.233541, 0.173933], [0.868051, 0.164091, 0.143714],
  [0.797186, 0.095194, 0.114187], [0.722445, 0.076678, 0.098224], [0.647643, 0.058962, 0.082476],
  [0.525967, 0.029527, 0.066728], [0.403922, 0.0, 0.05098],
];
export const REDS: Stops = REDS_RGB.map((c, i) => [i / (REDS_RGB.length - 1), c[0], c[1], c[2]] as const);

// The yellow-green-to-ink ramp the figure gave its depth contours, shallow to
// deep: light near the surface, dark where the slab is deepest.
const DEPTH_RGB: readonly (readonly [number, number, number])[] = [
  [1.0, 0.984314, 0.901961], [0.960784, 0.94902, 0.670588], [0.886275, 0.921569, 0.533333],
  [0.784314, 0.878431, 0.396078], [0.666667, 0.839216, 0.294118], [0.556863, 0.8, 0.239216],
  [0.431373, 0.760784, 0.160784], [0.317647, 0.721569, 0.113725], [0.211765, 0.678431, 0.082353],
  [0.109804, 0.631373, 0.05098], [0.082353, 0.588235, 0.082353], [0.109804, 0.54902, 0.152941],
  [0.113725, 0.521569, 0.203922], [0.117647, 0.490196, 0.243137], [0.117647, 0.45098, 0.270588],
  [0.113725, 0.4, 0.278431], [0.109804, 0.34902, 0.278431], [0.094118, 0.278431, 0.25098],
  [0.086275, 0.231373, 0.219608], [0.07451, 0.172549, 0.180392], [0.054902, 0.109804, 0.121569],
];
export const DEPTH: Stops = DEPTH_RGB.map((c, i) => [i / (DEPTH_RGB.length - 1), c[0], c[1], c[2]] as const);

export function sample(stops: Stops, f: number): [number, number, number] {
  const x = Math.min(1, Math.max(0, f));
  let i = 1;
  while (i < stops.length - 1 && stops[i][0] < x) i++;
  const a = stops[i - 1];
  const b = stops[i];
  const w = (x - a[0]) / (b[0] - a[0] || 1);
  return [a[1] + (b[1] - a[1]) * w, a[2] + (b[2] - a[2]) * w, a[3] + (b[3] - a[3]) * w];
}

/** A 256×1 lookup texture; sRGB-tagged, so shaders read it back linear. */
export function lutTexture(stops: Stops): DataTexture {
  const data = new Uint8Array(256 * 4);
  for (let i = 0; i < 256; i++) {
    const c = sample(stops, i / 255);
    data[4 * i] = Math.round(c[0] * 255);
    data[4 * i + 1] = Math.round(c[1] * 255);
    data[4 * i + 2] = Math.round(c[2] * 255);
    data[4 * i + 3] = 255;
  }
  const t = new DataTexture(data, 256, 1, RGBAFormat);
  t.colorSpace = SRGBColorSpace;
  t.minFilter = t.magFilter = LinearFilter;
  t.wrapS = t.wrapT = ClampToEdgeWrapping;
  t.needsUpdate = true;
  return t;
}

export function cssGradient(stops: Stops, direction = "to right"): string {
  const parts = [];
  for (let i = 0; i <= 16; i++) {
    const c = sample(stops, i / 16).map((v) => Math.round(v * 255));
    parts.push(`rgb(${c[0]} ${c[1]} ${c[2]}) ${((i / 16) * 100).toFixed(1)}%`);
  }
  return `linear-gradient(${direction}, ${parts.join(", ")})`;
}
