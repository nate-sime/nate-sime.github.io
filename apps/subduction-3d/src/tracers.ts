// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// The tracers: a sphere at each one's position now and its pathline over the
// last `tail` steps, both linear between the exporter's per-step positions.
// The tail fades toward the page's white with age, which reads as a fade
// without the sorting a transparent line would need.

import { Color, InstancedMesh, Matrix4, MeshStandardMaterial, SphereGeometry } from "three";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";
import { type Dataset, trail } from "./data";

/** ParaView's (0, 2/3, 1), the published animation's tracer colour. */
export const TRACER_COLOUR = new Color().setRGB(0, 2 / 3, 1);

export class Tracers {
  readonly spheres: InstancedMesh;
  readonly tails: LineSegments2;
  private readonly tailPos: Float32Array;
  private readonly tailCol: Float32Array;
  private readonly points: Float32Array;
  private readonly matrix = new Matrix4();

  constructor(
    private readonly data: Dataset,
    radius: number,
    private readonly maxTail: number,
  ) {
    const n = data.tracerCount;
    const mat = new MeshStandardMaterial({
      color: TRACER_COLOUR,
      emissive: TRACER_COLOUR,
      emissiveIntensity: 0.55,
      roughness: 0.35,
    });
    this.spheres = new InstancedMesh(new SphereGeometry(radius, 16, 12), mat, n);
    this.spheres.frustumCulled = false;
    this.spheres.count = 0;

    const maxSegments = n * (maxTail + 2);
    this.tailPos = new Float32Array(6 * maxSegments);
    this.tailCol = new Float32Array(6 * maxSegments);
    this.points = new Float32Array(3 * (maxTail + 3));
    const g = new LineSegmentsGeometry();
    g.setPositions(this.tailPos);
    g.setColors(this.tailCol);
    g.instanceCount = 0;
    const lm = new LineMaterial({ linewidth: 2.5, vertexColors: true });
    lm.toneMapped = false;
    this.tails = new LineSegments2(g, lm);
    this.tails.frustumCulled = false;
  }

  /** Show the tracers at fractional step t, with tails `tail` steps long. */
  update(t: number, tail: number): void {
    const d = this.data;
    const steps = d.manifest.steps;
    const P = d.tracers;
    const pts = this.points;
    const L = Math.min(tail, this.maxTail);
    const c0 = TRACER_COLOUR;
    let shown = 0;
    let seg = 0;
    for (let i = 0; i < d.tracerCount; i++) {
      const n = trail(P, steps, i, t, L, pts, 0);
      if (n === 0) continue;
      this.matrix.makeTranslation(pts[0], pts[1], pts[2]);
      this.spheres.setMatrixAt(shown++, this.matrix);
      for (let q = 0; q + 1 < n; q++) {
        this.tailPos.set(pts.subarray(3 * q, 3 * q + 6), 6 * seg);
        // Age of each end in steps: the head is now; point q >= 1 is the
        // step ceil(t) - q, except the last, which may be cut at t - L.
        for (let e = 0; e < 2; e++) {
          const p = q + e;
          const age = p === 0 ? 0 : Math.min(L, t - (Math.ceil(t) - p));
          const f = L > 0 ? 0.85 * (age / L) : 0;
          const o = 6 * seg + 3 * e;
          this.tailCol[o] = c0.r + (1 - c0.r) * f;
          this.tailCol[o + 1] = c0.g + (1 - c0.g) * f;
          this.tailCol[o + 2] = c0.b + (1 - c0.b) * f;
        }
        seg++;
      }
    }
    this.spheres.count = shown;
    this.spheres.instanceMatrix.needsUpdate = true;
    const g = this.tails.geometry;
    (g.getAttribute("instanceStart") as unknown as { data: { needsUpdate: boolean } }).data.needsUpdate = true;
    (g.getAttribute("instanceColorStart") as unknown as { data: { needsUpdate: boolean } }).data.needsUpdate = true;
    g.instanceCount = seg;
  }

  setResolution(w: number, h: number): void {
    (this.tails.material as LineMaterial).resolution.set(w, h);
  }

  setVisible(on: boolean): void {
    this.spheres.visible = on;
    this.tails.visible = on;
  }
}
