// This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

// The domain's feature edges. The domain follows the slab and is remeshed
// every step, so its outline differs in shape and segment count from step to
// step and cannot be blended: the exporter keeps every step instead, and this
// shows the nearest one.

import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";
import type { Dataset } from "./data";

export class Outline {
  readonly lines: LineSegments2;
  private readonly buf: Float32Array;
  private shown = -1;

  constructor(private readonly data: Dataset) {
    let most = 0;
    for (let s = 0; s < data.manifest.steps; s++) most = Math.max(most, data.outlineStart[s + 1] - data.outlineStart[s]);
    this.buf = new Float32Array(6 * most);
    const g = new LineSegmentsGeometry();
    g.setPositions(this.buf);
    g.instanceCount = 0;
    const m = new LineMaterial({ color: 0x2a2d33, linewidth: 1.6 });
    m.toneMapped = false;
    this.lines = new LineSegments2(g, m);
    this.lines.frustumCulled = false;
  }

  update(step: number): void {
    const s = Math.max(0, Math.min(this.data.manifest.steps - 1, Math.round(step)));
    if (s === this.shown) return;
    this.shown = s;
    const a = this.data.outlineStart[s];
    const b = this.data.outlineStart[s + 1];
    this.buf.set(this.data.outline.subarray(6 * a, 6 * b));
    const g = this.lines.geometry;
    (g.getAttribute("instanceStart") as unknown as { data: { needsUpdate: boolean } }).data.needsUpdate = true;
    g.instanceCount = b - a;
  }

  setResolution(w: number, h: number): void {
    (this.lines.material as LineMaterial).resolution.set(w, h);
  }
}
