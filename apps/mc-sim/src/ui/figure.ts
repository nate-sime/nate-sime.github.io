// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The figure: every graph a view draws, side by side, in one CSS grid.
 *
 * Each panel is its own canvas with its own `Plot`, so hover and redraw stay
 * per graph. A view asks for the panels it needs; canvases are made on first
 * use and kept, the spare ones hidden. On a narrow screen the grid falls back
 * to one column of fixed-height panels that scrolls (see index.html).
 */

import { Plot, type PlotSpec } from "./plot";

export interface Layout {
  /** Columns of the grid; rows follow. */
  readonly cols: number;
  /** The first panel spans the whole first row. */
  readonly wideFirst?: boolean;
}

export class Figure {
  private readonly plots: Plot[] = [];
  private layout = "";

  constructor(private readonly el: HTMLElement) {}

  /** `n` panels, laid out and sized to their cells. */
  panels(n: number, layout: Layout = { cols: Math.min(n, 2) }): Plot[] {
    while (this.plots.length < n) {
      const canvas = document.createElement("canvas");
      canvas.className = "panel";
      this.el.appendChild(canvas);
      this.plots.push(new Plot(canvas));
    }
    const { cols } = layout, wide = !!layout.wideFirst && n > 1;
    const key = `${n}/${cols}/${wide}`;
    if (key !== this.layout) {
      this.layout = key;
      const rows = Math.ceil((n + (wide ? cols - 1 : 0)) / cols);
      this.el.style.setProperty("--cols", String(cols));
      this.el.style.setProperty("--rows", String(rows));
      this.plots.forEach((p, i) => {
        p.canvas.style.display = i < n ? "block" : "none";
        p.canvas.style.gridColumn = wide && i === 0 ? "1 / -1" : "";
      });
    }
    const used = this.plots.slice(0, n);
    for (const p of used) p.resize();
    return used;
  }

  /**
   * The view's panels, in its own grid, drawn empty — names and axes, no data
   * — while it waits for samples or cannot run: the page keeps its shape, and
   * the graphs fill in where they will be.
   */
  blank(specs: readonly Omit<PlotSpec, "series">[], layout?: Layout): void {
    this.panels(specs.length, layout).forEach((p, i) => p.draw({ ...specs[i], series: [] }));
  }
}
