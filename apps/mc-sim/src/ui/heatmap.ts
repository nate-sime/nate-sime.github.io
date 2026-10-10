// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * A scalar field over a rectangle, drawn inside a `Plot`: colour for value,
 * contour lines, a colour bar, and limits that keep the rectangle's shape.
 *
 * Two scales, by the job the colour does, both stepped off the reference
 * palette's ramps and validated against this page's #05050c surface
 * (`validate_palette.js --mode dark --ordinal`, each arm on its own: monotone
 * lightness, one hue, every step clear of the surface):
 *
 * - sequential, for a magnitude (a deflection, a variance): blue, from dark
 *   (near zero, receding toward the surface) to light;
 * - diverging, for a signed field (a mode shape, a fluctuation): red below
 *   zero, blue above, the arms at matched lightness, and a neutral gray at zero
 *   — "nothing" carries no hue.
 *
 * Contours are drawn in ink, not in a scale colour: the nodal lines of a mode
 * (its zero contour) are where Chladni's sand collects.
 */

import type { Axes, Plot } from "./plot";

type RGB = readonly [number, number, number];
const hex = (h: string): RGB => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16)) as unknown as RGB;

/** Stops, low to high, interpolated linearly between neighbours (each pair is one ramp step apart). */
const SEQUENTIAL = ["#104281", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef", "#cde2fb"].map(hex);
const NEG = ["#383835", "#9e3a39", "#e66767", "#f5b3b1"].map(hex);
const POS = ["#383835", "#1c5cab", "#3987e5", "#9ec5f4"].map(hex);

function ramp(stops: readonly RGB[], t: number): RGB {
  const u = Math.min(1, Math.max(0, t)) * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(u)), f = u - i;
  const a = stops[i], b = stops[i + 1];
  return [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]), a[2] + f * (b[2] - a[2])];
}

export type Scale = "sequential" | "diverging";

/** Colour of v on a scale spanning [0, top] (sequential) or [−top, top] (diverging). */
export function colourOf(scale: Scale, v: number, top: number): RGB {
  if (scale === "sequential") return ramp(SEQUENTIAL, v / top);
  return v < 0 ? ramp(NEG, -v / top) : ramp(POS, v / top);
}

/** A field sampled at the nodes of an nx × ny grid over [x0, x1] × [y0, y1], flat [iy·nx + ix]. */
export interface Grid {
  readonly nx: number;
  readonly ny: number;
  readonly x0: number;
  readonly x1: number;
  readonly y0: number;
  readonly y1: number;
  readonly v: ArrayLike<number>;
}

export function gridOf(nx: number, ny: number, x1: number, y1: number, v: ArrayLike<number>): Grid {
  return { nx, ny, x0: 0, x1, y0: 0, y1, v };
}

export const linspace = (a: number, b: number, n: number): Float64Array =>
  Float64Array.from({ length: n }, (_, i) => a + ((b - a) * i) / (n - 1));

let scratch: HTMLCanvasElement | null = null;

/** Paint the grid, one pixel per node, smoothed up to the axes. `factor` multiplies every value (an animated mode). */
export function paint(ctx: CanvasRenderingContext2D, a: Axes, g: Grid, scale: Scale, top: number, factor = 1): void {
  scratch ??= document.createElement("canvas");
  scratch.width = g.nx;
  scratch.height = g.ny;
  const sc = scratch.getContext("2d")!, img = sc.createImageData(g.nx, g.ny);
  for (let iy = 0; iy < g.ny; iy++)
    for (let ix = 0; ix < g.nx; ix++) {
      const [r, gr, b] = colourOf(scale, factor * g.v[iy * g.nx + ix], top), o = ((g.ny - 1 - iy) * g.nx + ix) * 4;
      img.data[o] = r; img.data[o + 1] = gr; img.data[o + 2] = b; img.data[o + 3] = 255;
    }
  sc.putImageData(img, 0, 0);
  // Pixel centres sit on the nodes: the image reaches half a cell past the rectangle, then is clipped to it.
  const dx = (g.x1 - g.x0) / (g.nx - 1), dy = (g.y1 - g.y0) / (g.ny - 1);
  const X0 = a.sx(g.x0 - dx / 2), X1 = a.sx(g.x1 + dx / 2), Y0 = a.sy(g.y1 + dy / 2), Y1 = a.sy(g.y0 - dy / 2);
  ctx.save();
  ctx.beginPath();
  ctx.rect(a.sx(g.x0), a.sy(g.y1), a.sx(g.x1) - a.sx(g.x0), a.sy(g.y0) - a.sy(g.y1));
  ctx.clip();
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(scratch, X0, Y0, X1 - X0, Y1 - Y0);
  ctx.restore();
}

/** The `level` contour by marching squares, as line segments in screen space. */
export function contour(ctx: CanvasRenderingContext2D, a: Axes, g: Grid, level: number, color: string, width = 1.5, dash: number[] = []): void {
  const { nx, ny, v } = g, dx = (g.x1 - g.x0) / (nx - 1), dy = (g.y1 - g.y0) / (ny - 1);
  const at = (ix: number, iy: number) => v[iy * nx + ix] - level;
  const cross = (x0: number, y0: number, f0: number, x1: number, y1: number, f1: number): [number, number] => {
    const t = f0 / (f0 - f1);
    return [a.sx(g.x0 + (x0 + t * (x1 - x0)) * dx), a.sy(g.y0 + (y0 + t * (y1 - y0)) * dy)];
  };
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.setLineDash(dash);
  ctx.beginPath();
  for (let iy = 0; iy < ny - 1; iy++)
    for (let ix = 0; ix < nx - 1; ix++) {
      const f = [at(ix, iy), at(ix + 1, iy), at(ix + 1, iy + 1), at(ix, iy + 1)];
      const c = [[ix, iy], [ix + 1, iy], [ix + 1, iy + 1], [ix, iy + 1]];
      const pts: [number, number][] = [];
      for (let e = 0; e < 4; e++) {
        const f0 = f[e], f1 = f[(e + 1) % 4];
        if ((f0 < 0) !== (f1 < 0)) pts.push(cross(c[e][0], c[e][1], f0, c[(e + 1) % 4][0], c[(e + 1) % 4][1], f1));
      }
      // Two crossings: one segment. Four (a saddle; crossings bottom, right, top,
      // left): the cell's mean decides which diagonal pair of corners connects,
      // and the segments cut off the other two.
      if (pts.length === 2) { ctx.moveTo(...pts[0]); ctx.lineTo(...pts[1]); }
      else if (pts.length === 4) {
        const centre = (f[0] + f[1] + f[2] + f[3]) / 4, isolate0 = (centre < 0) !== (f[0] < 0);
        const [p, q, r, s] = isolate0 ? [pts[0], pts[3], pts[1], pts[2]] : [pts[0], pts[1], pts[2], pts[3]];
        ctx.moveTo(...p); ctx.lineTo(...q); ctx.moveTo(...r); ctx.lineTo(...s);
      }
    }
  ctx.stroke();
  ctx.setLineDash([]);
}

/**
 * Plot limits that show [0, w] × [0, h] undistorted in this plot's box, with a
 * strip on the right for the colour bar. Mirrors `Plot`'s box margins.
 */
export function shapeLimits(plot: Plot, w: number, h: number, titled = true): { xlim: [number, number]; ylim: [number, number] } {
  const bw = Math.max(40, plot.canvas.clientWidth - 74 - 16), bh = Math.max(40, plot.canvas.clientHeight - (titled ? 30 : 14) - 46);
  const pad = 0.04, bar = 0.24;
  // Data units per pixel: whichever direction is tighter sets the scale.
  const need = { x: w * (1 + 2 * pad + bar), y: h * (1 + 2 * pad) };
  const k = Math.max(need.x / bw, need.y / bh);
  const xs = bw * k, ys = bh * k;
  const x0 = -w * pad - (xs - need.x) / 2, y0 = -h * pad - (ys - need.y) / 2;
  return { xlim: [x0, x0 + xs], ylim: [y0, y0 + ys] };
}

/** A vertical colour bar at the right of the box, with its end values in ink. */
export function colourBar(ctx: CanvasRenderingContext2D, a: Axes, scale: Scale, top: number, label: string, fmt: (v: number) => string): void {
  const { box } = a, w = 10, x = box.r - 64, y0 = box.t + 18, y1 = box.b - 18, n = 64;
  for (let i = 0; i < n; i++) {
    const t = 1 - (i + 0.5) / n, v = scale === "sequential" ? t * top : (2 * t - 1) * top;
    const [r, g, b] = colourOf(scale, v, top);
    ctx.fillStyle = `rgb(${r | 0}, ${g | 0}, ${b | 0})`;
    ctx.fillRect(x, y0 + ((y1 - y0) * i) / n, w, (y1 - y0) / n + 0.5);
  }
  ctx.fillStyle = "rgba(207, 238, 255, 0.92)";
  ctx.font = "11px ui-monospace, monospace";
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  ctx.fillText(fmt(top), x + w + 5, y0);
  ctx.fillText(fmt(scale === "sequential" ? 0 : -top), x + w + 5, y1);
  if (scale === "diverging") ctx.fillText("0", x + w + 5, (y0 + y1) / 2);
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.fillText(label, x + w / 2, y0 - 4);
}

/** The rectangle's outline, with edge conditions marked: hatched for clamped, dashed for simply supported, plain for free. */
export function plateOutline(ctx: CanvasRenderingContext2D, a: Axes, w: number, h: number, edges: string): void {
  const ink = "rgba(207, 238, 255, 0.85)";
  const corners: [number, number][] = [[0, 0], [0, h], [w, h], [w, 0]];
  // Leissa's order: x = 0, y = 0, x = w, y = h — as segments with their outward normal.
  const sides: { p: [number, number]; q: [number, number]; n: [number, number] }[] = [
    { p: corners[0], q: corners[1], n: [-1, 0] },
    { p: corners[0], q: corners[3], n: [0, 1] },
    { p: corners[3], q: corners[2], n: [1, 0] },
    { p: corners[1], q: corners[2], n: [0, -1] },
  ];
  sides.forEach(({ p, q, n }, i) => {
    const e = edges[i], P = [a.sx(p[0]), a.sy(p[1])], Q = [a.sx(q[0]), a.sy(q[1])];
    ctx.strokeStyle = ink;
    ctx.lineWidth = e === "F" ? 1 : 2;
    ctx.setLineDash(e === "S" ? [6, 4] : []);
    ctx.beginPath(); ctx.moveTo(P[0], P[1]); ctx.lineTo(Q[0], Q[1]); ctx.stroke();
    ctx.setLineDash([]);
    if (e === "C") {
      // Screen normal: y flips.
      const nx = n[0], ny = -n[1], len = Math.hypot(Q[0] - P[0], Q[1] - P[1]), tx = (Q[0] - P[0]) / len, ty = (Q[1] - P[1]) / len;
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let s = 4; s < len; s += 8) {
        const X = P[0] + tx * s, Y = P[1] + ty * s;
        ctx.moveTo(X, Y);
        ctx.lineTo(X + 6 * nx - 4 * tx, Y + 6 * ny - 4 * ty);
      }
      ctx.stroke();
    }
  });
}

/** Element boundaries over a rectangle: the lines x = xb[i] and y = yb[j], interior ones only. */
export function meshLines(
  ctx: CanvasRenderingContext2D, a: Axes, xb: ArrayLike<number>, yb: ArrayLike<number>, color: string, width = 1,
): void {
  const x0 = xb[0], x1 = xb[xb.length - 1], y0 = yb[0], y1 = yb[yb.length - 1];
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.beginPath();
  for (let i = 1; i < xb.length - 1; i++) { ctx.moveTo(a.sx(xb[i]), a.sy(y0)); ctx.lineTo(a.sx(xb[i]), a.sy(y1)); }
  for (let j = 1; j < yb.length - 1; j++) { ctx.moveTo(a.sx(x0), a.sy(yb[j])); ctx.lineTo(a.sx(x1), a.sy(yb[j])); }
  ctx.stroke();
}

/** Element boundaries along a line, as ticks across y = y0 in screen space `half` pixels each way. */
export function meshTicks(ctx: CanvasRenderingContext2D, a: Axes, breaks: ArrayLike<number>, y0: number, color: string, half = 5): void {
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  for (let i = 0; i < breaks.length; i++) {
    const X = a.sx(breaks[i]), Y = a.sy(y0);
    ctx.moveTo(X, Y - half);
    ctx.lineTo(X, Y + half);
  }
  ctx.stroke();
}
