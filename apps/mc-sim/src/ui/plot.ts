// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * One 2-D chart on one canvas: linear or log axes, line series with optional
 * markers and dashes, a legend, and a nearest-point hover readout.
 *
 * Colours are roles, not choices. The series slots are the reference
 * categorical palette's dark steps, validated against this page's #05050c
 * surface (`validate_palette.js --mode dark`: CVD ΔE ≥ 8.4 adjacent, contrast
 * ≥ 3:1), and a series keeps its slot whatever else is on screen — p = 3 is
 * orange in every view. Text is in ink, never in a series colour.
 */

export const SLOT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9"] as const;

const INK = "rgba(207, 238, 255, 0.92)";
const MUTED = "rgba(207, 238, 255, 0.70)";
const GRID = "rgba(207, 238, 255, 0.10)";
const AXIS = "rgba(207, 238, 255, 0.28)";
export const SURFACE = "#05050c";
const FONT = "12px ui-monospace, monospace";
const SMALL = "11px ui-monospace, monospace";

export interface Series {
  readonly label: string;
  readonly x: ArrayLike<number>;
  readonly y: ArrayLike<number>;
  readonly color: string;
  readonly width?: number;
  readonly dash?: readonly number[];
  readonly markers?: boolean;
  readonly alpha?: number;
  /** Left out of the legend (still hoverable). */
  readonly unlisted?: boolean;
  /** Left out of hover (decoration, e.g. a reference line). */
  readonly inert?: boolean;
}

export interface Axes {
  readonly sx: (x: number) => number;
  readonly sy: (y: number) => number;
  readonly box: { l: number; t: number; r: number; b: number };
}

export interface PlotSpec {
  readonly title?: string;
  readonly xlabel: string;
  readonly ylabel: string;
  readonly xlog?: boolean;
  readonly ylog?: boolean;
  readonly xlim?: readonly [number, number];
  readonly ylim?: readonly [number, number];
  /** y grows downward — a deflection read the way a beam sags. */
  readonly yflip?: boolean;
  readonly series: readonly Series[];
  /** Tick labels; default numeric. */
  readonly xfmt?: (v: number) => string;
  readonly yfmt?: (v: number) => string;
  /** Legend corner; top-right by default. */
  readonly legend?: "tr" | "br" | "tl" | "bl";
  /** Hover text for a point. */
  readonly hover?: (s: Series, i: number) => string;
  /** Drawn after the grid, before the series, in data coordinates. */
  readonly under?: (ctx: CanvasRenderingContext2D, a: Axes) => void;
  /** Drawn over everything. */
  readonly over?: (ctx: CanvasRenderingContext2D, a: Axes) => void;
}

const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sup = (n: number) => (n < 0 ? "⁻" : "") + String(Math.abs(n)).split("").map((d) => SUP[+d]).join("");
export const decade = (v: number) => `10${sup(Math.round(Math.log10(v)))}`;

function niceTicks(lo: number, hi: number, count = 5): number[] {
  const raw = (hi - lo) / count, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => raw <= s) ?? 10 * mag;
  const out: number[] = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}

function logTicks(lo: number, hi: number): number[] {
  const a = Math.floor(Math.log10(lo) + 1e-9), b = Math.ceil(Math.log10(hi) - 1e-9);
  const every = Math.max(1, Math.ceil((b - a) / 7));
  const out: number[] = [];
  for (let k = a; k <= b; k += every) out.push(10 ** k);
  return out;
}

/** `text`, cut to `width` pixels with an ellipsis — a title in a small panel. */
function fit(ctx: CanvasRenderingContext2D, text: string, width: number): string {
  if (ctx.measureText(text).width <= width) return text;
  let lo = 0, hi = text.length;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (ctx.measureText(`${text.slice(0, mid)}…`).width <= width) lo = mid; else hi = mid - 1;
  }
  return `${text.slice(0, lo).trimEnd()}…`;
}

const tickText = (v: number) => {
  if (v === 0) return "0";
  const a = Math.abs(v);
  return a >= 1e4 || a < 1e-3 ? v.toExponential(0) : String(+v.toPrecision(4));
};

function extent(series: readonly Series[], pick: (s: Series) => ArrayLike<number>, log: boolean): [number, number] {
  let lo = Infinity, hi = -Infinity;
  for (const s of series)
    for (const v of Array.from(pick(s)))
      if (Number.isFinite(v) && (!log || v > 0)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  if (!Number.isFinite(lo)) return log ? [1e-3, 1] : [0, 1];
  if (lo === hi) return log ? [lo / 10, hi * 10] : [lo - 1, hi + 1];
  if (log) return [lo / 1.6, hi * 1.6];
  const pad = 0.06 * (hi - lo);
  return [lo - pad, hi + pad];
}

export class Plot {
  private readonly ctx: CanvasRenderingContext2D;
  private spec: PlotSpec | null = null;
  private mouse: { x: number; y: number } | null = null;
  private dpr = 1;
  private w = 0;
  private h = 0;

  constructor(readonly canvas: HTMLCanvasElement) {
    this.ctx = canvas.getContext("2d")!;
    canvas.addEventListener("pointermove", (e) => {
      const r = canvas.getBoundingClientRect();
      this.mouse = { x: e.clientX - r.left, y: e.clientY - r.top };
      this.redraw();
    });
    canvas.addEventListener("pointerleave", () => { this.mouse = null; this.redraw(); });
  }

  /** Match the backing store to the element's CSS size; true if it changed. */
  resize(): boolean {
    const dpr = window.devicePixelRatio || 1;
    const w = this.canvas.clientWidth, h = this.canvas.clientHeight;
    if (w === this.w && h === this.h && dpr === this.dpr) return false;
    this.w = w; this.h = h; this.dpr = dpr;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    return true;
  }

  draw(spec: PlotSpec): void {
    this.spec = spec;
    this.redraw();
  }

  redraw(): void {
    const spec = this.spec;
    if (!spec) return;
    const { ctx, w, h } = this;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = SURFACE;
    ctx.fillRect(0, 0, w, h);

    const [x0, x1] = spec.xlim ?? extent(spec.series, (s) => s.x, !!spec.xlog);
    const [y0, y1] = spec.ylim ?? extent(spec.series, (s) => s.y, !!spec.ylog);
    const box = { l: 74, t: spec.title ? 30 : 14, r: w - 16, b: h - 46 };
    if (box.r - box.l < 40 || box.b - box.t < 40) return;
    const tx = spec.xlog ? Math.log10 : (v: number) => v, ty = spec.ylog ? Math.log10 : (v: number) => v;
    const [ax0, ax1, ay0, ay1] = [tx(x0), tx(x1), ty(y0), ty(y1)];
    const sx = (x: number) => box.l + ((tx(x) - ax0) / (ax1 - ax0)) * (box.r - box.l);
    const syRaw = (y: number) => box.b - ((ty(y) - ay0) / (ay1 - ay0)) * (box.b - box.t);
    const sy = spec.yflip ? (y: number) => box.t + box.b - syRaw(y) : syRaw;
    const axes: Axes = { sx, sy, box };

    // Grid and ticks.
    ctx.font = SMALL;
    ctx.lineWidth = 1;
    const xt = spec.xlog ? logTicks(x0, x1) : niceTicks(x0, x1);
    const yt = spec.ylog ? logTicks(y0, y1) : niceTicks(y0, y1);
    const xf = spec.xfmt ?? (spec.xlog ? decade : tickText), yf = spec.yfmt ?? (spec.ylog ? decade : tickText);
    ctx.strokeStyle = GRID;
    ctx.fillStyle = MUTED;
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const v of xt) {
      const X = Math.round(sx(v)) + 0.5;
      if (X < box.l - 1 || X > box.r + 1) continue;
      ctx.beginPath(); ctx.moveTo(X, box.t); ctx.lineTo(X, box.b); ctx.stroke();
      ctx.fillText(xf(v), X, box.b + 6);
    }
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    for (const v of yt) {
      const Y = Math.round(sy(v)) + 0.5;
      if (Y < box.t - 1 || Y > box.b + 1) continue;
      ctx.beginPath(); ctx.moveTo(box.l, Y); ctx.lineTo(box.r, Y); ctx.stroke();
      ctx.fillText(yf(v), box.l - 8, Y);
    }
    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(box.l + 0.5, box.t); ctx.lineTo(box.l + 0.5, box.b + 0.5); ctx.lineTo(box.r, box.b + 0.5);
    ctx.stroke();

    // Labels.
    ctx.fillStyle = INK;
    ctx.font = FONT;
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    ctx.fillText(spec.xlabel, (box.l + box.r) / 2, h - 6);
    if (spec.title) {
      ctx.textBaseline = "top";
      const cx = (box.l + box.r) / 2;
      ctx.fillText(fit(ctx, spec.title, 2 * Math.min(cx - 8, w - 8 - cx)), cx, 8);
    }
    ctx.save();
    ctx.translate(14, (box.t + box.b) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = "middle";
    ctx.fillText(spec.ylabel, 0, 0);
    ctx.restore();

    ctx.save();
    ctx.beginPath();
    ctx.rect(box.l, box.t - 6, box.r - box.l + 6, box.b - box.t + 12);
    ctx.clip();
    spec.under?.(ctx, axes);
    for (const s of spec.series) this.line(s, axes, spec);
    ctx.restore();
    spec.over?.(ctx, axes);
    this.legend(spec, box);
    this.hoverLayer(spec, axes);
  }

  private visible(spec: PlotSpec, x: number, y: number): boolean {
    return Number.isFinite(x) && Number.isFinite(y) && (!spec.xlog || x > 0) && (!spec.ylog || y > 0);
  }

  private line(s: Series, a: Axes, spec: PlotSpec): void {
    const { ctx } = this;
    ctx.globalAlpha = s.alpha ?? 1;
    ctx.strokeStyle = s.color;
    ctx.lineWidth = s.width ?? 2;
    ctx.lineJoin = "round";
    ctx.setLineDash(s.dash ? [...s.dash] : []);
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < s.x.length; i++) {
      if (!this.visible(spec, s.x[i], s.y[i])) { pen = false; continue; }
      const X = a.sx(s.x[i]), Y = a.sy(s.y[i]);
      if (pen) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
      pen = true;
    }
    ctx.stroke();
    ctx.setLineDash([]);
    if (s.markers) {
      for (let i = 0; i < s.x.length; i++) {
        if (!this.visible(spec, s.x[i], s.y[i])) continue;
        ctx.beginPath();
        ctx.arc(a.sx(s.x[i]), a.sy(s.y[i]), 4, 0, 2 * Math.PI);
        ctx.fillStyle = s.color;
        ctx.fill();
        ctx.lineWidth = 2;
        ctx.strokeStyle = SURFACE;
        ctx.stroke();
      }
    }
    ctx.globalAlpha = 1;
  }

  private legend(spec: PlotSpec, box: Axes["box"]): void {
    const items = spec.series.filter((s) => !s.unlisted);
    if (items.length < 2) return;
    const { ctx } = this;
    ctx.font = SMALL;
    const rowH = 16, sw = 22;
    const width = Math.max(...items.map((s) => ctx.measureText(s.label).width)) + sw + 22;
    const h = items.length * rowH + 8, corner = spec.legend ?? "tr";
    const x = corner === "tl" || corner === "bl" ? box.l + 6 : box.r - width - 6;
    const y = corner === "br" || corner === "bl" ? box.b - h - 6 : box.t + 6;
    ctx.fillStyle = "rgba(5, 5, 12, 0.82)";
    ctx.fillRect(x, y, width, h);
    items.forEach((s, i) => {
      const Y = y + 4 + i * rowH + rowH / 2;
      ctx.strokeStyle = s.color;
      ctx.globalAlpha = s.alpha ?? 1;
      ctx.lineWidth = s.width ?? 2;
      ctx.setLineDash(s.dash ? [...s.dash] : []);
      ctx.beginPath(); ctx.moveTo(x + 8, Y); ctx.lineTo(x + 8 + sw, Y); ctx.stroke();
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      ctx.fillStyle = INK;
      ctx.textAlign = "left";
      ctx.textBaseline = "middle";
      ctx.fillText(s.label, x + sw + 14, Y);
    });
  }

  private hoverLayer(spec: PlotSpec, a: Axes): void {
    const m = this.mouse;
    if (!m) return;
    let best: { s: Series; i: number; d: number } | null = null;
    for (const s of spec.series) {
      if (s.inert) continue;
      for (let i = 0; i < s.x.length; i++) {
        if (!this.visible(spec, s.x[i], s.y[i])) continue;
        const d = Math.hypot(a.sx(s.x[i]) - m.x, a.sy(s.y[i]) - m.y);
        if (d < 14 && (!best || d < best.d)) best = { s, i, d };
      }
    }
    if (!best) return;
    const { ctx } = this, { s, i } = best;
    const X = a.sx(s.x[i]), Y = a.sy(s.y[i]);
    ctx.beginPath();
    ctx.arc(X, Y, 6, 0, 2 * Math.PI);
    ctx.lineWidth = 2;
    ctx.strokeStyle = SURFACE;
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(X, Y, 4.5, 0, 2 * Math.PI);
    ctx.fillStyle = s.color;
    ctx.fill();
    const lines = (spec.hover?.(s, i) ?? `${s.label}\nx = ${tickText(s.x[i])}\ny = ${s.y[i].toPrecision(5)}`).split("\n");
    ctx.font = SMALL;
    const tw = Math.max(...lines.map((l) => ctx.measureText(l).width)) + 16, th = lines.length * 15 + 8;
    let bx = X + 12, by = Y - th - 8;
    if (bx + tw > this.w - 4) bx = X - tw - 12;
    if (by < 4) by = Y + 12;
    ctx.fillStyle = "rgba(14, 16, 28, 0.95)";
    ctx.strokeStyle = AXIS;
    ctx.lineWidth = 1;
    ctx.fillRect(bx, by, tw, th);
    ctx.strokeRect(bx + 0.5, by + 0.5, tw - 1, th - 1);
    ctx.fillStyle = INK;
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    lines.forEach((l, k) => ctx.fillText(l, bx + 8, by + 5 + k * 15));
  }
}
