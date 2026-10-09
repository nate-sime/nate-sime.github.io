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
 *
 * A panel with a `heading` is drawn as an instrument panel instead: the
 * heading bold at the top left, a `note` (a fitted rate) at the top right,
 * the legend as a strip of keys beneath them (`legend: "top"`), and no rotated
 * y label — the heading names the quantity. Such a panel can also draw grouped
 * bars, square and hollow markers (a projection rather than a measurement),
 * and a log y axis ticked in powers of two.
 */

export const SLOT = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9"] as const;

const INK = "rgba(207, 238, 255, 0.92)";
const MUTED = "rgba(207, 238, 255, 0.70)";
const GRID = "rgba(207, 238, 255, 0.10)";
const AXIS = "rgba(207, 238, 255, 0.28)";
export const SURFACE = "#05050c";
const FONT = "12px ui-monospace, monospace";
const SMALL = "11px ui-monospace, monospace";
const HEAD = "600 12px ui-monospace, monospace";

export interface Series {
  readonly label: string;
  readonly x: ArrayLike<number>;
  readonly y: ArrayLike<number>;
  readonly color: string;
  readonly width?: number;
  readonly dash?: readonly number[];
  readonly markers?: boolean;
  /** Marker shape; round by default. */
  readonly marker?: "circle" | "square";
  /** Markers (or bars) drawn as outlines — all of them, or point by point. */
  readonly hollow?: boolean | readonly boolean[];
  /** Bars up from the axis, grouped side by side with the other bar series. */
  readonly bars?: boolean;
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
  /** Centred above the plot; with a `heading`, a muted message across it instead. */
  readonly title?: string;
  /** Bold at the top left: the panel's name, and the y quantity's. */
  readonly heading?: string;
  /** Muted at the top right, beside the heading: a fitted rate, a summary. */
  readonly note?: string;
  readonly xlabel: string;
  readonly ylabel: string;
  readonly xlog?: boolean;
  readonly ylog?: boolean;
  /** A log y axis ticked and labelled in powers of 2 rather than of 10. */
  readonly ybase?: 2 | 10;
  readonly xlim?: readonly [number, number];
  /** Where the x ticks go, when the default spacing would skip some (levels). */
  readonly xticks?: readonly number[];
  readonly ylim?: readonly [number, number];
  /** y grows downward — a deflection read the way a beam sags. */
  readonly yflip?: boolean;
  readonly series: readonly Series[];
  /** Tick labels; default numeric. */
  readonly xfmt?: (v: number) => string;
  readonly yfmt?: (v: number) => string;
  /** Legend corner, top-right by default; "top" is a strip of keys under the heading. */
  readonly legend?: "tr" | "br" | "tl" | "bl" | "top";
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

/** Powers of 2 across [lo, hi], every one or every few, at most about six. */
function pow2Ticks(lo: number, hi: number): number[] {
  const a = Math.floor(Math.log2(lo) + 1e-9), b = Math.ceil(Math.log2(hi) - 1e-9);
  const every = Math.max(1, Math.ceil((b - a) / 6));
  const out: number[] = [];
  for (let k = Math.ceil(a / every) * every; k <= b; k += every) out.push(2 ** k);
  return out;
}

export const octave = (v: number) => `2${sup(Math.round(Math.log2(v)))}`;

const isHollow = (s: Series, i: number) => (typeof s.hollow === "boolean" ? s.hollow : !!s.hollow?.[i]);

/** A marker with a ring of the surface colour, so overlapping ones stay legible. */
function marker(ctx: CanvasRenderingContext2D, x: number, y: number, s: Series, hollow: boolean): void {
  const r = 4;
  ctx.beginPath();
  if (s.marker === "square") ctx.rect(x - r, y - r, 2 * r, 2 * r);
  else ctx.arc(x, y, r, 0, 2 * Math.PI);
  if (hollow) {
    ctx.lineWidth = 2;
    ctx.strokeStyle = SURFACE;
    ctx.stroke();
    ctx.fillStyle = SURFACE;
    ctx.fill();
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = s.color;
    ctx.stroke();
  } else {
    ctx.fillStyle = s.color;
    ctx.fill();
    ctx.lineWidth = 2;
    ctx.strokeStyle = SURFACE;
    ctx.stroke();
  }
}

/** A bar's outline: square at the foot, rounded at the top. */
function barPath(ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number): void {
  const r = Math.max(0, Math.min(4, w / 2, h));
  ctx.beginPath();
  ctx.moveTo(x, y + h);
  ctx.lineTo(x, y + r);
  ctx.arcTo(x, y, x + r, y, r);
  ctx.lineTo(x + w - r, y);
  ctx.arcTo(x + w, y, x + w, y + r, r);
  ctx.lineTo(x + w, y + h);
  ctx.closePath();
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

/** A point or bar drawn, where the hover can find it; a bar also by its body. */
interface Hit {
  readonly s: Series;
  readonly i: number;
  readonly X: number;
  readonly Y: number;
  readonly bar?: { readonly half: number; readonly base: number };
}

/** The legend strip under a heading: each key's place, wrapped to the panel's width. */
interface Strip {
  readonly keys: readonly { readonly s: Series; readonly x: number; readonly y: number }[];
  readonly height: number;
}

const SWATCH = 16;
const HEAD_H = 24;
const ROW_H = 15;

export class Plot {
  private readonly ctx: CanvasRenderingContext2D;
  private spec: PlotSpec | null = null;
  private mouse: { x: number; y: number } | null = null;
  private hits: Hit[] = [];
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
    this.hits = [];

    const pow2 = !!spec.ylog && spec.ybase === 2;
    const [x0, x1] = spec.xlim ?? extent(spec.series, (s) => s.x, !!spec.xlog);
    let [y0, y1] = spec.ylim ?? extent(spec.series, (s) => s.y, !!spec.ylog);
    // Whole octaves at the ends, so even a narrow range has two labelled ticks.
    if (pow2 && !spec.ylim) [y0, y1] = [2 ** Math.floor(Math.log2(y0)), 2 ** Math.ceil(Math.log2(y1))];
    const xt = spec.xticks ?? (spec.xlog ? logTicks(x0, x1) : niceTicks(x0, x1));
    const yt = pow2 ? pow2Ticks(y0, y1) : spec.ylog ? logTicks(y0, y1) : niceTicks(y0, y1);
    const xf = spec.xfmt ?? (spec.xlog ? decade : tickText);
    const yf = spec.yfmt ?? (pow2 ? octave : spec.ylog ? decade : tickText);

    // A heading, and its legend strip, take the top; without a y label the
    // left margin is just the tick labels' width.
    ctx.font = SMALL;
    const strip = spec.heading && spec.legend === "top" ? this.strip(spec, w) : null;
    const top = spec.heading ? HEAD_H + (strip?.height ?? 0) : spec.title ? 30 : 14;
    const left = spec.ylabel ? 74 : 14 + Math.max(16, ...yt.map((v) => ctx.measureText(yf(v)).width));
    const box = { l: left, t: top, r: w - 16, b: h - 46 };
    if (box.r - box.l < 40 || box.b - box.t < 40) return;
    const tx = spec.xlog ? Math.log10 : (v: number) => v, ty = spec.ylog ? Math.log10 : (v: number) => v;
    const [ax0, ax1, ay0, ay1] = [tx(x0), tx(x1), ty(y0), ty(y1)];
    const sx = (x: number) => box.l + ((tx(x) - ax0) / (ax1 - ax0)) * (box.r - box.l);
    const syRaw = (y: number) => box.b - ((ty(y) - ay0) / (ay1 - ay0)) * (box.b - box.t);
    const sy = spec.yflip ? (y: number) => box.t + box.b - syRaw(y) : syRaw;
    const axes: Axes = { sx, sy, box };

    // Grid and ticks.
    ctx.lineWidth = 1;
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
    if (spec.heading) this.header(spec, strip);
    else if (spec.title) {
      ctx.textBaseline = "top";
      const cx = (box.l + box.r) / 2;
      ctx.fillText(fit(ctx, spec.title, 2 * Math.min(cx - 8, w - 8 - cx)), cx, 8);
    }
    if (spec.ylabel) {
      ctx.save();
      ctx.translate(14, (box.t + box.b) / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textBaseline = "middle";
      ctx.fillText(spec.ylabel, 0, 0);
      ctx.restore();
    }

    ctx.save();
    ctx.beginPath();
    ctx.rect(box.l, box.t - 6, box.r - box.l + 6, box.b - box.t + 12);
    ctx.clip();
    spec.under?.(ctx, axes);
    const group = spec.series.filter((s) => s.bars);
    for (const s of spec.series) {
      if (s.bars) this.bars(s, group, axes, spec, x1 - x0);
      else this.line(s, axes, spec);
    }
    ctx.restore();
    if (spec.heading && spec.title) {
      // A heading's panel says what it is waiting for across the plot itself.
      ctx.font = SMALL;
      ctx.fillStyle = MUTED;
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText(fit(ctx, spec.title, box.r - box.l - 16), (box.l + box.r) / 2, (box.t + box.b) / 2);
    }
    spec.over?.(ctx, axes);
    if (!strip) this.legend(spec, box);
    this.hoverLayer(spec);
  }

  private visible(spec: PlotSpec, x: number, y: number): boolean {
    return Number.isFinite(x) && Number.isFinite(y) && (!spec.xlog || x > 0) && (!spec.ylog || y > 0);
  }

  private line(s: Series, a: Axes, spec: PlotSpec): void {
    const { ctx } = this;
    ctx.globalAlpha = s.alpha ?? 1;
    if (s.width !== 0) {
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
    }
    for (let i = 0; i < s.x.length; i++) {
      if (!this.visible(spec, s.x[i], s.y[i])) continue;
      const X = a.sx(s.x[i]), Y = a.sy(s.y[i]);
      if (s.markers) marker(ctx, X, Y, s, isHollow(s, i));
      if (!s.inert) this.hits.push({ s, i, X, Y });
    }
    ctx.globalAlpha = 1;
  }

  /** One bar series of `group`, each bar beside the others' at the same x. */
  private bars(s: Series, group: readonly Series[], a: Axes, spec: PlotSpec, span: number): void {
    const { ctx } = this, k = group.indexOf(s), n = group.length;
    const slot = Math.min(22, (((a.box.r - a.box.l) / Math.max(1, span)) * 0.7) / n);
    const off = (k - (n - 1) / 2) * (slot + 2);
    const base = spec.ylog ? a.box.b : Math.max(a.box.t, Math.min(a.box.b, a.sy(0)));
    ctx.globalAlpha = s.alpha ?? 1;
    for (let i = 0; i < s.x.length; i++) {
      if (!this.visible(spec, s.x[i], s.y[i])) continue;
      const X = a.sx(s.x[i]) + off, top = Math.min(a.sy(s.y[i]), base - 1);
      barPath(ctx, X - slot / 2, top, slot, base - top);
      if (isHollow(s, i)) {
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = s.color;
        ctx.stroke();
      } else {
        ctx.fillStyle = s.color;
        ctx.fill();
      }
      if (!s.inert) this.hits.push({ s, i, X, Y: top, bar: { half: slot / 2, base } });
    }
    ctx.globalAlpha = 1;
  }

  private header(spec: PlotSpec, strip: Strip | null): void {
    const { ctx, w } = this;
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.font = HEAD;
    ctx.fillStyle = INK;
    const head = fit(ctx, spec.heading!, w - 16);
    ctx.fillText(head, 8, 7);
    const room = w - 16 - ctx.measureText(head).width - 16;
    if (spec.note && room > 30) {
      ctx.font = SMALL;
      ctx.fillStyle = MUTED;
      ctx.textAlign = "right";
      ctx.fillText(fit(ctx, spec.note, room), w - 8, 8);
    }
    for (const k of strip?.keys ?? []) this.key(k.s, k.x, k.y);
  }

  /** The legend's keys in a row under the heading, wrapped where the panel is narrow. */
  private strip(spec: PlotSpec, w: number): Strip {
    const { ctx } = this;
    const keys: { s: Series; x: number; y: number }[] = [];
    let x = 8, row = 0;
    for (const s of spec.series) {
      if (s.unlisted) continue;
      const width = SWATCH + 6 + ctx.measureText(s.label).width;
      if (x > 8 && x + width > w - 8) { x = 8; row++; }
      keys.push({ s, x, y: HEAD_H + row * ROW_H + ROW_H / 2 - 2 });
      x += width + 14;
    }
    return { keys, height: keys.length ? (row + 1) * ROW_H + 2 : 0 };
  }

  /** A legend key: the series' mark as the plot draws it, then its label. */
  private key(s: Series, x: number, y: number): void {
    const { ctx } = this, cx = x + SWATCH / 2;
    ctx.globalAlpha = s.alpha ?? 1;
    if (s.bars) {
      barPath(ctx, cx - 4, y - 4, 8, 8);
      if (s.hollow === true) { ctx.lineWidth = 1.5; ctx.strokeStyle = s.color; ctx.stroke(); }
      else { ctx.fillStyle = s.color; ctx.fill(); }
    } else if (s.markers) {
      marker(ctx, cx, y, s, s.hollow === true);
    } else {
      ctx.strokeStyle = s.color;
      ctx.lineWidth = Math.min(s.width ?? 2, 2);
      ctx.setLineDash(s.dash ? [...s.dash] : []);
      ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + SWATCH, y); ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.globalAlpha = 1;
    ctx.font = SMALL;
    ctx.fillStyle = MUTED;
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.fillText(s.label, x + SWATCH + 6, y);
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

  private hoverLayer(spec: PlotSpec): void {
    const m = this.mouse;
    if (!m) return;
    let best: { hit: Hit; d: number } | null = null;
    for (const hit of this.hits) {
      const inBar = hit.bar && Math.abs(m.x - hit.X) <= hit.bar.half && m.y >= hit.Y - 6 && m.y <= hit.bar.base;
      const d = inBar ? 0 : Math.hypot(hit.X - m.x, hit.Y - m.y);
      if (d < 14 && (!best || d < best.d)) best = { hit, d };
    }
    if (!best) return;
    const { ctx } = this, { s, i, X, Y } = best.hit;
    if (spec.heading) {
      // A ring, so a hollow marker or a bar's top keeps its look under the pointer.
      ctx.beginPath();
      ctx.arc(X, Y, 7, 0, 2 * Math.PI);
      ctx.lineWidth = 1;
      ctx.strokeStyle = INK;
      ctx.stroke();
    } else {
      ctx.beginPath();
      ctx.arc(X, Y, 6, 0, 2 * Math.PI);
      ctx.lineWidth = 2;
      ctx.strokeStyle = SURFACE;
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(X, Y, 4.5, 0, 2 * Math.PI);
      ctx.fillStyle = s.color;
      ctx.fill();
    }
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
