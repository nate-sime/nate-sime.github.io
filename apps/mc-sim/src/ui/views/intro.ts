// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * View 0, the introduction: the problem before the methods. Not a figure and a
 * readout like the other views but a page of its own, over the whole frame —
 * the pane is put away, and the page ends in the guided tours.
 *
 * Two parts come off a production line: a cantilever beam and a free square
 * plate (`manufacture.ts`). Each is drawn in 3-D, vibrating in one of its
 * lowest modes, inside a wireframe ghost of the deterministic solution — the
 * design, solved, in the same mode and struck at the same moment. The made
 * part's frequency is not the design's, so the two drift apart within a few
 * swings; and its shape is not the design's either, which the plate shows
 * best: for its second and third modes the made plate's nodal lines (white, on
 * its surface) part company with the ghost's (dashed). Colour is the part's
 * thickness against the drawing's.
 *
 * Below, every part made so far, as a dot at its ω₁ over the design's: a
 * distribution filling in, which is all a Monte Carlo estimate is — and a line
 * on what pinning its mean down would cost, at this mesh's measured solve time.
 */

import { MADE_MODES, PLATE_POINTS, beamX, design, make, plateX, type Kind, type Part } from "../../manufacture";
import { SLOT } from "../plot";
import { cellZeros, orbit, paintQuads, projector, type Camera, type Projector, type Quad, type Vec3 } from "../scene";
import { TOUR_NAMES, type TourName } from "../tours";
import type { ViewResult } from "./view";

export interface IntroActions {
  startTour(name: TourName): void;
  /** Leave the introduction for the views. */
  explore(): void;
  /** Ask for a frame. */
  request(): void;
}

/** One swing of the design, in any mode: slowed far below life, for the eye. */
const PERIOD_MS = 1500;
/** The parts are struck again every this many of the design's swings, and ring down in between. */
const RESTRIKE = 6;
/** Depth deviations are drawn this many times their size. */
const DEPTH_GAIN = 2;
/** "+ more": how many parts a click makes. */
const BATCH = 20;

const TOUR_BLURB: Record<TourName, string> = {
  "Monte Carlo and √N": "the estimate, its error bar, and why that bar shrinks so slowly",
  "discretisation hierarchies": "how fast the finite-element error falls as the mesh is refined",
  "the bias–variance split": "the point past which more samples on one mesh buy nothing",
  "multilevel Monte Carlo": "most samples on coarse meshes, a few on fine ones: the same accuracy for a fraction of the work",
};

const NAME: Record<Kind, string> = { beam: "a cantilever beam", plate: "a free square plate" };

type RGB = readonly [number, number, number];
const THIN: RGB = [217, 89, 38], STEEL: RGB = [150, 164, 180], THICK: RGB = [57, 135, 229];

/** Thinner than drawn toward orange, thicker toward blue, as drawn in steel. Full colour at ±25%. */
function depthColour(h: number): RGB {
  const t = Math.max(-1, Math.min(1, (h - 1) / 0.25)), to = t < 0 ? THIN : THICK, a = Math.abs(t);
  return [STEEL[0] + a * (to[0] - STEEL[0]), STEEL[1] + a * (to[1] - STEEL[1]), STEEL[2] + a * (to[2] - STEEL[2])];
}

interface Card {
  readonly kind: Kind;
  readonly canvas: HTMLCanvasElement;
  readonly stat: HTMLElement;
  readonly make: HTMLButtonElement;
  readonly more: HTMLButtonElement;
  readonly chips: HTMLButtonElement[];
  readonly cam: Camera;
  /** Every part made, in order; the last is on screen. */
  readonly made: Part[];
  mode: number;
  /** When the part on screen was struck. */
  t0: number;
}

const el = <K extends keyof HTMLElementTagNameMap>(tag: K, cls: string, parent?: HTMLElement, text?: string) => {
  const e = document.createElement(tag);
  e.className = cls;
  if (text !== undefined) e.textContent = text;
  parent?.appendChild(e);
  return e;
};

export class Intro {
  private readonly cards: Card[] = [];
  private readonly strip: HTMLCanvasElement;
  private readonly cost: HTMLElement;
  private busy = false;
  private now = 0;

  constructor(root: HTMLElement, private readonly actions: IntroActions) {
    const head = el("header", "in-head", root);
    el("div", "in-kicker", head, "MC-sim · introduction");
    el("h1", "in-title", head, "No two made the same");
    el("p", "in-lead", head,
      "Consider a plant that makes beams and plates. An engineer designs them with deterministic models: given a " +
      "geometry and a material, a mathematical model approximates how the part will behave. But no manufacturing " +
      "process is perfect. Each part's thickness, mass and stiffness come out a little different from the drawing, " +
      "and that can profoundly change how it behaves, for example how it vibrates. To design with this in mind, we " +
      "treat each part as random and estimate its statistics by simulation: its mean frequency, its spread, and the " +
      "chance it resonates with whatever shakes it.");

    const stage = el("div", "in-stage", root);
    for (const kind of ["beam", "plate"] as const) this.cards.push(this.card(stage, kind));

    const key = el("div", "in-key", root);
    const swatch = (rgb: RGB, label: string) => {
      const s = el("span", "in-swatch", key);
      el("i", "", s).style.background = `rgb(${rgb.join(",")})`;
      s.append(label);
    };
    swatch(THIN, "thinner than drawn");
    swatch(STEEL, "as drawn");
    swatch(THICK, "thicker");
    el("span", "in-swatch in-ghost", key, "ghost: the deterministic solution");
    const nodal = el("span", "in-swatch", key);
    el("i", "in-nodal", nodal);
    nodal.append("the made plate's nodal lines (the ghost's are dashed)");

    const batch = el("div", "in-batch", root);
    this.strip = el("canvas", "in-strip", batch);
    this.cost = el("p", "in-cost", batch);
    el("p", "in-note", root,
      "Exaggerated for the picture: the stiffness varies by about 50% (σ of log EI = 0.5), thickness deviations are " +
      "drawn twice their size, and the swings are huge and slow. Each part is solved live by the app's cubic-spline " +
      "finite elements and eigensolver, and is the same part every time you make it. Drag a part to turn it.");

    const next = el("footer", "in-next", root);
    el("h2", "in-ask", next, "So what can we do?");
    el("p", "in-lead", next,
      "Monte Carlo is the honest answer: make thousands of virtual parts, solve each one, average. But its error falls " +
      "only as 1/√N, so each extra digit costs a hundred times the solves, and every solve on a mesh fine enough to trust " +
      "is expensive. The guided tours build, one idea at a time, to a way out:");
    const tours = el("div", "in-tours", next);
    TOUR_NAMES.forEach((name, i) => {
      const b = el("button", "in-tour", tours);
      el("span", "in-tour-n", b, String(i + 1));
      const words = el("span", "in-tour-words", b);
      el("b", "", words, name);
      el("span", "", words, TOUR_BLURB[name]);
      b.addEventListener("click", () => actions.startTour(name));
    });
    el("button", "in-explore", next, "or explore the views on your own  →").addEventListener("click", actions.explore);

    // The first of each, made on arrival, so the page opens moving.
    for (const c of this.cards) this.manufacture(c);
  }

  private card(stage: HTMLElement, kind: Kind): Card {
    const box = el("section", "in-card", stage);
    const top = el("div", "in-card-head", box);
    el("span", "in-name", top, NAME[kind]);
    const chipRow = el("span", "in-chips", top, "mode ");
    const chips = Array.from({ length: MADE_MODES }, (_, n) => el("button", "in-chip", chipRow, String(n + 1)));
    const canvas = el("canvas", "in-view", box);
    const foot = el("div", "in-card-foot", box);
    const make = el("button", "in-make", foot, `manufacture a new ${kind}`);
    const more = el("button", "in-more", foot, `+ ${BATCH} more`);
    const stat = el("div", "in-stat", box);
    const cam: Camera = kind === "beam" ? { yaw: -0.42, pitch: 0.32 } : { yaw: 0.55, pitch: 0.62 };
    const c: Card = { kind, canvas, stat, make, more, chips, cam, made: [], mode: kind === "beam" ? 0 : 1, t0: 0 };
    chips.forEach((b, n) => b.addEventListener("click", () => { c.mode = n; c.t0 = this.now; this.actions.request(); }));
    make.addEventListener("click", () => this.manufacture(c));
    more.addEventListener("click", () => this.batch(c));
    orbit(canvas, cam, this.actions.request);
    return c;
  }

  private manufacture(c: Card): void {
    c.made.push(make(c.kind, c.made.length));
    c.t0 = this.now;
    this.actions.request();
  }

  /** Make a batch, one part a tick, so the page keeps moving while it fills. */
  private batch(c: Card): void {
    if (this.busy) return;
    this.busy = true;
    for (const k of this.cards) { k.make.disabled = true; k.more.disabled = true; }
    let left = BATCH;
    const step = () => {
      this.manufacture(c);
      if (--left > 0) { setTimeout(step, 0); return; }
      this.busy = false;
      for (const k of this.cards) { k.make.disabled = false; k.more.disabled = false; }
    };
    step();
  }

  render(t: number): ViewResult {
    this.now = t;
    for (const c of this.cards) this.drawCard(c, t);
    this.drawStrip();
    this.writeCost();
    return { readout: "", animate: true };
  }

  private drawCard(c: Card, t: number): void {
    const ctx = fit(c.canvas);
    if (!ctx) return;
    const part = c.made[c.made.length - 1], d = design(c.kind), n = c.mode;
    const ratio = part.omega[n] / d.omega[n];
    // Both struck at t0; the design swings once a period, the part `ratio` times as fast. Each rings down until struck again.
    const since = (Math.max(0, t - c.t0) / PERIOD_MS) % RESTRIKE;
    const ring = Math.exp(-0.28 * since);
    const zDesign = ring * Math.cos(2 * Math.PI * since), zPart = ring * Math.cos(2 * Math.PI * since * ratio);
    const dpr = window.devicePixelRatio || 1, W = c.canvas.width / dpr, H = c.canvas.height / dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    if (c.kind === "beam") drawBeam(ctx, c.cam, W, H, part, d, n, zPart, zDesign);
    else drawPlate(ctx, c.cam, W, H, part, d, n, zPart, zDesign);

    for (const [i, b] of c.chips.entries()) b.classList.toggle("on", i === n);
    const no = part.index + 1, dev = (ratio - 1) * 100;
    const html =
      `<b>${c.kind} #${no}</b> · mode ${n + 1} rings at <b>${ratio.toFixed(3)}×</b> the design's frequency ` +
      `<span class="${Math.abs(dev) < 0.05 ? "" : dev > 0 ? "up" : "down"}">(${dev >= 0 ? "+" : "−"}${Math.abs(dev).toFixed(1)}%)</span>`;
    setHtml(c.stat, html);
  }

  private drawStrip(): void {
    const ctx = fit(this.strip);
    if (!ctx) return;
    const dpr = window.devicePixelRatio || 1, W = this.strip.width / dpr, H = this.strip.height / dpr;
    ctx.save();
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, W, H);
    const lo = 0.7, hi = 1.3, left = 64, right = 14, X = (r: number) => left + ((r - lo) / (hi - lo)) * (W - left - right);
    const rows = H - 34, rowH = rows / 2, axisY = 10 + rows;
    ctx.font = "11px ui-monospace, monospace";
    ctx.textBaseline = "middle";
    // The design's line.
    ctx.strokeStyle = "rgba(207, 238, 255, 0.35)";
    ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(X(1), 4); ctx.lineTo(X(1), axisY); ctx.stroke();
    ctx.setLineDash([]);
    // Axis.
    ctx.strokeStyle = "rgba(207, 238, 255, 0.25)";
    ctx.beginPath(); ctx.moveTo(left, axisY); ctx.lineTo(W - right, axisY); ctx.stroke();
    ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
    ctx.textAlign = "center";
    for (let r = 0.7; r <= 1.3001; r += 0.1) {
      ctx.fillRect(X(r) - 0.5, axisY, 1, 4);
      ctx.fillText(r.toFixed(1), X(r), axisY + 12);
    }
    ctx.textAlign = "right";
    ctx.fillText("× design ω₁", W - right, axisY + 24 > H - 2 ? axisY - 8 : axisY + 24);
    this.cards.forEach((c, row) => {
      const y0 = 10 + row * rowH, mid = y0 + rowH / 2, colour = SLOT[row];
      ctx.textAlign = "left";
      ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
      ctx.fillText(`${c.kind}s`, 0, mid);
      ctx.fillText(`${c.made.length}`, 0, mid + 12);
      const r = c.made.map((p) => p.omega[0] / design(c.kind).omega[0]);
      // Mean and its 95% interval, once there are two to spread.
      if (r.length > 1) {
        const m = mean(r), half = (1.96 * sd(r)) / Math.sqrt(r.length);
        ctx.fillStyle = "rgba(207, 238, 255, 0.10)";
        ctx.fillRect(X(m - half), y0 + 2, Math.max(1, X(m + half) - X(m - half)), rowH - 4);
        ctx.fillStyle = "rgba(207, 238, 255, 0.9)";
        ctx.fillRect(X(m) - 0.75, y0 + 2, 1.5, rowH - 4);
      }
      ctx.globalAlpha = 0.75;
      ctx.fillStyle = colour;
      r.forEach((v, i) => {
        const y = mid + (hash(i + 7 * row) - 0.5) * (rowH - 10);
        ctx.beginPath(); ctx.arc(X(clamp(v, lo, hi)), y, 3.2, 0, 2 * Math.PI); ctx.fill();
      });
      ctx.globalAlpha = 1;
      if (r.length) {
        const i = r.length - 1, y = mid + (hash(i + 7 * row) - 0.5) * (rowH - 10);
        ctx.strokeStyle = "#fff"; ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.arc(X(clamp(r[i], lo, hi)), y, 5.5, 0, 2 * Math.PI); ctx.stroke();
      }
    });
    ctx.restore();
  }

  /** What pinning the mean down would cost, from the parts made so far. */
  private writeCost(): void {
    const c = [...this.cards].reverse().find((k) => k.made.length >= 5);
    if (!c) {
      setHtml(this.cost, "Each dot is one part made, at its fundamental frequency over the design's. Make a few more and watch the spread appear.");
      return;
    }
    const r = c.made.map((p) => p.omega[0] / design(c.kind).omega[0]), m = mean(r), s = sd(r);
    const pm = (100 * 1.96 * s) / (m * Math.sqrt(r.length));
    const need = Math.ceil(((1.96 * s) / (1e-3 * m)) ** 2);
    const ms = median(c.made.map((p) => p.ms));
    setHtml(this.cost,
      `From ${r.length} ${c.kind}s, the mean ω₁ is <b>${m.toFixed(3)}×</b> the design's, give or take <b>${pm.toFixed(1)}%</b> (95%). ` +
      `To know it to ±0.1% would take about <b>${need.toLocaleString()}</b> ${c.kind}s. At ${ms.toFixed(1)} ms a solve on this ` +
      `coarse mesh, that is ${duration(need * ms)}; on a mesh fine enough to trust, each solve costs many times more.`);
  }
}

/** Rewritten only when it changes: every frame would undo a reader's text selection. */
function setHtml(e: HTMLElement, html: string): void {
  if (e.dataset.html !== html) { e.dataset.html = html; e.innerHTML = html; }
}

const GHOST = "220, 238, 255";

/**
 * The made part, and the deterministic one as a wireframe ghost around it. The
 * ghost is drawn twice: plainly before the part, so the part hides what lies
 * behind it, and faintly after, so what it hides still shows through. Its
 * dashed lines (the plate's nodal lines) are drawn plainly both times.
 */
function ghost(
  ctx: CanvasRenderingContext2D, pr: Projector, quads: readonly Quad[],
  wire: { readonly lines: readonly Vec3[][]; readonly dashed?: readonly Vec3[][] },
): void {
  const pass = (alpha: number, dashedAlpha: number) => {
    for (const [set, dash, a] of [[wire.lines, [], alpha], [wire.dashed ?? [], [4, 3], dashedAlpha]] as const) {
      ctx.strokeStyle = `rgba(${GHOST}, ${a})`;
      ctx.setLineDash(dash);
      ctx.lineWidth = dash.length ? 1.6 : 1.1;
      ctx.beginPath();
      for (const line of set) line.forEach((v, i) => (i ? ctx.lineTo(...pr.to(v)) : ctx.moveTo(...pr.to(v))));
      ctx.stroke();
    }
    ctx.setLineDash([]);
  };
  pass(0.7, 0.85);
  paintQuads(ctx, pr, quads);
  // The ghost's nodal lines stay plain through the part: they are what its own are compared with.
  pass(0.22, 0.85);
}

/** Size a canvas's backing store to its box; null while it has none. */
function fit(canvas: HTMLCanvasElement): CanvasRenderingContext2D | null {
  const dpr = window.devicePixelRatio || 1, w = Math.round(canvas.clientWidth * dpr), h = Math.round(canvas.clientHeight * dpr);
  if (!w || !h) return null;
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
  return canvas.getContext("2d");
}

function drawBeam(ctx: CanvasRenderingContext2D, cam: Camera, W: number, H: number, part: Part, d: Part, n: number, zp: number, zd: number): void {
  const A = 0.24, H0 = 0.06, B = 0.05, N = beamX.length;
  const bounds: Vec3[] = [];
  for (const x of [-0.08, 1.02]) for (const y of [-0.18, 0.18]) for (const z of [-0.3, 0.3]) bounds.push([x, y, z]);
  const pr = projector(cam, bounds, W, H);
  const shape = part.shapes[n], z = (i: number) => A * zp * shape[i];
  const half = (i: number) => (H0 / 2) * (1 + DEPTH_GAIN * (part.depth[i] - 1));
  const quads: Quad[] = [];
  // The wall it is clamped into.
  const wall: RGB = [52, 58, 70], wx = -0.07, s = 0.16;
  const box = (x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, rgb: RGB) => {
    quads.push(
      { p: [[x1, y0, z0], [x1, y1, z0], [x1, y1, z1], [x1, y0, z1]], rgb },
      { p: [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]], rgb },
      { p: [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]], rgb },
      { p: [[x0, y1, z0], [x1, y1, z0], [x1, y1, z1], [x0, y1, z1]], rgb },
    );
  };
  box(wx, 0, -s, s, -s, s, wall);
  for (let i = 0; i < N - 1; i++) {
    const x0 = beamX[i], x1 = beamX[i + 1], a = z(i), b = z(i + 1), ha = half(i), hb = half(i + 1);
    const rgb = depthColour((part.depth[i] + part.depth[i + 1]) / 2);
    quads.push(
      { p: [[x0, -B, a + ha], [x1, -B, b + hb], [x1, B, b + hb], [x0, B, a + ha]], rgb },
      { p: [[x0, -B, a - ha], [x1, -B, b - hb], [x1, B, b - hb], [x0, B, a - ha]], rgb },
      { p: [[x0, -B, a - ha], [x1, -B, b - hb], [x1, -B, b + hb], [x0, -B, a + ha]], rgb },
      { p: [[x0, B, a - ha], [x1, B, b - hb], [x1, B, b + hb], [x0, B, a + ha]], rgb },
    );
  }
  const e = N - 1, ze = z(e), he = half(e);
  quads.push({ p: [[1, -B, ze - he], [1, B, ze - he], [1, B, ze + he], [1, -B, ze + he]], rgb: depthColour(part.depth[e]) });

  // The deterministic beam, struck at the same moment: its four long edges, and a section every eighth of the span.
  const g = (i: number) => A * zd * d.shapes[n][i], h0 = H0 / 2;
  const corners: [number, number][] = [[-B, -h0], [B, -h0], [B, h0], [-B, h0]];
  const lines: Vec3[][] = corners.map(([y, dz]) => Array.from(beamX, (x, i): Vec3 => [x, y, g(i) + dz]));
  for (let i = 0; i < N; i += (N - 1) / 8)
    lines.push([...corners, corners[0]].map(([y, dz]): Vec3 => [beamX[i], y, g(i) + dz]));
  ghost(ctx, pr, quads, { lines });
}

function drawPlate(ctx: CanvasRenderingContext2D, cam: Camera, W: number, H: number, part: Part, d: Part, n: number, zp: number, zd: number): void {
  const A = 0.17, T0 = 0.035, G = PLATE_POINTS;
  const bounds: Vec3[] = [];
  for (const x of [-0.5, 0.5]) for (const y of [-0.5, 0.5]) for (const z of [-0.22, 0.2]) bounds.push([x, y, z]);
  const pr = projector(cam, bounds, W, H);
  const phi = part.shapes[n], psi = d.shapes[n], h = part.depth;
  const at = (ix: number, iy: number): Vec3 => [plateX[ix] - 0.5, plateX[iy] - 0.5, A * zp * phi[iy * G + ix]];
  const thick = (ix: number, iy: number) => T0 * (1 + DEPTH_GAIN * (h[iy * G + ix] - 1));
  const quads: Quad[] = [];
  for (let iy = 0; iy < G - 1; iy++)
    for (let ix = 0; ix < G - 1; ix++) {
      const c = [at(ix, iy), at(ix + 1, iy), at(ix + 1, iy + 1), at(ix, iy + 1)] as const;
      const k = [iy * G + ix, iy * G + ix + 1, (iy + 1) * G + ix + 1, (iy + 1) * G + ix];
      const rgb = depthColour((h[k[0]] + h[k[1]] + h[k[2]] + h[k[3]]) / 4);
      const mine = cellZeros([phi[k[0]], phi[k[1]], phi[k[2]], phi[k[3]]]);
      // A point of the cell's unit square, on the surface as it is now (bilinear between the corners).
      const on = (u: number, v: number): Vec3 => {
        const lerp = (j: number) => (1 - u) * (1 - v) * c[0][j] + u * (1 - v) * c[1][j] + u * v * c[2][j] + (1 - u) * v * c[3][j];
        return [lerp(0), lerp(1), lerp(2)];
      };
      // The made plate's nodal lines, on its surface, hidden with it.
      const then = mine.length ? (cx: CanvasRenderingContext2D, to: (v: Vec3) => [number, number]) => {
        cx.strokeStyle = "#fff"; cx.lineWidth = 2;
        cx.beginPath();
        for (const [u0, v0, u1, v1] of mine) { cx.moveTo(...to(on(u0, v0))); cx.lineTo(...to(on(u1, v1))); }
        cx.stroke();
      } : undefined;
      quads.push({ p: c, rgb, then });
    }
  // Its edges, as thick as it was made.
  const skirt = (cells: [number, number][]) => {
    for (let i = 0; i < cells.length - 1; i++) {
      const [ax, ay] = cells[i], [bx, by] = cells[i + 1], a = at(ax, ay), b = at(bx, by);
      const ta = thick(ax, ay), tb = thick(bx, by);
      quads.push({
        p: [a, b, [b[0], b[1], b[2] - tb], [a[0], a[1], a[2] - ta]],
        rgb: depthColour((h[ay * G + ax] + h[by * G + bx]) / 2),
      });
    }
  };
  const idx = Array.from({ length: G }, (_, i) => i);
  skirt(idx.map((i) => [i, 0]));
  skirt(idx.map((i) => [i, G - 1]));
  skirt(idx.map((i) => [0, i]));
  skirt(idx.map((i) => [G - 1, i]));

  // The deterministic plate, struck at the same moment: every third grid line, and its nodal lines (which sit at z = 0).
  const gz = (ix: number, iy: number): Vec3 => [plateX[ix] - 0.5, plateX[iy] - 0.5, A * zd * psi[iy * G + ix]];
  const lines: Vec3[][] = [];
  for (let j = 0; j < G; j += 3) {
    lines.push(idx.map((i) => gz(i, j)));
    lines.push(idx.map((i) => gz(j, i)));
  }
  const nodal: Vec3[][] = [];
  const px = (ix: number, u: number) => plateX[ix] + u * (plateX[ix + 1] - plateX[ix]) - 0.5;
  for (let iy = 0; iy < G - 1; iy++)
    for (let ix = 0; ix < G - 1; ix++) {
      const k = [iy * G + ix, iy * G + ix + 1, (iy + 1) * G + ix + 1, (iy + 1) * G + ix];
      for (const [u0, v0, u1, v1] of cellZeros([psi[k[0]], psi[k[1]], psi[k[2]], psi[k[3]]]))
        nodal.push([[px(ix, u0), px(iy, v0), 0], [px(ix, u1), px(iy, v1), 0]]);
    }
  ghost(ctx, pr, quads, { lines, dashed: nodal });
}

const mean = (r: readonly number[]) => r.reduce((s, v) => s + v, 0) / r.length;
const sd = (r: readonly number[]) => {
  const m = mean(r);
  return Math.sqrt(r.reduce((s, v) => s + (v - m) ** 2, 0) / Math.max(1, r.length - 1));
};
const median = (r: readonly number[]) => [...r].sort((a, b) => a - b)[r.length >> 1];
const clamp = (v: number, lo: number, hi: number) => Math.min(hi, Math.max(lo, v));
/** A fixed scatter in [0, 1) for dot i, so dots keep their place as others arrive. */
const hash = (i: number) => {
  const s = Math.sin(i * 12.9898 + 4.1414) * 43758.5453;
  return s - Math.floor(s);
};

function duration(ms: number): string {
  const s = ms / 1000;
  if (s < 90) return `${s.toFixed(0)} s of computing`;
  if (s < 5400) return `${(s / 60).toFixed(0)} minutes of computing`;
  if (s < 2 * 86400) return `${(s / 3600).toFixed(1)} hours of computing`;
  return `${(s / 86400).toFixed(1)} days of computing`;
}
