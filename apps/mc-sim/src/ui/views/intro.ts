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
 * Below each, every part of its kind made so far, as a dot at its ω₁ over the
 * design's: a distribution filling in, which is all a Monte Carlo estimate is —
 * and a table of what the parts so far say, and what pinning the mean down would
 * cost at this mesh's measured
 * solve time.
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
  /** The kind's batch: a dot per part made, and a table of what they say and what pinning their mean down would cost. */
  readonly strip: HTMLCanvasElement;
  readonly table: HTMLElement;
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
      "treat each part as random and estimate its statistics by simulation: its mean frequency, the standard deviation of that frequency, and the " +
      "chance it resonates with whatever shakes it.");
    el("p", "in-lead", head,
      "Below, you run the factory. Each press of a button makes a new beam or plate, simulates it, and shows it " +
      "vibrating. Its colour shows where it came out thinner (orange) or thicker (blue) than drawn. The pale wireframe " +
      "around it is the part as designed, vibrating alongside: watch the two drift out of step, because the made part " +
      "rings at a different frequency. On the plate, the white lines mark where the surface stays still; compare them " +
      "with the design's dashed ones. The mode buttons switch between each part's simplest ways of vibrating, and you " +
      "can drag either part to turn it. The differences are exaggerated so they are easy to see.");

    const stage = el("div", "in-stage", root);
    // The batches sit in their own row, under the key and a line on how to read them, each in its part's column.
    const key = el("div", "in-key");
    const dots = el("p", "in-dots", undefined,
      "Each dot is one part made, placed by how far its lowest resonant frequency lands from the design's. " +
      "Make a few more and watch the spread appear: the table under each plot keeps track of what the parts made so " +
      "far tell us, and of what it would take to pin their average down.");
    const batches = el("div", "in-stage in-batches");
    for (const kind of ["beam", "plate"] as const) this.cards.push(this.card(stage, batches, kind));
    root.append(key, dots, batches);

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


    const next = el("footer", "in-next", root);
    const ask = el("h2", "in-ask", next,
      "Simulating a detailed model that accounts for the uncertainty in its parts can take a ");
    el("em", "", ask, "very");
    ask.append(" long time. So what can we do?");
    // Paragraphs with sub- and superscripts, so written as HTML; all of it is this file's own text.
    const para = (cls: string, html: string) => { el("p", cls, next).innerHTML = html; };
    para("in-lead",
      "Monte Carlo is the honest answer: simulate <i>N</i> virtual parts and average what they do. Two errors stand " +
      "between that average and the truth. The first is <b>statistical</b>: the average of <i>N</i> random parts " +
      "wanders by about σ/√<i>N</i>, where σ is the standard deviation of the resonant frequency from part to part, so each extra digit of accuracy costs a " +
      "hundred times as many parts. The second is the <b>model's own</b>: a simulation divides the part into a mesh " +
      "of small elements of size <i>h</i>, and its error is proportional to <i>h</i><sup>α</sup>, for some α &gt; 0. " +
      "Shrinking <i>h</i> makes each simulation more accurate, but also more expensive, and all <i>N</i> of them pay " +
      "that price. Needing both at " +
      "once is what makes plain Monte Carlo so slow.");
    para("in-lead",
      "Can we do better? Most of the difference between one part and the next already shows on a coarse, cheap mesh; " +
      "a finer mesh only adds a small correction. So build a <b>hierarchy of meshes</b>, coarse to fine, " +
      "<i>h</i><sub>0</sub> &gt; <i>h</i><sub>1</sub> &gt; … &gt; <i>h</i><sub><i>L</i></sub>, and write the answer " +
      "on the finest as the answer on the coarsest plus a chain of corrections:");
    para("in-eq",
      "E[<i>Q</i><sub><i>L</i></sub>] = E[<i>Q</i><sub>0</sub>] + E[<i>Q</i><sub>1</sub> − <i>Q</i><sub>0</sub>] + … + " +
      "E[<i>Q</i><sub><i>L</i></sub> − <i>Q</i><sub><i>L</i>−1</sub>]");
    para("in-lead",
      "Here <i>Q</i><sub>ℓ</sub> is what a part does when simulated on mesh ℓ, and E[·] its average over all parts. " +
      "Spend many samples on the cheap coarse mesh, and only a few on each correction: corrections vary little from " +
      "part to part, so a few are enough. This is <b>multilevel Monte Carlo</b>, and it reaches the same accuracy for " +
      "a fraction of the work. The guided tours build up to it, one idea at a time:");
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

  private card(stage: HTMLElement, batches: HTMLElement, kind: Kind): Card {
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
    const batch = el("div", "in-batch", batches);
    const strip = el("canvas", "in-strip", batch);
    const table = el("div", "in-table", batch);
    const cam: Camera = kind === "beam" ? { yaw: -0.42, pitch: 0.32 } : { yaw: 0.55, pitch: 0.62 };
    const c: Card = { kind, canvas, stat, make, more, chips, strip, table, cam, made: [], mode: kind === "beam" ? 0 : 1, t0: 0 };
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
    this.cards.forEach((c, i) => {
      this.drawCard(c, t);
      drawStrip(c, SLOT[i]);
      writeTable(c);
    });
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

}

/** One kind's batch: a dot per part made at its ω₁ over the design's, with the mean and its 95% interval. */
function drawStrip(c: Card, colour: string): void {
  const ctx = fit(c.strip);
  if (!ctx) return;
  const dpr = window.devicePixelRatio || 1, W = c.strip.width / dpr, H = c.strip.height / dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, W, H);
  const lo = 0.7, hi = 1.3, left = 6, right = 6, X = (r: number) => left + ((r - lo) / (hi - lo)) * (W - left - right);
  const top = 20, axisY = H - 16, rowH = axisY - top;
  ctx.font = "11px ui-monospace, monospace";
  ctx.textBaseline = "middle";
  ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
  ctx.textAlign = "left";
  ctx.fillText(`${c.made.length} ${c.kind}${c.made.length === 1 ? "" : "s"} made`, left, 8);
  ctx.textAlign = "right";
  ctx.fillText("resonant frequency vs. design", W - right, 8);
  // The design's line, and the axis.
  ctx.strokeStyle = "rgba(207, 238, 255, 0.35)";
  ctx.setLineDash([3, 3]);
  ctx.beginPath(); ctx.moveTo(X(1), top); ctx.lineTo(X(1), axisY); ctx.stroke();
  ctx.setLineDash([]);
  ctx.strokeStyle = "rgba(207, 238, 255, 0.25)";
  ctx.beginPath(); ctx.moveTo(left, axisY); ctx.lineTo(W - right, axisY); ctx.stroke();
  ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
  ctx.textAlign = "center";
  for (let r = lo; r <= hi + 1e-9; r += 0.1) {
    ctx.fillRect(X(r) - 0.5, axisY, 1, 4);
    const pc = Math.round((r - 1) * 100);
    ctx.fillText(pc === 0 ? "design" : `${pc > 0 ? "+" : "−"}${Math.abs(pc)}%`, clamp(X(r), left + 14, W - right - 14), axisY + 11);
  }
  const r = c.made.map((p) => p.omega[0] / design(c.kind).omega[0]);
  const yOf = (i: number) => top + 6 + hash(i) * (rowH - 12);
  // Mean and its 95% interval, once there are two to spread.
  if (r.length > 1) {
    const m = mean(r), half = (1.96 * sd(r)) / Math.sqrt(r.length);
    ctx.fillStyle = "rgba(207, 238, 255, 0.10)";
    ctx.fillRect(X(m - half), top, Math.max(1, X(m + half) - X(m - half)), rowH - 2);
    ctx.fillStyle = "rgba(207, 238, 255, 0.9)";
    ctx.fillRect(X(m) - 0.75, top, 1.5, rowH - 2);
  }
  ctx.globalAlpha = 0.75;
  ctx.fillStyle = colour;
  r.forEach((v, i) => { ctx.beginPath(); ctx.arc(X(clamp(v, lo, hi)), yOf(i), 3.2, 0, 2 * Math.PI); ctx.fill(); });
  ctx.globalAlpha = 1;
  if (r.length) {
    const i = r.length - 1;
    ctx.strokeStyle = "#fff"; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(X(clamp(r[i], lo, hi)), yOf(i), 5.5, 0, 2 * Math.PI); ctx.stroke();
  }
}

/**
 * What one kind's parts made so far say: their average resonant frequency
 * against the design's and how sure that average is, and what pinning it down
 * to ±0.1% would cost at the measured time of one simulation. Dashes until
 * there are two parts to spread.
 */
function writeTable(c: Card): void {
  const n = c.made.length, r = c.made.map((p) => p.omega[0] / design(c.kind).omega[0]);
  const ms = median(c.made.map((p) => p.ms));
  let avg = "—", pm = "—", need = "—", total = "—";
  if (n > 1) {
    const m = mean(r), s = sd(r), N = Math.ceil(((1.96 * s) / (1e-3 * m)) ** 2);
    avg = `${Math.abs(100 * (m - 1)).toFixed(1)}% ${m < 1 ? "below" : "above"}`;
    pm = `±${((100 * 1.96 * s) / (m * Math.sqrt(n))).toFixed(1)}%`;
    need = N.toLocaleString();
    total = duration(N * ms);
  }
  const rows: [string, string][] = [
    [`${c.kind}s made`, n.toLocaleString()],
    ["average resonant frequency, against the design's", avg],
    ["how sure that average is (95% confidence)", pm],
    [`${c.kind}s needed to know the average to ±0.1%`, need],
    ["time to simulate one, on this coarse model", `${ms.toFixed(1)} ms`],
    [`time to simulate all the ${c.kind}s needed`, total],
  ];
  setHtml(c.table,
    `<table>${rows.map(([k, v]) => `<tr><th>${k}</th><td>${v}</td></tr>`).join("")}</table>`);
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
  if (s < 90) return `${s.toFixed(0)} s`;
  if (s < 5400) return `${(s / 60).toFixed(0)} minutes`;
  if (s < 2 * 86400) return `${(s / 3600).toFixed(1)} hours`;
  return `${(s / 86400).toFixed(1)} days`;
}
