/**
 * The Nusselt trace, drawn.
 *
 * A 2-D canvas rather than a third WebGPU pass, and redrawn when a sample lands
 * rather than every frame: diagnostics arrive four times a second off the
 * frame's dependency chain (see `pollStats`), and putting a chart on the hot
 * loop's critical path to show a number that changed 250 ms ago would invert the
 * whole arrangement.
 *
 * The panel sits in the corner opposite the pane and below the readout, over the
 * *outer* boundary — which is the cold one, and therefore the dark end of the
 * colour map, so the one part of the field a panel can cover is the part with
 * nothing in it. It keeps 8% transparency so it still reads as an overlay on the
 * simulation, and no more than that: the series colours were validated for ≥ 3:1
 * contrast against the surface that transparency produces, including over field
 * brightnesses this corner cannot actually reach.
 *
 * There is no hover layer, which for an interactive chart is a deliberate
 * omission rather than an oversight. A crosshair would have to take pointer
 * events over the simulation to read a value that is already printed, live, in
 * the legend and again in the readout above — so the values are never gated
 * behind an interaction, which is what the tooltip rule is protecting.
 */

import { boundaryNames, type GeometryKind } from "../geometry";
import { dimensionalTime } from "./dimensional";
import {
  NU_COLOUR, NuTrace, axisDecimals, niceAxis, type NuSample, type NuSeries,
} from "./nusselt";

const TAU = 2 * Math.PI;

/** Series order: outer first, so the warm curve is drawn over the cool one. */
const DRAW: NuSeries[] = ["outer", "inner"];
/** Legend order: as the readout names them. */
const KEY: NuSeries[] = ["inner", "outer"];

/**
 * `l` is a floor — the y gutter is measured from the tick labels. `b` holds two
 * label rows, not one: the same two instants in nondimensional time above and in
 * years below. It is sized here rather than left to the canvas edge because a
 * container whose height excludes its own axis band is how a chart ends up with
 * its labels cropped.
 */
const PAD = { l: 8, r: 11, t: 5, b: 28 };
/** Centres of the two label rows, up from the bottom of the canvas. */
const ROW = { nondim: 19, dim: 7 };

/**
 * Stroke width and end-marker radius per series — deliberately *not* equal.
 *
 * These two curves converge on the same value, and that convergence is the
 * entire reading. Drawn at one width they render as a single line: the fluxes
 * agree to six digits at steady state, which on an 84 px axis is a ten-thousandth
 * of a pixel, so no honest rendering can separate them. Two 2 px lines therefore
 * leave the reader unable to tell a series that agrees from one that is not being
 * drawn at all — the plot's own subject, silently lost.
 *
 * So they are nested rather than laid side by side: `outer` wide, `inner` narrow
 * on top, which reads as a blue rim around an orange core exactly where the two
 * coincide and as two ordinary lines wherever they do not. The mean is the 2 px
 * line spec, and nothing is displaced — the same x, the same y, the same value.
 *
 * The gap is as small as it can be and still show. Narrowing the display window
 * separates the curves, and at that zoom the asymmetry has no job left to do —
 * it is just one series drawn heavier than the other, which reads as emphasis
 * nothing intends. 0.55 px of rim per side is what survives that trade: enough
 * to see against the core at any device pixel ratio, little enough that two
 * separated curves look like a pair rather than a hierarchy.
 */
const WIDTH: Record<NuSeries, number> = { outer: 2.8, inner: 1.7 };
const DOT: Record<NuSeries, number> = { outer: 5.5, inner: 4 };
const RING = 2;      // surface ring, so the ends stay legible over the curves

/**
 * Chart chrome, as the reference palette's *roles* — muted label, hairline grid,
 * baseline, surface — restated in this page's own ink. The HUD is `#cfe` on
 * near-black, and the palette's warm grays beside it read as a different design
 * pasted in. The muted step is held at 70% because the panel is translucent:
 * at 55% the tick labels fall under 4.5:1 wherever the field brightens behind
 * them.
 */
const INK = "rgba(207, 238, 255, 0.70)";    // tick labels
/**
 * The dimensional row, one step back from the nondimensional one above it — but
 * only one step. It is subordinate in *position* and by carrying its unit, not by
 * being faint: it is a small-text label like any other and has to clear 4.5:1
 * against the surface the panel's transparency can produce, which rules out the
 * 45% that would have read as properly secondary (3.8:1). 58% measures 5.4:1.
 */
const DIM = "rgba(207, 238, 255, 0.58)";
const GRID = "rgba(207, 238, 255, 0.12)";   // gridlines
const AXIS = "rgba(207, 238, 255, 0.22)";   // baseline
const SURFACE = "#07070e";                  // the panel's own ground, for the rings
const FONT = "10px ui-monospace, monospace";

const el = <K extends keyof HTMLElementTagNameMap>(
  tag: K, cls?: string,
): HTMLElementTagNameMap[K] => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  return e;
};

export class NusseltPlot {
  readonly trace = new NuTrace();
  private readonly canvas = el("canvas");
  private readonly ctx: CanvasRenderingContext2D | null;
  private readonly value: Record<NuSeries, HTMLElement>;
  private readonly name: Record<NuSeries, HTMLElement>;
  private window = Infinity;
  /** Which domain the trace is of: it sets both the clock and the two names. */
  private kind: GeometryKind = "annulus";
  private dpr = 1;
  private w = 0;
  private h = 0;

  /**
   * Builds the figure inside `root` — caption, canvas, legend — rather than
   * having the page declare it, so `NU_COLOUR` is the single source the legend
   * keys and the curves both read. `#pane` is mounted the same way.
   */
  constructor(root: HTMLElement, window = Infinity, kind: GeometryKind = "annulus") {
    this.window = window;
    this.kind = kind;
    const caption = el("figcaption");
    caption.textContent = "Nusselt number vs time";

    const key = el("ul", "nu-key");
    const value = {} as Record<NuSeries, HTMLElement>;
    const name = {} as Record<NuSeries, HTMLElement>;
    for (const s of KEY) {
      const row = el("li");
      const swatch = el("i");
      swatch.style.background = NU_COLOUR[s];
      // The key is a sample of the mark, weight included: the two curves are
      // drawn at different widths on purpose, and a legend that showed them
      // equal would leave the wide blue rim around the orange core unexplained.
      swatch.style.height = `${WIDTH[s].toFixed(1)}px`;
      name[s] = el("span");
      value[s] = el("b");
      value[s].textContent = "—";
      row.append(swatch, name[s], value[s]);
      key.append(row);
    }
    this.value = value;
    this.name = name;
    this.label();

    // The curve is not reachable through a canvas whatever is announced about
    // it, and everything it *quantifies* is in the legend below as live text —
    // which is the reason those values sit there rather than as labels on the
    // curve ends. So the figure is the caption plus that legend, and the canvas
    // is hidden rather than announced as an unnamed graphic.
    this.canvas.setAttribute("aria-hidden", "true");
    root.append(caption, this.canvas, key);
    this.ctx = this.canvas.getContext("2d");

    new ResizeObserver(() => this.measure()).observe(this.canvas);
    this.measure();
  }

  /** Record one poll; redraw only if it added a sample. */
  push(sample: NuSample): void {
    if (this.trace.push(sample)) this.draw();
  }

  /**
   * Span to display, in solver steps; `Infinity` for the whole trace. Only the
   * view changes — the buffer keeps everything either way, so widening the window
   * again brings the earlier history back rather than starting from now.
   */
  setWindow(steps: number): void {
    if (steps === this.window) return;
    this.window = steps;
    this.draw();
  }

  /**
   * Adopt a rebuilt run's geometry: it renames the two series and rescales the
   * dimensional axis row, since a box's unit of time is a diffusion time across
   * its depth and the annulus' is one across its outer radius (`dimensional.ts`).
   * Paired with `clear` by the caller — a trace does not survive a rebuild.
   */
  setGeometry(kind: GeometryKind): void {
    if (kind === this.kind) return;
    this.kind = kind;
    this.label();
    this.draw();
  }

  private label(): void {
    const n = boundaryNames(this.kind);
    for (const s of KEY) this.name[s].textContent = n[s];
  }

  /** Drop the window — a new run's samples do not continue the old run's curve. */
  clear(): void {
    this.trace.clear();
    this.draw();
  }

  /** Match the backing store to the element, then redraw at the new size. */
  private measure(): void {
    this.dpr = Math.min(devicePixelRatio || 1, 2);
    this.w = this.canvas.clientWidth;
    this.h = this.canvas.clientHeight;
    this.canvas.width = Math.max(1, Math.round(this.w * this.dpr));
    this.canvas.height = Math.max(1, Math.round(this.h * this.dpr));
    this.draw();
  }

  /** Hairlines land on device pixels, or a 1px rule renders as a 2px smear. */
  private snap(v: number): number {
    return (Math.round(v * this.dpr) + 0.5) / this.dpr;
  }

  /**
   * One row of the x axis: the window's start at the left, its end at the right.
   *
   * The pair is measured before it is drawn, and the *left* label is dropped if
   * the two would not both fit — this panel can be 150 px wide, and "1.29 Tyr"
   * twice does not go into it. The right-hand one survives because it is the
   * current instant, which is the value being read; the left is the window's age,
   * recoverable from the window setting. Neither is ever clipped, and nothing is
   * nudged apart to make room.
   *
   * It is also dropped when the two *round to the same string*, which the short
   * windows do: three significant figures cannot separate the ends of a 500-step
   * window, and the same value printed at both ends of an axis reads as a fault
   * rather than as a narrow span.
   */
  private pair(
    ctx: CanvasRenderingContext2D, x0: number, x1: number, y: number,
    left: string, right: string, colour: string, wide: boolean,
  ): void {
    ctx.fillStyle = colour;
    ctx.textAlign = "right";
    ctx.fillText(right, x1, y);
    if (!wide || left === right) return;
    const room = x1 - ctx.measureText(right).width - 8;
    if (x0 + ctx.measureText(left).width <= room) {
      ctx.textAlign = "left";
      ctx.fillText(left, x0, y);
    }
  }

  draw(): void {
    const ctx = this.ctx;
    // Below this there is no plot left under the axis band, only the band.
    if (!ctx || this.w < 40 || this.h < PAD.b + PAD.t + 12) return;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, this.w, this.h);
    ctx.font = FONT;
    ctx.textBaseline = "middle";

    const last = this.trace.last;
    for (const s of KEY)
      this.value[s].textContent = last ? last[s].toFixed(4) : "—";

    // The display window, resolved against the trace before anything is scaled:
    // every extent, tick and coordinate below is over the visible slice only.
    const n = this.trace.length;
    const from = this.trace.first(this.window);
    const ex = this.trace.extent(from);
    const hair = 1 / this.dpr;
    const y1 = this.h - PAD.b;
    if (!ex) {
      // Nothing polled yet — draw the baseline the plot will hang from and stop.
      // It lasts a quarter of a second; an empty frame reads better than a
      // caption explaining that.
      ctx.strokeStyle = AXIS;
      ctx.lineWidth = hair;
      ctx.beginPath();
      ctx.moveTo(PAD.l, this.snap(y1));
      ctx.lineTo(this.w - PAD.r, this.snap(y1));
      ctx.stroke();
      return;
    }

    // The gutter is measured, not assumed: Nu runs from ~1 to well past 100
    // across the Ra range, and "112.5" is twice the width of "1.0".
    const ax = niceAxis(ex.lo, ex.hi);
    const labels = ax.ticks.map((v) => v.toFixed(ax.decimals));
    const gutter = labels.reduce((m, s) => Math.max(m, ctx.measureText(s).width), 0);
    const x0 = PAD.l + (gutter > 0 ? gutter + 5 : 0);
    const x1 = this.w - PAD.r;
    const y0 = PAD.t;
    if (x1 - x0 < 24 || y1 - y0 < 12) return;

    const Y = (v: number) => y1 - ((v - ax.lo) / (ax.hi - ax.lo)) * (y1 - y0);
    // A single sample, or a paused run, has no time span: pin it to the right
    // edge, where the newest sample always is.
    const span = ex.t1 - ex.t0;
    const X = (t: number) => (span > 0 ? x0 + ((t - ex.t0) / span) * (x1 - x0) : x1);

    ctx.strokeStyle = GRID;
    ctx.lineWidth = hair;
    ctx.fillStyle = INK;
    ctx.textAlign = "right";
    ax.ticks.forEach((v, i) => {
      const y = this.snap(Y(v));
      ctx.beginPath();
      ctx.moveTo(x0, y);
      ctx.lineTo(x1, y);
      ctx.stroke();
      ctx.fillText(labels[i], x0 - 5, Y(v));
    });

    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(x0, this.snap(y1));
    ctx.lineTo(x1, this.snap(y1));
    ctx.stroke();

    // The two ends of the window, twice: nondimensional above, dimensional below.
    // The same instants in both, one under the other, rather than one unit on the
    // axis and the other in a footer — the reader is converting a *coordinate* on
    // this axis, and the conversion should be legible without arithmetic.
    //
    // The nondimensional values resolve to the span rather than to the gridline
    // step: these are endpoint values, and the step that suits a tick column is
    // far too coarse for them — a window over 0.70 … 4.80 would print as
    // "t 1 … 5". See `axisDecimals` for the one case where the span alone is not
    // enough.
    const td = axisDecimals(ex.t0, span);
    this.pair(ctx, x0, x1, this.h - ROW.nondim,
      `t ${ex.t0.toFixed(td)}`, ex.t1.toFixed(td), INK, span > 0);
    this.pair(ctx, x0, x1, this.h - ROW.dim,
      dimensionalTime(ex.t0), dimensionalTime(ex.t1),
      DIM, span > 0);

    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    for (const s of DRAW) {
      ctx.lineWidth = WIDTH[s];
      ctx.strokeStyle = NU_COLOUR[s];
      ctx.beginPath();
      for (let i = from; i < n; i++) {
        const p = this.trace.at(i);
        const x = X(p.t), y = Y(p[s]);
        if (i === from) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    // End markers, ringed in the surface colour so they stay legible over the
    // curves they terminate. Every ring goes down *before* any fill: ringing each
    // marker as it is drawn would put a surface ring between the two when they
    // coincide, erasing the very rim that shows both are there.
    if (!last) return;
    const xe = X(last.t);
    ctx.lineWidth = RING;
    ctx.strokeStyle = SURFACE;
    for (const s of DRAW) {
      ctx.beginPath();
      ctx.arc(xe, Y(last[s]), DOT[s] + RING / 2, 0, TAU);
      ctx.stroke();
    }
    for (const s of DRAW) {
      ctx.beginPath();
      ctx.arc(xe, Y(last[s]), DOT[s], 0, TAU);
      ctx.fillStyle = NU_COLOUR[s];
      ctx.fill();
    }
  }
}
