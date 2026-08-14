/**
 * The RMS velocity trace, drawn — `nuplot.ts`'s twin, with `v`/`vs` in place
 * of Nu's inner/outer.
 *
 * It sits in the opposite corner from the Nusselt panel: bottom-right, clear
 * of the pane (top-right) and the readout (top-left) for the same reason `#nu`
 * is where it is — a corner is the only place on a frame this square that
 * is not covering the field. The two panels are deliberately not merged into
 * one wider one along the bottom: Nu's whole reading is two curves converging,
 * and giving that plot a neighbour on the same axis would invite comparing
 * curves that share nothing but a time axis, which is exactly the reading a
 * second y scale on the *same* panel would produce (see `nusselt.ts`'s header
 * on why Nu keeps one axis for its own two series).
 *
 * The two curves drawn *here* do share a *scale*, and for the same reason
 * Nu's two series do: `v` and `vs` are the same reduction, in the same
 * code-unit velocity, taken over the domain and over the top boundary alone
 * (`rms.ts`'s header) — one pixel-per-code-unit mapping, so the picture never
 * claims the surface is faster or slower than the interior than it actually
 * is. What differs from Nu is that this shared quantity has a unit a reader
 * outside the code actually reasons in — cm/yr, a real plate speed — and the
 * left axis's code units do not name it. So the *right* axis is not a second
 * scale: it is the same linear mapping, relabelled in cm/yr via
 * `velocityUnitCmPerYear`, for the one curve (`vs`) that unit belongs to.
 * Its own ticks fall at different pixels than the left axis's, exactly
 * because "nice" round numbers in one unit are not round numbers in the
 * other — the gridlines stay the left axis's alone, or the panel would carry
 * two disagreeing grids for one plotted line.
 *
 * A 2-D canvas, redrawn on a new sample, off the frame's dependency chain —
 * same reasoning as `NusseltPlot`, same poll.
 */

import { dimensionalTime, dimensionalVelocity, velocityUnitCmPerYear } from "./dimensional";
import { decimalsFor, niceStep } from "./nusselt";
import {
  axisDecimals, niceAxis, RMS_COLOUR, RMS_SURFACE_COLOUR, RmsTrace, type RmsSample,
} from "./rms";

const TAU = 2 * Math.PI;

/** Which series a leg of the legend or a curve belongs to. */
type RmsSeries = "domain" | "surface";
const SERIES: RmsSeries[] = ["domain", "surface"];
const COLOUR: Record<RmsSeries, string> = { domain: RMS_COLOUR, surface: RMS_SURFACE_COLOUR };
/** Fixed rather than named from `boundaryNames`: "outer"/"top" is a detail of
 * which geometry `vs` happens to be read on, not the reading itself. */
const NAME: Record<RmsSeries, string> = { domain: "v_rms", surface: "surface v_rms" };

const PAD = { l: 8, r: 11, t: 5, b: 28 };
const ROW = { nondim: 19, dim: 7 };

const WIDTH = 2.4;
const DOT = 5;
const RING = 2;      // surface ring, so the end markers stay legible over the curves

const INK = "rgba(207, 238, 255, 0.70)";    // tick labels
/** Right axis's own tick labels — `RMS_SURFACE_COLOUR` (#8b5cf6) at `INK`'s
 * alpha, so the scale reads as the violet curve's own rather than a second
 * neutral axis a reader has to match to a curve by position alone. */
const INK_SURFACE = "rgba(139, 92, 246, 0.70)";
const DIM = "rgba(207, 238, 255, 0.58)";    // dimensional row
const GRID = "rgba(207, 238, 255, 0.12)";   // gridlines
const AXIS = "rgba(207, 238, 255, 0.22)";   // baseline
const SURFACE = "#07070e";                  // the panel's own ground, for the ring
const FONT = "10px ui-monospace, monospace";

const el = <K extends keyof HTMLElementTagNameMap>(
  tag: K, cls?: string,
): HTMLElementTagNameMap[K] => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  return e;
};

export class RmsPlot {
  readonly trace = new RmsTrace();
  private readonly canvas = el("canvas");
  private readonly ctx: CanvasRenderingContext2D | null;
  private readonly value: Record<RmsSeries, HTMLElement>;
  private window = Infinity;
  private dpr = 1;
  private w = 0;
  private h = 0;

  constructor(root: HTMLElement, window = Infinity) {
    this.window = window;
    const caption = el("figcaption");
    caption.textContent = "RMS velocity vs time";

    const key = el("ul", "rms-key");
    const value = {} as Record<RmsSeries, HTMLElement>;
    for (const s of SERIES) {
      const row = el("li");
      const swatch = el("i");
      swatch.style.background = COLOUR[s];
      swatch.style.height = `${WIDTH.toFixed(1)}px`;
      const name = el("span");
      name.textContent = NAME[s];
      value[s] = el("b");
      value[s].textContent = "—";
      row.append(swatch, name, value[s]);
      key.append(row);
    }
    this.value = value;

    // As in `NusseltPlot`: the canvas carries no information the legend below
    // does not already say as live text, so it is hidden rather than announced
    // as an unnamed graphic.
    this.canvas.setAttribute("aria-hidden", "true");
    root.append(caption, this.canvas, key);
    this.ctx = this.canvas.getContext("2d");

    new ResizeObserver(() => this.measure()).observe(this.canvas);
    this.measure();
  }

  /** Record one poll; redraw only if it added a sample. */
  push(sample: RmsSample): void {
    if (this.trace.push(sample)) this.draw();
  }

  /** Span to display, in solver steps; `Infinity` for the whole trace. */
  setWindow(steps: number): void {
    if (steps === this.window) return;
    this.window = steps;
    this.draw();
  }

  /** Drop the window — a new run's samples do not continue the old run's curve. */
  clear(): void {
    this.trace.clear();
    this.draw();
  }

  private measure(): void {
    this.dpr = Math.min(devicePixelRatio || 1, 2);
    this.w = this.canvas.clientWidth;
    this.h = this.canvas.clientHeight;
    this.canvas.width = Math.max(1, Math.round(this.w * this.dpr));
    this.canvas.height = Math.max(1, Math.round(this.h * this.dpr));
    this.draw();
  }

  private snap(v: number): number {
    return (Math.round(v * this.dpr) + 0.5) / this.dpr;
  }

  /** One row of the x axis — see `NusseltPlot.pair` for the rules this follows. */
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
    if (!ctx || this.w < 40 || this.h < PAD.b + PAD.t + 12) return;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, this.w, this.h);
    ctx.font = FONT;
    ctx.textBaseline = "middle";

    const last = this.trace.last;
    this.value.domain.textContent = last ? last.v.toFixed(4) : "—";
    this.value.surface.textContent = last ? dimensionalVelocity(last.vs) : "—";

    const n = this.trace.length;
    const from = this.trace.first(this.window);
    const ex = this.trace.extent(from);
    const hair = 1 / this.dpr;
    const y1 = this.h - PAD.b;
    if (!ex) {
      ctx.strokeStyle = AXIS;
      ctx.lineWidth = hair;
      ctx.beginPath();
      ctx.moveTo(PAD.l, this.snap(y1));
      ctx.lineTo(this.w - PAD.r, this.snap(y1));
      ctx.stroke();
      return;
    }

    const ax = niceAxis(ex.lo, ex.hi);
    const labels = ax.ticks.map((v) => v.toFixed(ax.decimals));
    const gutter = labels.reduce((m, s) => Math.max(m, ctx.measureText(s).width), 0);
    const x0 = PAD.l + (gutter > 0 ? gutter + 5 : 0);

    // The right axis: `ax`'s own bounds converted to cm/yr, not a second
    // extent of its own — the two axes are one linear scale wearing two
    // labels, so both cover exactly the range the curves are drawn over.
    // `niceStep` (not `niceAxis`) picks the tick spacing without padding the
    // range a second time; `ax` already carries the padding this range needs.
    const cmLo = ax.lo * velocityUnitCmPerYear, cmHi = ax.hi * velocityUnitCmPerYear;
    const cmStep = niceStep(cmHi - cmLo, 3);
    const cmDecimals = decimalsFor(cmStep);
    const cmTicks: number[] = [];
    for (let k = Math.ceil(cmLo / cmStep); cmTicks.length < 64 && k * cmStep <= cmHi; k++)
      cmTicks.push(k * cmStep);
    // No unit on the axis itself, matching the left axis's own bare numbers —
    // it is named where the value is, in the legend's `dimensionalVelocity`
    // readout below, and the violet ink already says which curve this scale
    // belongs to.
    const cmLabels = cmTicks.map((v) => v.toFixed(cmDecimals));
    const cmGutter = cmLabels.reduce((m, s) => Math.max(m, ctx.measureText(s).width), 0);
    const x1 = this.w - PAD.r - (cmGutter > 0 ? cmGutter + 5 : 0);
    const y0 = PAD.t;
    if (x1 - x0 < 24 || y1 - y0 < 12) return;

    const Y = (v: number) => y1 - ((v - ax.lo) / (ax.hi - ax.lo)) * (y1 - y0);
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

    ctx.fillStyle = INK_SURFACE;
    ctx.textAlign = "left";
    cmTicks.forEach((v, i) => ctx.fillText(cmLabels[i], x1 + 5, Y(v / velocityUnitCmPerYear)));

    ctx.strokeStyle = AXIS;
    ctx.beginPath();
    ctx.moveTo(x0, this.snap(y1));
    ctx.lineTo(x1, this.snap(y1));
    ctx.stroke();

    const td = axisDecimals(ex.t0, span);
    this.pair(ctx, x0, x1, this.h - ROW.nondim,
      `t ${ex.t0.toFixed(td)}`, ex.t1.toFixed(td), INK, span > 0);
    this.pair(ctx, x0, x1, this.h - ROW.dim,
      dimensionalTime(ex.t0), dimensionalTime(ex.t1),
      DIM, span > 0);

    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    ctx.lineWidth = WIDTH;
    const at = (s: RmsSeries, p: RmsSample) => (s === "domain" ? p.v : p.vs);
    for (const s of SERIES) {
      ctx.strokeStyle = COLOUR[s];
      ctx.beginPath();
      for (let i = from; i < n; i++) {
        const p = this.trace.at(i);
        const x = X(p.t), y = Y(at(s, p));
        if (i === from) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    if (!last) return;
    const xe = X(last.t);
    ctx.lineWidth = RING;
    ctx.strokeStyle = SURFACE;
    for (const s of SERIES) {
      ctx.beginPath();
      ctx.arc(xe, Y(at(s, last)), DOT + RING / 2, 0, TAU);
      ctx.stroke();
    }
    for (const s of SERIES) {
      ctx.beginPath();
      ctx.arc(xe, Y(at(s, last)), DOT, 0, TAU);
      ctx.fillStyle = COLOUR[s];
      ctx.fill();
    }
  }
}
