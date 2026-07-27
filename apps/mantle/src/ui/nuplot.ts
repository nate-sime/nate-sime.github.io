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

import { NU_COLOUR, NuTrace, decimalsFor, niceAxis, type NuSeries } from "./nusselt";

const TAU = 2 * Math.PI;

/** Series order: outer first, so the warm curve is drawn over the cool one. */
const DRAW: NuSeries[] = ["outer", "inner"];
/** Legend order: as the readout names them. */
const KEY: NuSeries[] = ["inner", "outer"];

const PAD = { l: 8, r: 11, t: 5, b: 14 };   // `l` is a floor; the gutter is measured

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
 */
const WIDTH: Record<NuSeries, number> = { outer: 3.2, inner: 1.5 };
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
  private dpr = 1;
  private w = 0;
  private h = 0;

  /**
   * Builds the figure inside `root` — caption, canvas, legend — rather than
   * having the page declare it, so `NU_COLOUR` is the single source the legend
   * keys and the curves both read. `#pane` is mounted the same way.
   */
  constructor(root: HTMLElement) {
    const caption = el("figcaption");
    caption.textContent = "Nusselt number vs time";

    const key = el("ul", "nu-key");
    const value = {} as Record<NuSeries, HTMLElement>;
    for (const s of KEY) {
      const row = el("li");
      const swatch = el("i");
      swatch.style.background = NU_COLOUR[s];
      // The key is a sample of the mark, weight included: the two curves are
      // drawn at different widths on purpose, and a legend that showed them
      // equal would leave the wide blue rim around the orange core unexplained.
      swatch.style.height = `${Math.round(WIDTH[s])}px`;
      const name = el("span");
      name.textContent = s;
      value[s] = el("b");
      value[s].textContent = "—";
      row.append(swatch, name, value[s]);
      key.append(row);
    }
    this.value = value;

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
  push(t: number, inner: number, outer: number): void {
    if (this.trace.push(t, inner, outer)) this.draw();
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

  draw(): void {
    const ctx = this.ctx;
    if (!ctx || this.w < 40 || this.h < 24) return;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.clearRect(0, 0, this.w, this.h);
    ctx.font = FONT;
    ctx.textBaseline = "middle";

    const last = this.trace.last;
    for (const s of KEY)
      this.value[s].textContent = last ? last[s].toFixed(4) : "—";

    const ex = this.trace.extent();
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

    // The window's own bounds, named on the left so the axis says what it is.
    // Resolved to a tenth of the span, not to the gridline step: these are
    // *endpoint values*, and the step that suits a tick column is far too coarse
    // for them — a window over 0.70 … 4.80 would print as "t 1 … 5". Scaling to
    // the span rather than to a fixed width covers t of a few 1e-3 just after a
    // reseed and of order 1 in a settled run.
    const td = decimalsFor(span / 10);
    ctx.fillStyle = INK;
    ctx.textAlign = "left";
    ctx.fillText(`t ${ex.t0.toFixed(td)}`, x0, this.h - PAD.b / 2);
    if (span > 0) {
      ctx.textAlign = "right";
      ctx.fillText(ex.t1.toFixed(td), x1, this.h - PAD.b / 2);
    }

    ctx.lineJoin = "round";
    ctx.lineCap = "round";
    const n = this.trace.length;
    for (const s of DRAW) {
      ctx.lineWidth = WIDTH[s];
      ctx.strokeStyle = NU_COLOUR[s];
      ctx.beginPath();
      for (let i = 0; i < n; i++) {
        const p = this.trace.at(i);
        const x = X(p.t), y = Y(p[s]);
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
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
