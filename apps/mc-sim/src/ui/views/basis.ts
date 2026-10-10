// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 1 on screen: every basis function of the chosen space (or its first or
 * second derivative), with the knots and their multiplicities under the axis;
 * and, to its right, the plate's basis built from it.
 *
 * Each element is drawn on its own, one-sided at both ends, and the pen lifts
 * between elements. A derivative that jumps at a knot therefore shows as two
 * ends that do not meet — which is what C^k means — rather than as a vertical
 * segment the space does not contain. All functions share one hue: what is on
 * show is their shape and overlap, not which is which (hover names one).
 *
 * The plate's space is the tensor product of this one with itself, on ne × ne
 * elements: n² functions B_i(x) B_j(y), p + 1 squared of them nonzero on each
 * element, C^k across every element edge. Too many to overlay, so the right
 * panel shows one, with its support, the mesh, and the whole control net —
 * the Greville points (ξ_i, ξ_j), one per function. Hovering a point of the
 * net picks its function, and the choice stays when the pointer leaves. The
 * derivative control differentiates in x: ∂ʳ/∂xʳ = B_i⁽ʳ⁾(x) B_j(y).
 */

import { basisOnElement, basisRow, firstDof, greville, splineSpace, type SplineSpace } from "../../spline";
import type { Figure } from "../figure";
import { colourBar, contour, linspace, meshLines, paint, shapeLimits, type Grid } from "../heatmap";
import { SLOT, type Plot, type PlotSpec, type Series } from "../plot";
import { continuityK, continuityName, type State } from "../state";
import { ro, type Block } from "../readout";
import type { ViewResult } from "./view";

const PER = 32;

export function renderBasis(fig: Figure, st: State): ViewResult {
  const [plot, plane] = fig.panels(2, { cols: 2 });
  const p = Math.max(1, st.p), k = continuityK(st.continuity, p), r = Math.min(st.derivative, p);
  const s = splineSpace(p, st.ne, k);
  drawPlane(plane, s, r);
  const P1 = p + 1;
  const xs: number[][] = Array.from({ length: s.n }, () => []);
  const ys: number[][] = Array.from({ length: s.n }, () => []);
  const sumX: number[] = [], sumY: number[] = [];
  for (let e = 0; e < s.ne; e++) {
    const a = s.breaks[e], b = s.breaks[e + 1], f0 = firstDof(s, e);
    for (let j = 0; j < P1; j++) {
      // A gap from the previous element this function lived on.
      if (xs[f0 + j].length) { xs[f0 + j].push(NaN); ys[f0 + j].push(NaN); }
    }
    sumX.push(NaN); sumY.push(NaN);
    for (let i = 0; i <= PER; i++) {
      const x = a + ((b - a) * i) / PER, N = basisOnElement(s, e, x, r);
      let sum = 0;
      for (let j = 0; j < P1; j++) {
        xs[f0 + j].push(x);
        ys[f0 + j].push(N[r * P1 + j]);
        sum += N[r * P1 + j];
      }
      sumX.push(x); sumY.push(sum);
    }
  }
  const prime = ["", "′", "″"][r];
  const series: Series[] = xs.map((x, i) => ({
    label: `B${sub(i)}${prime}`, x, y: ys[i], color: SLOT[0], width: 1.75, unlisted: true,
  }));
  series.push({
    label: r === 0 ? "Σ Bᵢ = 1" : `Σ Bᵢ${prime} = 0`, x: sumX, y: sumY, color: "rgba(207, 238, 255, 0.55)",
    width: 1.25, dash: [5, 4], unlisted: true, inert: true,
  });

  plot.draw({
    title: `${continuityName(k)} splines of degree ${p} — ${s.n} functions on ${s.ne} element${s.ne > 1 ? "s" : ""}`,
    xlabel: "x / L",
    ylabel: `Bᵢ${prime}(x)`,
    xlim: [0, 1],
    ylim: r === 0 ? [-0.06, 1.1] : undefined,
    series,
    hover: (se, i) => `${se.label}\nx/L = ${se.x[i].toFixed(4)}\nvalue = ${se.y[i].toPrecision(5)}`,
    over: (ctx, a) => {
      // Knot multiplicities, along the foot of the plot.
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      for (let i = 0; i <= s.ne; i++) {
        const X = a.sx(s.breaks[i]), mult = i === 0 || i === s.ne ? p + 1 : s.m;
        ctx.fillStyle = "#c98500";
        ctx.beginPath();
        ctx.arc(X, a.box.b, 3.5, 0, 2 * Math.PI);
        ctx.fill();
        if (s.ne <= 24 || i === 0 || i === s.ne) {
          ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
          ctx.fillText(`×${mult}`, X, a.box.b - 6);
        }
      }
    },
  });

  const knots = Array.from(s.U, (u) => +u.toFixed(4));
  const shown = knots.length > 40 ? `${knots.slice(0, 18).join(" ")} … ${knots.slice(-18).join(" ")}` : knots.join(" ");
  const line: Block[] = [
    ro.tiles([
      { label: "degree, continuity", value: `p = ${p}, ${continuityName(k)}`, detail: "across every interior knot" },
      { label: "interior knot multiplicity", value: `m = ${s.m}`, detail: "m = p − k" },
      { label: "basis functions", value: `n = ${s.n}`, detail: `p + 1 + (ne − 1)·m = ${p + 1} + ${s.ne - 1}·${s.m}` },
      { label: "nonzero on each element", value: String(P1), detail: "p + 1" },
    ], 2),
    ro.note(`knot vector [${shown}]`),
  ];
  if (st.continuity === "c1" && p < 2) line.push(ro.warn("degree 1 cannot be C¹: shown at C⁰"));
  if (st.derivative > p) line.push(ro.warn(`degree ${p} has no derivative of order ${st.derivative}: shown at ${p}`));
  if (k === 0) line.push(ro.warn("C⁰: fine for a string or a bar — not for a beam, whose energy needs w″ (see the beam and plate view)."));
  const plate: Block[] = [
    ro.lead(`B_i(x) B_j(y) on ${s.ne} × ${s.ne} elements, C${sup(k)} across every element edge`),
    ro.tiles([
      { label: "basis functions", value: `n² = ${s.n * s.n}`, detail: `${s.n}²` },
      { label: "nonzero on each element", value: String(P1 * P1), detail: "(p + 1)²" },
    ], 2),
    ro.note("one is shown, dashed round its support; hover a point of the control net to show its function"),
  ];
  return {
    readout: { sections: [{ title: "SPLINE BASIS ON THE LINE", blocks: line }, { title: "PLATE: TENSOR PRODUCT", blocks: plate }] },
    animate: false,
  };
}

const SUB = "₀₁₂₃₄₅₆₇₈₉";
const sub = (i: number) => String(i).split("").map((d) => SUB[+d]).join("");
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const sup = (i: number) => String(i).split("").map((d) => SUP[+d]).join("");

const RES = 129;
const INK = "rgba(207, 238, 255, 0.85)";
const NET = "rgba(207, 238, 255, 0.45)";
const MESH_INK = "rgba(5, 5, 12, 0.45)";
/** Above this many functions a side the net is hoverable but not drawn: its dots would cover the field. */
const NET_DRAWN = 24;

/** The function shown, as (i, j); kept across redraws and clamped when the space shrinks. */
let pick: { i: number; j: number } | null = null;

/** One function of the plate's basis, ∂ʳ/∂xʳ [B_i(x) B_j(y)], over the unit square. */
function drawPlane(plot: Plot, s: SplineSpace, r: number): void {
  const n = s.n, mid = Math.floor((n - 1) / 2);
  pick = pick ? { i: Math.min(pick.i, n - 1), j: Math.min(pick.j, n - 1) } : { i: mid, j: mid };
  const xs = linspace(0, 1, RES), xi = greville(s);
  // Every function along each axis, once: along(d)[a][f] = B_f⁽ᵈ⁾(x_a).
  const along = (d: number) => Array.from(xs, (x) => basisRow(s, x, d));
  const bx = along(r), by = r === 0 ? bx : along(0);
  const net: Series = {
    label: "control net", x: Array.from({ length: n * n }, (_, f) => xi[f % n]), y: Array.from({ length: n * n }, (_, f) => xi[Math.floor(f / n)]),
    color: NET, width: 0, markers: n <= NET_DRAWN, unlisted: true,
  };
  const spec = (i: number, j: number): PlotSpec => {
    const v = new Float64Array(RES * RES);
    let top = 0;
    for (let iy = 0; iy < RES; iy++)
      for (let ix = 0; ix < RES; ix++) {
        const w = bx[ix][i] * by[iy][j];
        v[iy * RES + ix] = w;
        top = Math.max(top, Math.abs(w));
      }
    top ||= 1;
    const g: Grid = { nx: RES, ny: RES, x0: 0, x1: 1, y0: 0, y1: 1, v };
    const scale = r === 0 ? "sequential" : "diverging";
    const d = ["", "∂/∂x ", "∂²/∂x² "][r];
    const [ux0, ux1, uy0, uy1] = [s.U[i], s.U[i + s.p + 1], s.U[j], s.U[j + s.p + 1]];
    return {
      title: `${d}B${sub(i)}(x) B${sub(j)}(y) — one of ${n * n} on ${s.ne} × ${s.ne} elements`,
      xlabel: "x / L", ylabel: "y / L", ...shapeLimits(plot, 1, 1),
      series: [net, { label: "shown", x: [xi[i]], y: [xi[j]], color: SLOT[1], width: 0, markers: true, unlisted: true, inert: true }],
      under: (ctx, a) => {
        paint(ctx, a, g, scale, top);
        meshLines(ctx, a, s.breaks, s.breaks, MESH_INK);
        if (r > 0) contour(ctx, a, g, 0, INK, 1.5);
        // The support: where the function is not zero.
        ctx.strokeStyle = INK;
        ctx.lineWidth = 1.5;
        ctx.setLineDash([5, 4]);
        ctx.strokeRect(a.sx(ux0), a.sy(uy1), a.sx(ux1) - a.sx(ux0), a.sy(uy0) - a.sy(uy1));
        ctx.setLineDash([]);
        ctx.strokeStyle = INK;
        ctx.lineWidth = 1;
        ctx.strokeRect(a.sx(0), a.sy(1), a.sx(1) - a.sx(0), a.sy(0) - a.sy(1));
      },
      over: (ctx, a) => colourBar(ctx, a, scale, top, r === 0 ? "B" : `${d.trim()}`, (u) => +u.toPrecision(2) + ""),
      hover: (se, f) => {
        const fi = f % n, fj = Math.floor(f / n);
        if (fi !== pick!.i || fj !== pick!.j) {
          pick = { i: fi, j: fj };
          // Redraw after this frame's hover is done, with the new function.
          queueMicrotask(() => plot.draw(spec(fi, fj)));
        }
        const span = (t: number) => `[${+s.U[t].toFixed(3)}, ${+s.U[t + s.p + 1].toFixed(3)}]`;
        return `B${sub(fi)}(x) B${sub(fj)}(y), Greville point (${se.x[f].toFixed(3)}, ${se.y[f].toFixed(3)})\nnonzero on ${span(fi)} × ${span(fj)}`;
      },
    };
  };
  plot.draw(spec(pick.i, pick.j));
}
