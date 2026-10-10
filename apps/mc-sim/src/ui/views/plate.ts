// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 7 on screen: one Kirchhoff plate, solved — its static deflection
 * under the chosen load on the left, and one of its modes, vibrating, on the
 * right.
 *
 * The deflection is a magnitude (the sequential scale) with contours at tenths
 * of its peak. The mode is signed (the diverging scale) and swings as
 * φ cos ωt, so its colours flip each half period while its nodal lines — the
 * zero contour, drawn in ink — stay still: they are where sand collects on a
 * bowed plate, Chladni's figures. The free plate (FFFF) is his own setting;
 * its three rigid modes, at ω = 0, are skipped, so "mode 1" is the first that
 * bends. On a square, modes come in degenerate pairs; any combination of a
 * pair is a mode, and the solver hands back one of them.
 *
 * Instead of a mode, the right panel can show the steady response to the load
 * applied harmonically at Ω (stage 8): Re(u e^{iΩt}) over the plate. Damping
 * puts its points out of phase with one another, so its zero contour is not
 * still — the "nodal lines" of a forced, damped plate travel.
 */

import type { Harmonic } from "../../hierarchy";
import { REFERENCE, exactPlateEigenvalues, navierResponse, navierStatic } from "../../plate/exact";
import { Plate, admissible, type PlateSpec } from "../../plate/plate";
import { evaluatePlateQoI, plateHarmonicOf, plateLoadOf, plateQoiPoint } from "../../plate/qoi";
import type { Modes } from "../../eig";
import type { Figure } from "../figure";
import { colourBar, contour, linspace, meshLines, paint, plateOutline, shapeLimits, type Grid } from "../heatmap";
import type { Axes, Plot, Series } from "../plot";
import { continuityK, continuityName, displayOf, plateCaseOf, type State } from "../state";
import { axisUnit, deflectionScale, fmt, plain, referenceNote, type Display } from "../units";
import { memo, table, type ViewResult } from "./view";

/** Elastic modes offered; the free plate's three rigid ones come on top. */
export const PLATE_MODES = 12;
const PERIOD_MS = 1600;
const INK = "rgba(207, 238, 255, 0.90)";
const RES = 129;

interface Solved {
  plate: Plate;
  c: Float64Array | null;
  F: Float64Array | null;
  why: string | null;
  modes: Modes;
  rigid: number;
  ms: number;
}

const solved = memo<Solved | string>();
const fields = memo<{ w: Grid | null; phi: Grid }>();
const responses = memo<{ h: Harmonic; re: Grid; im: Grid; at: number[] }>();
const MESH_INK = "rgba(5, 5, 12, 0.45)";

export function renderPlate(fig: Figure, st: State, t: number): ViewResult {
  const p = st.p, k = continuityK(st.continuity, p), pc = plateCaseOf(st);
  const spec: PlateSpec = { p, k, ne: st.plateNe, edges: st.edges, aspect: st.aspect, nu: st.nu };
  const key = JSON.stringify([p, k, st.plateNe, st.edges, st.aspect, st.nu, st.load]);
  const r = solved(key, () => {
    const why = admissible(spec, true);
    if (why) return why;
    const t0 = performance.now();
    const plate = new Plate(spec);
    const stat = admissible(spec);
    let c: Float64Array | null = null, F: Float64Array | null = null;
    if (!stat) {
      F = plate.load(plateLoadOf(pc));
      c = plate.solve(F);
    }
    const rigid = st.edges === "FFFF" ? 3 : 0;
    const modes = plate.modes(Math.min(PLATE_MODES + rigid, plate.dofs));
    return { plate, c, F, why: stat, modes, rigid, ms: performance.now() - t0 };
  });
  const [left, right] = fig.panels(2, { cols: 2 });
  if (typeof r === "string") {
    left.draw({ xlabel: "x / L", ylabel: "y / L", series: [] });
    right.draw({ xlabel: "x / L", ylabel: "y / L", series: [] });
    return { readout: `${continuityName(Math.max(k, 0))}, degree ${p}: not a plate discretisation.\n${r}`, animate: false };
  }

  const n = Math.min(st.mode, r.modes.values.length - r.rigid) - 1 + r.rigid;
  const a = st.aspect, nx = RES, ny = Math.max(9, Math.round(RES / a) | 1);
  const xs = linspace(0, a, nx), ys = linspace(0, 1, ny);
  const g = fields(`${key}/${n}`, () => ({
    w: r.c ? { nx, ny, x0: 0, x1: a, y0: 0, y1: 1, v: r.plate.grid(r.c, xs, ys) } : null,
    phi: { nx, ny, x0: 0, x1: a, y0: 0, y1: 1, v: r.plate.grid(r.modes.vectors[n], xs, ys) },
  }));

  const d = displayOf(st, "plate"), s = r.plate;
  const header = `${continuityName(k)} splines, degree ${p}, ${st.plateNe} × ${st.plateNe} elements on a ${+a.toFixed(3)} × 1 plate: ` +
    `${s.n ** 2} coefficients, ${s.dofs} free after the ${st.edges} edges; half-bandwidth ${s.K.bw}, ν = ${st.nu}`;
  let motion: string[];
  if (st.motion === "response") {
    const u = responses(`${key}/${st.forceRatio}/${st.zeta}`, () => {
      const h = plateHarmonicOf(pc), at = plateQoiPoint(pc);
      const { c, ci } = evaluatePlateQoI(r.plate, "response", plateLoadOf(pc), at, h);
      const grid = (v: Float64Array): Grid => ({ nx, ny, x0: 0, x1: a, y0: 0, y1: 1, v: r.plate.grid(v, xs, ys) });
      return { h, re: grid(c), im: grid(ci!), at: [r.plate.evaluate(c, at.x, at.y), r.plate.evaluate(ci!, at.x, at.y)] };
    });
    motion = response(right, st, d, r, u, t);
  } else motion = modes(right, st, d, r, g.phi, n, t);
  const lines = [header, "", ...statics(left, st, d, r, g.w), "", ...motion];
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return { readout: lines.join("\n"), animate: true };
}

/** Invisible points to hover, on a coarse sub-grid of the field: the series, and the node each point is. */
function hoverOf(g: Grid, L: number, stride = 8): { series: Series; node: (i: number) => number } {
  const nodes: number[] = [], x: number[] = [], y: number[] = [];
  for (let iy = 0; iy < g.ny; iy += stride)
    for (let ix = 0; ix < g.nx; ix += stride) {
      nodes.push(iy * g.nx + ix);
      x.push((g.x0 + ((g.x1 - g.x0) * ix) / (g.nx - 1)) * L);
      y.push((g.y0 + ((g.y1 - g.y0) * iy) / (g.ny - 1)) * L);
    }
  return { series: { label: "value", x, y, color: "rgba(0, 0, 0, 0)", width: 0.001, unlisted: true }, node: (i) => nodes[i] };
}

function statics(plot: Plot, st: State, d: Display, r: Solved, g: Grid | null): string[] {
  const L = d.dimensional ? d.ref.L : 1, a = st.aspect, at = plateQoiPoint(st);
  const lim = shapeLimits(plot, a * L, L);
  const axis = d.dimensional ? ["x [m]", "y [m]"] : ["x / L", "y / L"];
  if (!g || !r.c || !r.F) {
    plot.draw({
      title: "static deflection — none: a free plate has no static solution", xlabel: axis[0], ylabel: axis[1], ...lim, series: [],
      over: (ctx, ax) => plateOutline(ctx, scaled(ax, L), a, 1, st.edges),
    });
    return [`static: ${r.why}`];
  }
  const W = d.dimensional ? deflectionScale(d.ref, st.load) : 1;
  let top = 0;
  for (let i = 0; i < g.v.length; i++) top = Math.max(top, Math.abs(g.v[i]));
  const u = d.dimensional ? axisUnit(top * W, "m") : { factor: 1, label: "" };
  const k = W * u.factor;
  const hp = hoverOf(g, L);
  plot.draw({
    title: `static deflection — ${st.load === "uniform" ? "uniform pressure" : "point load"}, ${st.edges}`,
    xlabel: axis[0], ylabel: axis[1], ...lim, series: [hp.series],
    under: (ctx, ax) => {
      const sax = scaled(ax, L);
      paint(ctx, sax, g, "sequential", top);
      for (let j = 1; j < 10; j++) contour(ctx, sax, g, (j / 10) * top, "rgba(5, 5, 12, 0.55)", 1);
      if (st.mesh) meshLines(ctx, sax, r.plate.sx.breaks, r.plate.sy.breaks, MESH_INK);
    },
    over: (ctx, ax) => {
      const sax = scaled(ax, L);
      plateOutline(ctx, sax, a, 1, st.edges);
      marker(ctx, sax, at.x, at.y, st.load === "point");
      colourBar(ctx, ax, "sequential", top * k, d.dimensional ? `w [${u.label}]` : "w", (v) => plain(v, 3));
    },
    hover: (se, i) => `x = ${se.x[i].toPrecision(3)}, y = ${se.y[i].toPrecision(3)}\nw = ${(g.v[hp.node(i)] * k).toPrecision(4)}`,
  });

  const w = r.plate.evaluate(r.c, at.x, at.y);
  let C = 0;
  for (let j = r.plate.loy; j < r.plate.hiy; j++)
    for (let i = r.plate.lox; i < r.plate.hix; i++) C += r.F[r.plate.dof(i, j)] * r.c[j * r.plate.n + i];
  const ns = st.edges === "SSSS" ? navierStatic(a, st.load === "uniform" ? "uniform" : at) : null;
  const rel = (v: number, e: number | undefined) => (e === undefined ? "—" : Math.abs(v / e - 1).toExponential(2));
  const where = st.edges === "CFFF" ? "free-edge midpoint" : "centre";
  const lines = [
    table(["static", "spline", st.edges === "SSSS" ? "Navier" : "exact", "rel. error"], [
      [`w at the ${where}`, fmt.deflection(d, w, st.load), ns ? fmt.deflection(d, ns.w(at.x, at.y), st.load) : "—", rel(w, ns?.w(at.x, at.y))],
      ["compliance ℓ(w)", fmt.compliance(d, C, st.load), ns ? fmt.compliance(d, ns.compliance, st.load) : "—", rel(C, ns?.compliance)],
    ]),
    `assembled, solved, and ${r.modes.values.length} modes found in ${r.ms.toFixed(0)} ms`,
  ];
  if (!ns) lines.push(`${st.edges}: no closed form for the deflection — the convergence view measures it against itself.`);
  if (st.load === "point") lines.push("point load: w ~ r² log r under it, so its curvature is unbounded there and convergence is slow.");
  return lines;
}

function modes(plot: Plot, st: State, d: Display, r: Solved, g: Grid, n: number, t: number): string[] {
  const L = d.dimensional ? d.ref.L : 1, a = st.aspect, amp = Math.cos((2 * Math.PI * t) / PERIOD_MS);
  let top = 0;
  for (let i = 0; i < g.v.length; i++) top = Math.max(top, Math.abs(g.v[i]));
  const lim = shapeLimits(plot, a * L, L);
  const shown = n - r.rigid + 1;
  const hp = hoverOf(g, L);
  plot.draw({
    title: `mode ${shown} of the ${st.edges} plate, and its nodal lines`,
    xlabel: d.dimensional ? "x [m]" : "x / L", ylabel: d.dimensional ? "y [m]" : "y / L", ...lim, series: [hp.series],
    under: (ctx, ax) => {
      const sax = scaled(ax, L);
      paint(ctx, sax, g, "diverging", top, amp);
      if (st.mesh) meshLines(ctx, sax, r.plate.sx.breaks, r.plate.sy.breaks, MESH_INK);
      contour(ctx, sax, g, 0, INK, 2);
    },
    over: (ctx, ax) => {
      plateOutline(ctx, scaled(ax, L), a, 1, st.edges);
      colourBar(ctx, ax, "diverging", 1, "φ", (v) => plain(v, 2));
    },
    hover: (se, i) => `x = ${se.x[i].toPrecision(3)}, y = ${se.y[i].toPrecision(3)}\nφ = ${(g.v[hp.node(i)] / top).toPrecision(3)} (of its peak)`,
  });

  const count = r.modes.values.length - r.rigid;
  const exact = exactPlateEigenvalues(st.edges, a, count);
  const ref = a === 1 ? REFERENCE[st.edges] : undefined;
  const refOmega = (i: number) => (exact ? Math.sqrt(exact[i]) : ref && i + r.rigid < ref.omega.length ? ref.omega[i + r.rigid] : NaN);
  const unit = d.dimensional ? "f [Hz]" : "ω̂ = ω L² √(ρt/D)";
  const rows = Array.from({ length: count }, (_, i) => {
    const w = Math.sqrt(Math.max(0, r.modes.values[i + r.rigid])), we = refOmega(i);
    return [
      `${i + 1}${i + r.rigid === n ? " ◂" : ""}`,
      fmt.frequency(d, w),
      Number.isFinite(we) ? fmt.frequency(d, we) : "—",
      Number.isFinite(we) ? ((w - we) / we).toExponential(2) : "—",
    ];
  });
  const source = exact ? (st.edges === "SSSS" ? "Navier" : "Lévy") : ref ? ref.source : "reference";
  const lines = [`lowest ${count} natural frequencies, ${unit}`, table(["mode", "spline", source, "rel. error"], rows)];
  if (!exact && !ref) lines.push(`${st.edges}${a === 1 ? "" : `, aspect ${+a.toFixed(3)}`}: no closed form or table here — the convergence view measures it against itself.`);
  if (r.rigid) lines.push("free plate: three rigid modes at ω = 0 (two rotations, one translation) come first and are skipped; the solver shifts past them.");
  const lam = r.modes.values, gap = (i: number) => i >= r.rigid && i < lam.length && Math.abs(lam[i] / lam[n] - 1) < 1e-6;
  if (gap(n - 1) || gap(n + 1))
    lines.push("this mode is one of a degenerate pair: any combination of the two is a mode, so its nodal lines are one choice among many — as on a real square plate, where the pair's mix depends on how it is bowed.");
  lines.push("ink: the nodal lines, which stay still while the plate swings — where sand gathers in Chladni's figures.");
  return lines;
}

function response(plot: Plot, st: State, d: Display, r: Solved, u: { h: Harmonic; re: Grid; im: Grid; at: number[] }, t: number): string[] {
  const L = d.dimensional ? d.ref.L : 1, a = st.aspect, th = (2 * Math.PI * t) / PERIOD_MS, c = Math.cos(th), sn = Math.sin(th);
  let top = 0;
  for (let i = 0; i < u.re.v.length; i++) top = Math.max(top, Math.hypot(u.re.v[i], u.im.v[i]));
  const now: Grid = { ...u.re, v: Float64Array.from(u.re.v, (v, i) => v * c - u.im.v[i] * sn) };
  const W = d.dimensional ? deflectionScale(d.ref, st.load) : 1;
  const un = d.dimensional ? axisUnit(top * W, "m") : { factor: 1, label: "" };
  const k = W * un.factor, hp = hoverOf(u.re, L);
  plot.draw({
    title: `forced response at Ω = ${st.forceRatio} ω₁, ζ = ${st.zeta} — Re(u e^{iΩt}), and its moving zero line`,
    xlabel: d.dimensional ? "x [m]" : "x / L", ylabel: d.dimensional ? "y [m]" : "y / L", ...shapeLimits(plot, a * L, L),
    series: [hp.series],
    under: (ctx, ax) => {
      const sax = scaled(ax, L);
      paint(ctx, sax, now, "diverging", top);
      if (st.mesh) meshLines(ctx, sax, r.plate.sx.breaks, r.plate.sy.breaks, MESH_INK);
      contour(ctx, sax, now, 0, INK, 2);
    },
    over: (ctx, ax) => {
      const sax = scaled(ax, L), at = plateQoiPoint(st);
      plateOutline(ctx, sax, a, 1, st.edges);
      marker(ctx, sax, at.x, at.y, st.load === "point");
      colourBar(ctx, ax, "diverging", top * k, d.dimensional ? `w [${un.label}]` : "w", (v) => plain(v, 3));
    },
    hover: (se, i) => {
      const j = hp.node(i);
      return `x = ${se.x[i].toPrecision(3)}, y = ${se.y[i].toPrecision(3)}\n|u| = ${(Math.hypot(u.re.v[j], u.im.v[j]) * k).toPrecision(4)}`;
    },
  });

  const at = plateQoiPoint(st), aq = Math.hypot(u.at[0], u.at[1]);
  const exact = st.edges === "SSSS" ? navierResponse(a, st.load, at, u.h) : null;
  const wq = r.c ? r.plate.evaluate(r.c, at.x, at.y) : NaN;
  // As a lag in [0°, 360°).
  const lag = ((((Math.atan2(-u.at[1], u.at[0]) * 180) / Math.PI) % 360) + 360) % 360;
  const where = st.edges === "CFFF" ? "free-edge midpoint" : "centre";
  return [
    `forced at Ω = ${plain(u.h.Omega, 5)} (${st.forceRatio} × ω₁ of the uniform plate), Rayleigh damping C = aM + bK with ` +
      `a = ${plain(u.h.a, 4)}, b = ${plain(u.h.b, 4)}: ζ = ${st.zeta} on its first two ${st.edges === "FFFF" ? "elastic " : ""}modes`,
    table([`at the ${where}`, "spline", exact === null ? "exact" : "Navier series", "rel. error"], [
      ["amplitude |u|", fmt.deflection(d, aq, st.load), exact === null ? "—" : fmt.deflection(d, exact, st.load),
        exact === null ? "—" : Math.abs(aq / exact - 1).toExponential(2)],
      ["static w", Number.isFinite(wq) ? fmt.deflection(d, wq, st.load) : "— (a free plate has none)", "", ""],
    ]),
    Number.isFinite(wq)
      ? `dynamic amplification |u| / w_static = ${plain(aq / Math.abs(wq), 4)}; the response lags the load by ${plain(lag, 3)}°`
      : `the response lags the load by ${plain(lag, 3)}° at the ${where}`,
    "ink: Re(u e^{iΩt}) = 0. With damping the points of the plate are out of phase with each other, so these lines travel — a mode's would stand still.",
  ];
}

/** Axes in plate units when the plot is in metres: everything drawn in plate coordinates scales by L. */
function scaled(a: Axes, L: number): Axes {
  return { ...a, sx: (x: number) => a.sx(x * L), sy: (y: number) => a.sy(y * L) };
}

/** Where the deflection is read (and the point load acts). */
function marker(ctx: CanvasRenderingContext2D, a: { sx: (x: number) => number; sy: (y: number) => number }, x: number, y: number, load: boolean): void {
  const X = a.sx(x), Y = a.sy(y);
  ctx.strokeStyle = load ? "#c98500" : INK;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(X, Y, 5, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.beginPath(); ctx.moveTo(X - 9, Y); ctx.lineTo(X + 9, Y); ctx.moveTo(X, Y - 9); ctx.lineTo(X, Y + 9); ctx.stroke();
}
