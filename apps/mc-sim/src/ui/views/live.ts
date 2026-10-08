// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 8 on screen: multilevel Monte Carlo at work, the solver beside its
 * statistics — live, redrawn as every batch of samples arrives.
 *
 * Across the top, one sample of the run as the workers solved it: the random
 * stiffness it drew, its solution on level ℓ and on level ℓ − 1 from the same
 * ω — the coupled pair whose difference Y_ℓ the estimator averages — and the
 * mesh of every level of the hierarchy, the pair's two picked out. The view
 * steps through the level's samples every few seconds; whatever the quantity
 * is, the structure moves as it does — the first mode swinging for ω₁, the
 * forced response for the response, still for a static quantity. Each pair is
 * re-solved here from its index alone and checked against the number the
 * workers returned: they agree to the last bit, which is what "a sample is a
 * pure function of (seed, level, i)" means.
 *
 * Below, the run's statistics from every sample each level has so far — not
 * the survey's fixed prefix the MLMC view reads, so they move as the run does:
 * the distribution of Q₀ with the shown sample and the running MLMC estimate
 * (the telescoping sum over the levels the chosen tolerance uses, ± 95% of its
 * sampling error); the variance of Q_ℓ and Y_ℓ against level; and the samples
 * in hand on each level, rising toward what the chosen tolerance asks for.
 * "run again from zero" (pane) throws the run away so it can be watched
 * arriving: the same samples again, since a run is a function of its seed.
 *
 * "compare cost with MC" swaps those three for the cost comparison against
 * plain Monte Carlo at the chosen tolerance (`cost.ts`): samples per level for
 * each method, the total costs drawn to scale, samples against ε, and ε² × cost
 * against ε — and adds its table to the readout. It is the MLMC view's run — the
 * same samples, kept in the pool — so switching between the two loses nothing.
 */

import { THETA, slope } from "../../mc/mlmc";
import { Sampler, type Inspected } from "../../mc/sampler";
import { Z95, histogram } from "../../mc/stats";
import { FieldAt, draw } from "../../random/field";
import { FieldOnGrid } from "../../random/field2d";
import type { Figure } from "../figure";
import { colourBar, linspace, meshLines, meshTicks, paint, plateOutline, shapeLimits, type Grid, type Scale } from "../heatmap";
import { SLOT, type Axes, type Plot, type Series } from "../plot";
import { continuityK, continuityName, displayOf, type State } from "../state";
import { plain, referenceNote } from "../units";
import { QOI_NAME, valueOf } from "./convergence";
import { KERNEL_NAME } from "./field";
import {
  MLMC_INK, bar, costLines, drawPerLevel, drawSamplesVsEps, drawTotal, label, pairOf, type Pair as CostPair,
} from "./cost";
import { SWEEP, drawCost, fmt2, mlmcSession, type MlmcSession } from "./mlmc";
import { fmtCount, runStatus, vline } from "./montecarlo";
import { isPlate, structureText, workUnit } from "./structure";
import { memo, table, type ViewResult } from "./view";
import { workers } from "./workers";

/** Each pair stays on screen this long. */
const DWELL_MS = 2500;
/** The view cycles through this many of a level's first samples. */
const CYCLE = 48;
const PERIOD_MS = 1600;
const INK = "rgba(207, 238, 255, 0.85)";
const MUTED = "rgba(207, 238, 255, 0.35)";
const RES = 97;

const samplers = memo<Sampler>();

/** A shape that moves as a cos θ − b sin θ (b absent: a cos θ; still: a alone). */
interface Motion {
  readonly a: Float64Array;
  readonly b: Float64Array | null;
  readonly still: boolean;
}

function motionOf(ins: Inspected): Motion {
  if (ins.mode) return { a: ins.mode, b: null, still: false };
  if (ins.u) return { a: ins.u.re, b: ins.u.im, still: false };
  return { a: ins.w, b: null, still: true };
}

interface Pair {
  readonly i: number;
  readonly fine: Inspected;
  readonly coarse: Inspected;
  /** The shapes, sampled for drawing: along the beam, or on the plate's grid. */
  readonly fa: Float64Array; readonly fb: Float64Array | null;
  readonly ca: Float64Array; readonly cb: Float64Array | null;
  readonly still: boolean;
  /** log(e/e₀) of the sample at the drawing points. */
  readonly loge: Float64Array;
  readonly top: number;
}

const pairs = memo<Pair>();

export function renderLive(fig: Figure, st: State, t: number): ViewResult {
  const c = mlmcSession(st);
  if (typeof c === "string") {
    fig.panels(1)[0].draw({ xlabel: "", ylabel: "", series: [] });
    return { readout: c, animate: false };
  }
  const panels = st.liveCompare ? fig.panels(5, { cols: 2, wideFirst: true }) : fig.panels(4, { cols: 3, wideFirst: true });
  const [solver, ...below] = panels;
  const lv = Math.max(1, Math.min(st.liveLevel, c.levels - 1));
  const acc = c.run.streams[lv].acc;
  const live = liveLevels(c), tol = costAt(c);
  if (acc.n < 1) {
    solver.draw({ xlabel: "", ylabel: "", series: [], title: `waiting for the first samples of level ${lv}…` });
    for (const pl of below) pl.draw({ xlabel: "", ylabel: "", series: [] });
    return { readout: summary(c, lv, null, tol, live), animate: true };
  }

  const i = Math.floor(t / DWELL_MS) % Math.min(acc.n, CYCLE);
  const sampler = samplers(JSON.stringify([c.spec, c.M]), () => new Sampler(c.spec, c.kl, c.klY));
  const pair = pairs(JSON.stringify([c.spec, c.M, c.neOf(0), lv, i]), () => solvePair(c, sampler, lv, i));
  const stored = { Q: acc.values[i], Qc: acc.coarseValues[i] };
  const th = (2 * Math.PI * t) / PERIOD_MS;

  if (isPlate(st)) drawPlates(solver, c, lv, pair, th);
  else drawBeams(solver, c, lv, pair, th);
  const waiting = (pl: Plot, what: string) => pl.draw({ xlabel: "", ylabel: "", series: [], title: what });
  const settling = `waiting for ε = ${(SWEEP[tolIndex(st)] * st.mlEps).toExponential(1)} to settle on its levels…`;
  if (st.liveCompare) {
    const [perLevel, total, samples, cost] = below;
    if (tol) {
      drawPerLevel(perLevel, tol, {
        highlight: lv,
        title: `samples per level, ε = ${tol.eps.toExponential(1)}: plain MC ${plain(tol.mcCost / tol.mlCost, 3)}× the cost${tol.converged ? "" : " (settling)"}`,
      });
      drawTotal(total, c, tol);
    } else for (const pl of [perLevel, total]) waiting(pl, settling);
    drawSamplesVsEps(samples, c, SWEEP.map((_, k) => pairOf(c, k)));
    drawCost(cost, c);
  } else {
    const [hist, variance, progress] = below;
    drawHistogram(hist, c, pair, live, tol);
    drawLiveVariance(variance, c, live, lv);
    drawProgress(progress, c, live, lv, tol);
  }
  return { readout: summary(c, lv, { pair, stored }, tol, live), animate: true };
}

const tolIndex = (st: State) => Math.max(0, Math.min(SWEEP.length - 1, st.cmpTol));

/** The cost comparison at the tolerance the pane picks. */
const costAt = (c: MlmcSession): CostPair | null => pairOf(c, tolIndex(c.st));

/** What a level's stream says now, from every sample it has — the live counterpart of the survey's fixed prefix. */
interface LiveLevel {
  readonly n: number;
  readonly meanQ: number;
  readonly varQ: number;
  readonly meanY: number;
  readonly varY: number;
}

function liveLevels(c: MlmcSession): LiveLevel[] {
  return c.run.streams.slice(0, c.levels).map(({ acc }) =>
    ({ n: acc.n, meanQ: acc.q.mean, varQ: acc.q.variance, meanY: acc.y.mean, varY: acc.y.variance }));
}

/** |E[Q]| to make things relative to: the survey's, once it is in; level 0's running mean before. */
const scaleOf = (c: MlmcSession, live: readonly LiveLevel[]) =>
  Number.isFinite(c.scale) ? c.scale : Math.abs(live[0]?.meanY) || 1;

/**
 * The running telescoping sum Σ_{ℓ ≤ L} mean(Y_ℓ) over every sample in hand,
 * and its standard error √(Σ V_ℓ/n_ℓ), on the levels the chosen tolerance uses
 * — or every level with two samples, before it has settled. The adaptive
 * algorithm's own estimate reads prefixes and changes only when a round
 * completes; this one moves with every batch.
 */
function liveEstimate(live: readonly LiveLevel[], tol: CostPair | null): { L: number; est: number; se: number } | null {
  let L = tol ? Math.min(tol.L, live.length - 1) : live.length - 1;
  while (L >= 0 && !(live[L].n >= 2)) L--;
  if (L < 0) return null;
  let est = 0, v = 0;
  for (let l = 0; l <= L; l++) { est += live[l].meanY; v += live[l].varY / live[l].n; }
  return { L, est, se: Math.sqrt(v) };
}

function drawLiveVariance(plot: Plot, c: MlmcSession, live: readonly LiveLevel[], lv: number): void {
  const s2 = scaleOf(c, live) ** 2, ok = live.map((v) => v.n >= 2);
  const l = live.map((_, k) => k);
  const series: Series[] = [
    { label: "V[Q_ℓ]", x: l, y: live.map((v, k) => (ok[k] ? v.varQ / s2 : NaN)), color: SLOT[0], width: 2, markers: true },
    { label: "V[Y_ℓ] = V[Q_ℓ − Q_ℓ₋₁]", x: l, y: live.map((v, k) => (ok[k] && k > 0 ? v.varY / s2 : NaN)), color: SLOT[1], width: 2, markers: true },
  ];
  const ys = series.flatMap((s) => Array.from(s.y)).filter((v) => v > 0 && Number.isFinite(v));
  plot.draw({
    title: `variance against level, from every sample so far`,
    xlabel: "level ℓ", ylabel: "relative to E[Q]²", ylog: true, legend: "bl",
    xlim: [-0.3, c.levels - 0.7], ylim: ys.length ? [Math.min(...ys) / 4, Math.max(...ys) * 4] : [1e-6, 1],
    xfmt: (v) => (Number.isInteger(Math.round(v * 1e6) / 1e6) ? String(Math.round(v)) : ""),
    series,
    over: (ctx, a) => {
      // The level whose pair is on screen.
      ctx.strokeStyle = "rgba(207, 238, 255, 0.35)";
      ctx.setLineDash([3, 4]);
      ctx.beginPath(); ctx.moveTo(a.sx(lv), a.box.t); ctx.lineTo(a.sx(lv), a.box.b); ctx.stroke();
      ctx.setLineDash([]);
    },
    hover: (s, k) => `${s.label}, level ${k}\n${s.y[k].toExponential(3)} (relative)\nfrom ${live[k].n.toLocaleString()} samples`,
  });
}

/**
 * The samples in hand on each level, rising as batches arrive, against what
 * the chosen tolerance asks of each (ticks) and the survey's share of every
 * level (the dashed line). The shown level's bar is outlined.
 */
function drawProgress(plot: Plot, c: MlmcSession, live: readonly LiveLevel[], lv: number, tol: CostPair | null): void {
  const n = live.map((v) => v.n), survey = c.st.mlSurvey, w = 0.3;
  const want = tol ? tol.N : [];
  const ys = [...n, ...want, survey].filter((v) => v > 0);
  const top = Math.max(...ys) * 30, lo = Math.max(1, Math.min(...ys) / 4);
  const total = n.reduce((s, v) => s + v, 0);
  plot.draw({
    title: `samples in hand: ${total.toLocaleString()}${tol ? `, toward ε = ${tol.eps.toExponential(1)}` : ""}`,
    xlabel: "level ℓ", ylabel: "samples", ylog: true, legend: "tr",
    xlim: [-0.6, c.levels - 0.4], ylim: [lo, top],
    xfmt: (v) => (Number.isInteger(Math.round(v * 1e6) / 1e6) && v >= 0 ? String(Math.round(v)) : ""),
    series: [
      { label: "in hand", x: [], y: [], color: MLMC_INK, width: 8, inert: true },
      ...(tol ? [{ label: `asked for at ε = ${tol.eps.toExponential(1)}`, x: [], y: [], color: INK, width: 2, inert: true }] : []),
      { label: `survey, ${survey.toLocaleString()} per level`, x: [], y: [], color: MUTED, width: 1.5, dash: [4, 4], inert: true },
      { label: "level", x: n.map((_, l) => l), y: n.map((v) => Math.max(v, lo)), color: "rgba(0, 0, 0, 0)", width: 0.001, unlisted: true },
    ],
    under: (ctx, a) => {
      n.forEach((v, l) => { if (v > 0) bar(ctx, a, l, w, v, MLMC_INK); });
      ctx.strokeStyle = MUTED;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      ctx.beginPath(); ctx.moveTo(a.box.l, a.sy(survey)); ctx.lineTo(a.box.r, a.sy(survey)); ctx.stroke();
      ctx.setLineDash([]);
    },
    over: (ctx, a) => {
      want.forEach((v, l) => {
        ctx.strokeStyle = INK;
        ctx.lineWidth = 2;
        ctx.beginPath(); ctx.moveTo(a.sx(l - w - 0.08), a.sy(v)); ctx.lineTo(a.sx(l + w + 0.08), a.sy(v)); ctx.stroke();
      });
      n.forEach((v, l) => { if (v > 0) label(ctx, a.sx(l), Math.max(a.box.t + 12, a.sy(Math.max(v, want[l] ?? 0)) - 4), fmtCount(v)); });
      if (n[lv] > 0) {
        const X0 = a.sx(lv - w) - 3, X1 = a.sx(lv + w) + 3, Y = Math.max(a.box.t, a.sy(n[lv])) - 3;
        ctx.strokeStyle = INK;
        ctx.lineWidth = 2;
        ctx.beginPath(); ctx.roundRect(X0, Y, X1 - X0, a.box.b - Y, [6, 6, 0, 0]); ctx.stroke();
      }
    },
    hover: (_, l) => `level ${l}: ${n[l].toLocaleString()} samples in hand` +
      (want[l] !== undefined ? `\n${want[l].toLocaleString()} asked for at ε = ${tol!.eps.toExponential(1)}` : "") +
      (l === lv ? "\nthe level shown above" : ""),
  });
}

function solvePair(c: MlmcSession, sampler: Sampler, lv: number, i: number): Pair {
  const fine = sampler.inspect(i, c.neOf(lv), lv), coarse = sampler.inspect(i, c.neOf(lv - 1), lv);
  const mf = motionOf(fine), mc = motionOf(coarse), st = c.st;
  const normals = draw(c.spec.field, sampler.M, c.spec.seed, i, lv).xi;
  let fa: Float64Array, fb: Float64Array | null, ca: Float64Array, cb: Float64Array | null, loge: Float64Array;
  if (fine.plate) {
    const a = st.aspect, xs = linspace(0, a, RES), ys = linspace(0, 1, Math.max(9, Math.round(RES / a) | 1));
    const g = (ins: Inspected, v: Float64Array | null) => (v ? (ins.plate!).grid(v, xs, ys) : null);
    fa = g(fine, mf.a)!; fb = g(fine, mf.b); ca = g(coarse, mc.a)!; cb = g(coarse, mc.b);
    loge = new FieldOnGrid(c.kl, c.klY ?? c.kl, st.terms, xs, ys, a).lognormal(st.sigma, normals).map(Math.log);
  } else {
    const xs = linspace(0, 1, 361);
    const ev = (ins: Inspected, v: Float64Array | null) => (v ? xs.map((x) => ins.beam!.evaluate(v, x)[0]) : null);
    fa = ev(fine, mf.a)!; fb = ev(fine, mf.b); ca = ev(coarse, mc.a)!; cb = ev(coarse, mc.b);
    loge = new FieldAt(c.kl, xs).lognormal(st.sigma, normals).map(Math.log);
  }
  // A fixed scale for the pair's whole swing: the larger envelope of the two.
  let top = 0;
  for (const [A, B] of [[fa, fb], [ca, cb]] as const)
    for (let k = 0; k < A.length; k++) top = Math.max(top, Math.hypot(A[k], B ? B[k] : 0));
  return { i, fine, coarse, fa, fb, ca, cb, still: mf.still, loge, top: top || 1 };
}

/** The shape at phase θ. */
const at = (A: Float64Array, B: Float64Array | null, still: boolean, th: number): Float64Array =>
  still ? A : A.map((v, k) => v * Math.cos(th) - (B ? B[k] * Math.sin(th) : 0));

function what(st: State): string {
  return st.qoi === "omega1" ? "first mode φ₁, swinging" : st.qoi === "response" ? "forced response Re(u e^{iΩt})" : "static deflection";
}

function drawBeams(plot: Plot, c: MlmcSession, lv: number, p: Pair, th: number): void {
  const st = c.st, d = displayOf(st), L = d.dimensional ? d.ref.L : 1;
  const xs = linspace(0, L, p.fa.length), top = p.top;
  const signed = !p.still || p.fa.some((v) => v < -1e-9 * top);
  const lo = signed ? -1.15 * top : -0.2 * top, hi = 1.15 * top;
  // Below the curves, in screen space: the field as a strip, then one row of mesh ticks per level
  // (`beamDecor`). Room for them is reserved in pixels — the panel is short when the cost comparison
  // shares the figure — as the fraction of the plot box they need.
  const boxH = Math.max(60, plot.canvas.clientHeight - 30 - 46);
  const need = Math.min(0.72, (14 + 12 + 14 + 12 * c.levels + 6) / boxH);
  const extra = (need * (hi - lo)) / (1 - need);
  plot.draw({
    title: `sample ${p.i} of level ${lv}: ${what(st)} on ${c.neOf(lv)} and on ${c.neOf(lv - 1)} elements, same ω`,
    xlabel: d.dimensional ? "x [m]" : "x / L", ylabel: st.qoi === "omega1" ? "φ (M-normalised)" : "w, nondimensional  (downward)",
    xlim: [0, L], ylim: [lo, hi + extra], yflip: true, legend: "tr",
    yfmt: (v) => (v > hi ? "" : String(+v.toPrecision(3))),
    series: [
      { label: `level ${lv}: ${c.neOf(lv)} elements`, x: xs, y: at(p.fa, p.fb, p.still, th), color: SLOT[0], width: 2.5 },
      { label: `level ${lv - 1}: ${c.neOf(lv - 1)} elements`, x: xs, y: at(p.ca, p.cb, p.still, th), color: SLOT[1], width: 1.75, dash: [6, 4] },
    ],
    over: (ctx, a) => beamDecor(ctx, a, c, lv, p, L, hi),
  });
}

function beamDecor(ctx: CanvasRenderingContext2D, a: Axes, c: MlmcSession, lv: number, p: Pair, L: number, hi: number): void {
  const y0 = a.sy(hi) + 14, sh = 12, sigma = c.st.sigma || 1;
  const strip: Axes = { ...a, sy: (y) => y0 + sh - y * sh };
  const g: Grid = { nx: p.loge.length, ny: 2, x0: 0, x1: L, y0: 0, y1: 1, v: Float64Array.from({ length: 2 * p.loge.length }, (_, k) => p.loge[k % p.loge.length]) };
  paint(ctx, strip, g, "diverging", 3 * sigma);
  ctx.fillStyle = INK;
  ctx.font = "11px ui-monospace, monospace";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  ctx.fillText("log e", a.box.l - 6, y0 + sh / 2);
  for (let l = 0; l < c.levels; l++) {
    const Y = y0 + sh + 14 + 12 * l, ne = c.neOf(l);
    const color = l === lv ? SLOT[0] : l === lv - 1 ? SLOT[1] : MUTED;
    const row: Axes = { ...a, sy: () => Y };
    meshTicks(ctx, row, Array.from({ length: ne + 1 }, (_, k) => (L * k) / ne), 0, color, 4);
    ctx.fillStyle = color;
    ctx.fillText(`ℓ=${l}`, a.box.l - 6, Y);
  }
}

function drawPlates(plot: Plot, c: MlmcSession, lv: number, p: Pair, th: number): void {
  const st = c.st, a = st.aspect, gap = 0.12 * a, W = 3 * a + 2 * gap;
  const ny = p.loge.length / RES;
  const grid = (x0: number, v: Float64Array): Grid => ({ nx: RES, ny, x0, x1: x0 + a, y0: 0, y1: 1, v });
  const fine = grid(a + gap, at(p.fa, p.fb, p.still, th)), coarse = grid(2 * (a + gap), at(p.ca, p.cb, p.still, th));
  const field = grid(0, p.loge);
  const scale: Scale = p.still && !p.fa.some((v) => v < -1e-9 * p.top) ? "sequential" : "diverging";
  const sigma = st.sigma || 1;
  const pf = p.fine.plate!, pc = p.coarse.plate!;
  plot.draw({
    title: `sample ${p.i} of level ${lv}: its stiffness, and its ${what(st)} on ${c.neOf(lv)}² and ${c.neOf(lv - 1)}² elements, same ω`,
    xlabel: "", ylabel: "y / L", ...shapeLimits(plot, W, 1), xfmt: () => "", series: [],
    under: (ctx, ax) => {
      paint(ctx, ax, field, "diverging", 3 * sigma);
      paint(ctx, ax, fine, scale, p.top);
      paint(ctx, ax, coarse, scale, p.top);
      const shift = (x0: number): Axes => ({ ...ax, sx: (x: number) => ax.sx(x + x0) });
      meshLines(ctx, shift(a + gap), pf.sx.breaks, pf.sy.breaks, "rgba(5, 5, 12, 0.5)");
      meshLines(ctx, shift(2 * (a + gap)), pc.sx.breaks, pc.sy.breaks, "rgba(5, 5, 12, 0.5)");
    },
    over: (ctx, ax) => {
      const shift = (x0: number): Axes => ({ ...ax, sx: (x: number) => ax.sx(x + x0) });
      plateOutline(ctx, shift(0), a, 1, "FFFF");
      plateOutline(ctx, shift(a + gap), a, 1, st.edges);
      plateOutline(ctx, shift(2 * (a + gap)), a, 1, st.edges);
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      const label = (x0: number, text: string, color: string) => { ctx.fillStyle = color; ctx.fillText(text, ax.sx(x0 + a / 2), ax.sy(0) + 6); };
      label(0, `log(D/D₀) of sample ${p.i}`, INK);
      label(a + gap, `level ${lv}: Q = ${plain(p.fine.Q, 6)}`, SLOT[0]);
      label(2 * (a + gap), `level ${lv - 1}: Q = ${plain(p.coarse.Q, 6)}`, SLOT[1]);
      colourBar(ctx, ax, scale, 1, st.qoi === "omega1" ? "φ / max" : "w / max", (v) => plain(v, 2));
    },
  });
}

function drawHistogram(plot: Plot, c: MlmcSession, p: Pair, live: readonly LiveLevel[], tol: CostPair | null): void {
  const acc = c.run.streams[0].acc, { mean, sd } = acc.q;
  if (acc.n < 2) { plot.draw({ xlabel: "Q", ylabel: "", series: [] }); return; }
  const h = histogram(acc.values, acc.min, acc.max, sd), edges = h.edges, dens = h.density;
  const now = liveEstimate(live, tol), best = now !== null;
  const est = now ? now.est : NaN, half = now ? Z95 * now.se : NaN;
  const pad = 0.04 * (edges[edges.length - 1] - edges[0] || 1);
  const lo = Math.min(edges[0], p.fine.Q, p.coarse.Q) - pad, hi = Math.max(edges[edges.length - 1], p.fine.Q, p.coarse.Q) + pad;
  const topD = Math.max(...dens) * 1.2;
  const series: Series[] = [
    { label: `Q₀, ${acc.n.toLocaleString()} samples`, x: [], y: [], color: "rgba(57, 135, 229, 0.55)", width: 8, inert: true },
    { label: "bin", x: Array.from(dens, (_, k) => 0.5 * (edges[k] + edges[k + 1])), y: Array.from(dens), color: "rgba(0, 0, 0, 0)", width: 0.001, unlisted: true },
    { label: "E[Q₀]", x: [], y: [], color: MUTED, width: 2, inert: true },
    { label: "the shown sample, Q_ℓ", x: [], y: [], color: SLOT[0], width: 2, inert: true },
  ];
  if (best) series.push({ label: `MLMC estimate now, levels 0 … ${now!.L} ± 95%`, x: [], y: [], color: INK, width: 2, inert: true });
  plot.draw({
    title: "Q₀: its histogram, the sample, the estimate",
    xlabel: "Q", ylabel: "probability density", xlim: [lo, hi], ylim: [0, topD], legend: "tl", series,
    under: (ctx, a) => {
      ctx.fillStyle = "rgba(57, 135, 229, 0.55)";
      for (let k = 0; k < dens.length; k++) {
        const x0 = a.sx(edges[k]), x1 = a.sx(edges[k + 1]), y = a.sy(dens[k]);
        ctx.fillRect(x0 + 0.5, y, Math.max(1, x1 - x0 - 1), a.sy(0) - y);
      }
      if (best) {
        ctx.fillStyle = "rgba(207, 238, 255, 0.12)";
        ctx.fillRect(a.sx(est - half), a.box.t, a.sx(est + half) - a.sx(est - half), a.box.b - a.box.t);
      }
    },
    over: (ctx, a) => {
      vline(ctx, a, mean, MUTED, []);
      vline(ctx, a, p.fine.Q, SLOT[0], []);
      if (best) vline(ctx, a, est, INK, [5, 4]);
    },
    hover: (_, k) => `bin [${edges[k].toPrecision(4)}, ${edges[k + 1].toPrecision(4)}]\n${Math.round(dens[k] * (edges[k + 1] - edges[k]) * acc.n)} samples`,
  });
}

function summary(
  c: MlmcSession, lv: number, shown: { pair: Pair; stored: { Q: number; Qc: number } } | null, tol: CostPair | null,
  live: readonly LiveLevel[],
): string {
  const st = c.st, d = displayOf(st), v = (x: number) => valueOf(d, st.qoi, x, st.load);
  const elems = (n: number) => (isPlate(st) ? `${n} × ${n}` : String(n));
  const lines = [
    `${structureText(st)}, p = ${st.p} ${continuityName(continuityK(st.continuity, st.p))}; Q = ${QOI_NAME[st.qoi]}; ` +
      `${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L, σ = ${plain(st.sigma, 3)}, M = ${c.M}; seed ${st.seed}`,
  ];
  const n = c.run.streams.map((s) => s.acc.n);
  lines.push(`${runStatus(c.run)}: ${n.reduce((x, y) => x + y, 0).toLocaleString()} samples (${n.map(fmtCount).join(" · ")} by level) on ${workers.size} workers — the MLMC view's run`);
  if (shown) {
    const { pair: p, stored } = shown, Y = p.fine.Q - p.coarse.Q;
    const same = p.fine.Q === stored.Q && p.coarse.Q === stored.Qc;
    lines.push(
      "",
      `sample ${p.i} of level ${lv} (random stream ${lv}): ${elems(c.neOf(lv))} elements and its parent ${elems(c.neOf(lv - 1))}, the same ω — ` +
        `the view steps to the next every ${DWELL_MS / 1000} s, through the first ${CYCLE}`,
      table(["", "Q", "|Y| / |Q|"], [
        [`level ${lv}`, v(p.fine.Q), ""],
        [`level ${lv - 1}`, v(p.coarse.Q), ""],
        ["Y_ℓ = Q_ℓ − Q_ℓ₋₁", v(Y), (Math.abs(Y) / Math.abs(p.fine.Q)).toExponential(2)],
      ]),
      same
        ? "re-solved here from (seed, level, i) alone: the same two numbers the workers returned, to the last bit."
        : `re-solved here: differs from the workers' by ${(p.fine.Q - stored.Q).toExponential(2)}, ${(p.coarse.Q - stored.Qc).toExponential(2)} — the sample is not a pure function of its index`,
      `${isPlate(st) ? "the left plate" : "the strip under the beam"} is this ω, log(e/e₀); both meshes read it at their own quadrature points — ` +
        "the same function of x on both — which is why the two solutions are so close, and V[Y] so much smaller than V[Q].",
    );
  }
  if (tol) {
    const sumN = tol.N.reduce((s, n) => s + n, 0);
    lines.push("", lv <= tol.L
      ? `at ε = ${tol.eps.toExponential(1)} this pair is one of MLMC's ${fmtCount(tol.N[lv])} samples on level ${lv}, each costing ${fmtCount(tol.C[lv])} ` +
        `(${workUnit(st)}; fine and coarse solve) — ${fmtCount(sumN)} samples in all across levels 0 … ${tol.L}. ` +
        `Plain MC would need ${fmtCount(tol.nMC)}, every one on level ${tol.L} at ${fmtCount(tol.cMC)}: ${plain(tol.mcCost / tol.mlCost, 3)}× MLMC's cost${tol.note}.`
      : `at ε = ${tol.eps.toExponential(1)} MLMC stops at level ${tol.L}: level ${lv}'s samples are the survey's alone. ` +
        `Plain MC would need ${fmtCount(tol.nMC)} on level ${tol.L}: ${plain(tol.mcCost / tol.mlCost, 3)}× MLMC's cost${tol.note}.`);
    if (!st.liveCompare) lines.push(`"compare cost with MC" (pane) draws the comparison to scale.`);
  }
  const have = live.filter((v) => v.n >= 2);
  if (have.length > 1) {
    const alpha = -slope(live.map((v) => (v.n >= 2 ? Math.abs(v.meanY) : NaN)));
    const beta = -slope(live.map((v) => (v.n >= 2 ? v.varY : NaN)));
    const gamma = slope(live.map((_, l) => c.cost(l)));
    const now = liveEstimate(live, tol), scale = scaleOf(c, live);
    lines.push("", `live, from every sample so far: α ≈ ${fmt2(alpha)}, β ≈ ${fmt2(beta)}, γ = ${fmt2(gamma)}; ` +
      `V[Y_ℓ]/V[Q_ℓ] = ${live.slice(1).map((v) => (v.n >= 2 ? (v.varY / v.varQ).toExponential(1) : "—")).join(", ")}`);
    if (now) lines.push(`running MLMC estimate over levels 0 … ${now.L}: ${v(now.est)} ± ${(Z95 * now.se / scale).toExponential(1)} (95%, relative)`);
  }
  if (st.liveCompare) lines.push("", ...costLines(c, tolIndex(st)));
  const runs = c.sweep?.runs ?? [];
  if (runs.length) {
    lines.push("", `tolerances (θ = ${THETA}): ` + runs.map((r, k) => `ε = ${(SWEEP[k] * st.mlEps).toExponential(1)} ${r.status}` +
      (r.status === "converged" ? ` (L = ${r.L}, ${v(r.estimate)} ± ${(Z95 * Math.sqrt(r.sampleVariance) / c.scale).toExponential(1)})` : "")).join("; "));
  }
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}
