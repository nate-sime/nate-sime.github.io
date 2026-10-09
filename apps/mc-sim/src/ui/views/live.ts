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
 * Below, three of the MLMC view's panels, from every sample each level has so
 * far — not the survey's fixed prefix — so they move as the run does: log₂
 * variance and log₂ |mean| of Q_ℓ (dashed) and Y_ℓ (solid) against level, and
 * N_ℓ: the samples in hand on each level against what the chosen tolerance
 * asks for. A dotted line marks the level whose pair is shown. "run again from
 * zero" (pane) throws the run away so it can be watched arriving: the same
 * samples again, since a run is a function of its seed. It is the MLMC view's
 * run — the same samples, kept in the pool — so switching loses nothing.
 */

import { THETA, slope } from "../../mc/mlmc";
import { Sampler, type Inspected } from "../../mc/sampler";
import { Z95 } from "../../mc/stats";
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
  SWEEP, drawMean, drawVariance, fmt2, levelAxis, mcCost, mlmcSession, multi, single, tolIndex, type LevelStats, type MlmcSession,
} from "./mlmc";
import { fmtCount, runStatus } from "./montecarlo";
import { isPlate, structureText } from "./structure";
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
  const [solver, variance, mean, samples] = fig.panels(4, { cols: 3, wideFirst: true });
  const lv = Math.max(1, Math.min(st.liveLevel, c.levels - 1));
  const acc = c.run.streams[lv].acc;
  const live = liveLevels(c), tol = tolAt(c);
  const scale = scaleOf(c, live), rows = live.map(statsOf);
  drawVariance(variance, c.levels, rows, scale, lv);
  drawMean(mean, c.levels, rows, scale, lv);
  drawLiveSamples(samples, c, live, lv, tol);
  if (acc.n < 1) {
    solver.draw({ xlabel: "", ylabel: "", series: [], title: `waiting for the first samples of level ${lv}…` });
    return { readout: summary(c, lv, null, tol, live), animate: true };
  }

  const i = Math.floor(t / DWELL_MS) % Math.min(acc.n, CYCLE);
  const sampler = samplers(JSON.stringify([c.spec, c.M]), () => new Sampler(c.spec, c.kl, c.klY));
  const pair = pairs(JSON.stringify([c.spec, c.M, c.neOf(0), lv, i]), () => solvePair(c, sampler, lv, i));
  const stored = { Q: acc.values[i], Qc: acc.coarseValues[i] };
  const th = (2 * Math.PI * t) / PERIOD_MS;

  if (isPlate(st)) drawPlates(solver, c, lv, pair, th);
  else drawBeams(solver, c, lv, pair, th);
  return { readout: summary(c, lv, { pair, stored }, tol, live), animate: true };
}

/** The tolerance the pane picks, as it stands: its N_ℓ (settled, or asked for so far) and standard MC's count and cost beside them. */
interface Tol {
  readonly k: number;
  readonly eps: number;
  readonly N: readonly number[];
  readonly L: number;
  readonly mlCost: number;
  readonly nMC: number;
  readonly mcCost: number;
  readonly settled: boolean;
}

function tolAt(c: MlmcSession): Tol | null {
  const k = tolIndex(c.st), alg = c.sweep?.runs[k];
  if (!alg) return null;
  const settled = alg.status === "converged";
  const N = (settled ? alg.N : alg.wants()).slice(), L = N.length - 1;
  if (L >= c.S.length) return null;
  const mc = mcCost(c, alg);
  return {
    k, eps: SWEEP[k] * c.st.mlEps, N, L, settled,
    mlCost: N.reduce((s, n, l) => s + n * c.cost(l), 0), nMC: Math.ceil(mc / c.dofs(L)), mcCost: mc,
  };
}

/** What a level's stream says now, from every sample it has — the live counterpart of the survey's fixed prefix. */
interface LiveLevel extends LevelStats {
  readonly n: number;
}

function liveLevels(c: MlmcSession): LiveLevel[] {
  return c.run.streams.slice(0, c.levels).map(({ acc }) =>
    ({ n: acc.n, meanQ: acc.q.mean, varQ: acc.q.variance, meanY: acc.y.mean, varY: acc.y.variance }));
}

/** A level's moments for the plots: none until it has two samples. */
const statsOf = (v: LiveLevel): LevelStats => (v.n >= 2 ? v : { meanQ: NaN, varQ: NaN, meanY: NaN, varY: NaN });

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
function liveEstimate(live: readonly LiveLevel[], tol: Tol | null): { L: number; est: number; se: number } | null {
  let L = tol ? Math.min(tol.L, live.length - 1) : live.length - 1;
  while (L >= 0 && !(live[L].n >= 2)) L--;
  if (L < 0) return null;
  let est = 0, v = 0;
  for (let l = 0; l <= L; l++) { est += live[l].meanY; v += live[l].varY / live[l].n; }
  return { L, est, se: Math.sqrt(v) };
}

/**
 * Samples per level, as in the MLMC view's panel: what the chosen tolerance
 * asks for as blue bars (hollow while settling), standard MC's computed count
 * on level L as a hollow orange bar, and the samples in hand as a line over them.
 */
function drawLiveSamples(plot: Plot, c: MlmcSession, live: readonly LiveLevel[], lv: number, tol: Tol | null): void {
  const l = live.map((_, k) => k), n = live.map((v) => (v.n > 0 ? v.n : NaN));
  const series: Series[] = [];
  if (tol) {
    series.push({
      label: tol.settled ? "MLMC asks" : "MLMC asks (settling)", ...multi, bars: true, hollow: !tol.settled,
      x: tol.N.map((_, k) => k), y: tol.N.slice(),
    });
    if (Number.isFinite(tol.nMC)) series.push({ label: "standard MC (computed)", ...single, bars: true, hollow: true, x: [tol.L], y: [tol.nMC] });
  }
  series.push({ label: "in hand", x: l, y: n, color: INK, width: 2, markers: true });
  series.push({ label: "survey", x: [-0.3, c.levels - 0.7], y: [c.st.mlSurvey, c.st.mlSurvey], color: MUTED, width: 1, dash: [3, 4], inert: true });
  const total = live.reduce((s, v) => s + v.n, 0);
  plot.draw({
    ...levelAxis(c.levels), ylog: true, series,
    heading: tol ? `samples per level, ε = ${tol.eps.toExponential(1)}` : "samples per level",
    note: `${fmtCount(total)} in hand`,
    over: (ctx, a) => {
      ctx.strokeStyle = MUTED;
      ctx.lineWidth = 1;
      ctx.setLineDash([3, 4]);
      ctx.beginPath(); ctx.moveTo(a.sx(lv), a.box.t); ctx.lineTo(a.sx(lv), a.box.b); ctx.stroke();
      ctx.setLineDash([]);
    },
    hover: (s, i) => `${s.label}\nlevel ${s.x[i]}: N = ${s.y[i].toLocaleString()}` + (s.x[i] === lv ? "\nthe level shown above" : ""),
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

function summary(
  c: MlmcSession, lv: number, shown: { pair: Pair; stored: { Q: number; Qc: number } } | null, tol: Tol | null,
  live: readonly LiveLevel[],
): string {
  const st = c.st, d = displayOf(st), v = (x: number) => valueOf(d, st.qoi, x, st.load);
  const elems = (n: number) => (isPlate(st) ? `${n} × ${n}` : String(n));
  const lines = [
    `${structureText(st)}, p = ${st.p} ${continuityName(continuityK(st.continuity, st.p))}; Q = ${QOI_NAME[st.qoi]}; ` +
      `${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L, σ = ${plain(st.sigma, 3)}, M = ${c.M}; seed ${st.seed}`,
  ];
  const n = c.run.streams.map((s) => s.acc.n);
  lines.push(`${runStatus(c.run)}: ${n.reduce((x, y) => x + y, 0).toLocaleString()} samples (${n.map(fmtCount).join(" · ")} by level) on ${workers.size} workers`);
  if (shown) {
    const { pair: p, stored } = shown, Y = p.fine.Q - p.coarse.Q;
    const same = p.fine.Q === stored.Q && p.coarse.Q === stored.Qc;
    lines.push(
      "",
      `shown: sample ${p.i} of level ${lv}, ${elems(c.neOf(lv))} and ${elems(c.neOf(lv - 1))} elements, the same ω` +
        (same ? " (re-solved here: matches the workers bit for bit)" : ` (re-solved here: differs from the workers by ${(p.fine.Q - stored.Q).toExponential(2)}, ${(p.coarse.Q - stored.Qc).toExponential(2)})`),
      table(["", "Q", "|Y| / |Q|"], [
        [`Q_${lv}`, v(p.fine.Q), ""],
        [`Q_${lv - 1}`, v(p.coarse.Q), ""],
        [`Y_${lv} = Q_${lv} − Q_${lv - 1}`, v(Y), (Math.abs(Y) / Math.abs(p.fine.Q)).toExponential(2)],
      ]),
    );
  }
  if (tol) {
    const ratio = `${plain(tol.mcCost / tol.mlCost, 3)}×`;
    lines.push("", `ε = ${tol.eps.toExponential(1)}${tol.settled ? "" : " (settling)"}: MLMC N_ℓ = ${tol.N.map(fmtCount).join(" · ")}; ` +
      `Std MC would need ${fmtCount(tol.nMC)} on level ${tol.L}, ${ratio} the cost`);
  }
  if (live.filter((v) => v.n >= 2).length > 1) {
    const alpha = -slope(live.map((v) => (v.n >= 2 ? Math.abs(v.meanY) : NaN)));
    const beta = -slope(live.map((v) => (v.n >= 2 ? v.varY : NaN)));
    const gamma = slope(live.map((_, l) => c.cost(l)));
    const now = liveEstimate(live, tol), scale = scaleOf(c, live);
    lines.push(`from every sample so far: α = ${fmt2(alpha)}  β = ${fmt2(beta)}  γ = ${fmt2(gamma)}` +
      (now ? `; estimate over levels 0 … ${now.L}: ${v(now.est)} ± ${(Z95 * now.se / scale).toExponential(1)} (95%, relative)` : ""));
  }
  const runs = c.sweep?.runs ?? [];
  if (runs.length) {
    lines.push(`tolerances (θ = ${THETA}): ` + runs.map((r, k) => `ε = ${(SWEEP[k] * st.mlEps).toExponential(1)} ${r.status}` +
      (r.status === "converged" ? ` (L = ${r.L})` : "")).join("; "));
  }
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}
