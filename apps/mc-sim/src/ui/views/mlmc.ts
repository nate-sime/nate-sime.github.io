// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 6 on screen: multilevel Monte Carlo, as Giles' `mlmc_test` reports it.
 *
 * One experiment runs in the worker pool: a survey of every level of the
 * hierarchy at a fixed number of samples, and the adaptive algorithm at five
 * tolerances, ε_min·{16, 8, 4, 2, 1}. All of it reads one store of samples per
 * level — level ℓ's samples come from random stream ℓ, solved on ℓ and on its
 * parent with the same ω — so the five runs cost hardly more than the finest.
 * Tolerances are relative to |E[Q]| as the survey's level 0 estimates it: a
 * number fixed by the seed before any tolerance starts, so the runs stay
 * reproducible. (Q of the mean beam would not do: for the forced response it is
 * half of E[Q].)
 *
 * The figure is `mlmc_plot`'s six panels, drawn as instrument panels: a
 * heading, the rate it measures at the top right, and a strip of keys. Colour
 * and marker are the estimator, in every panel: orange squares are standard
 * Monte Carlo and everything it sees (Q_ℓ), blue circles are MLMC and
 * everything it sees (Y_ℓ = Q_ℓ − Q_ℓ₋₁). Dashed lines are predictions, and
 * hollow marks are computed or still settling rather than measured.
 *
 * - (a) variance and (b) |mean| against level, of Q_ℓ and of Y_ℓ, from the
 *   survey, on a log₂ axis: the blue lines fall with slopes −β and −α, the
 *   orange ones stay put;
 * - (c) the consistency check and (d) the kurtosis of Y_ℓ, Giles' two
 *   diagnostics of the coupling and of how far V_ℓ can be trusted;
 * - (e) samples per level for the tolerance the pane shows: MLMC's N_ℓ as
 *   bars, beside the one bar standard MC would need on level L, with the other
 *   tolerances' N_ℓ faint behind;
 * - (f) cost against ε: both estimators' costs predicted from the survey
 *   across the whole range, dashed, with every tolerance's run on them. With
 *   β > γ the MLMC line falls as ε⁻², as if the problem had no mesh at all,
 *   while standard MC's falls faster.
 *
 * The readout is `mlmc_test`'s printout: the survey table, α, β, γ, and one
 * row per tolerance.
 */

import { MlmcSweep, THETA, slope, survey, type Mlmc, type SurveyLevel } from "../../mc/mlmc";
import type { McRun, StreamSpec } from "../../mc/pool";
import type { McSpec } from "../../mc/sampler";
import type { KLData } from "../../random/kl";
import { Z95 } from "../../mc/stats";
import type { Figure } from "../figure";
import { SLOT, decade, type Plot, type PlotSpec, type Series } from "../plot";
import { continuityK, continuityName, displayOf, type State } from "../state";
import { plain, referenceNote } from "../units";
import { QOI_NAME, valueOf } from "./convergence";
import { KERNEL_NAME } from "./field";
import {
  PLATE_MAX_NE, admissibleAt, isPlate, klsOf, levelsWithin, mcSpecOf, structureText, workAt, workUnit,
} from "./structure";
import { fmtCount, runStatus } from "./montecarlo";
import { memo, table, type ViewResult } from "./view";
import { workers } from "./workers";

/** The tolerances of a sweep, as multiples of the finest. */
export const SWEEP = [16, 8, 4, 2, 1];
/** Warm-up samples on levels 0, 1, 2 before the first allocation. */
const N0 = 200;
/** Most samples a tolerance may ask of one level: 2·10⁶ is 32 MB of Q and Q_c, and minutes of work. */
const MAX_N = 2e6;

const sweeps = memo<MlmcSweep>();

/**
 * The multilevel experiment for the pane's settings: its run in the pool
 * (opened, its demand set), its sweep, and the survey so far — or why there
 * can be none. The MLMC view and the live view (stage 8) both call it, so they
 * drive the same run and the same store of samples.
 */
export function mlmcSession(st: State): MlmcSession | string {
  const p = st.p, k = continuityK(st.continuity, p), ne0 = st.ne0;
  const levels = isPlate(st) ? levelsWithin(ne0, st.levels, PLATE_MAX_NE) : st.levels;
  const why = admissibleAt(st, ne0);
  if (why) return `${continuityName(Math.max(k, 0))}, degree ${p}, ${ne0} elements: ${why}`;
  if (levels < 3)
    return "Multilevel Monte Carlo fits its rates on levels 1 … L and starts from levels 0, 1, 2: raise the number of levels to at least 3.";

  const { kl, klY, M } = klsOf(st);
  const spec = mcSpecOf(st);
  const neOf = (l: number) => ne0 * 2 ** l;
  /** One solve on level ℓ, in the cost model's unit (`workUnit`). */
  const dofs = (l: number) => (l < 0 ? 0 : workAt(st, neOf(l)));
  /** Cost of one sample of level ℓ: its fine and its coarse solve. */
  const cost = (l: number) => dofs(l) + dofs(l - 1);
  const key = JSON.stringify(["mlmc", spec, M, ne0, levels, st.mlSurvey, st.mlEps]);
  const run = workers.pool.open(key, () => ({ spec, kl, klY }));
  const streamOf = (l: number): StreamSpec => ({ ne: neOf(l), coarse: l > 0, wantW: false, stream: l });
  for (let l = 0; l < levels; l++) run.stream(streamOf(l));
  // The scale ε is relative to: |E[Q₀]| over the survey's samples of level 0, once they are in.
  const level0 = run.streams[0].acc;
  const scale = level0.n >= st.mlSurvey ? Math.abs(level0.prefix(st.mlSurvey).q.mean) || 1 : NaN;
  const sweep = Number.isFinite(scale)
    ? sweeps(`${run.id}`, () => new MlmcSweep(levels, st.mlSurvey, SWEEP.map((f) => f * st.mlEps * scale),
      { N0, Lmin: 2, cost }, MAX_N))
    : null;
  const demand = sweep ? sweep.step((l) => run.streams[l].acc) : run.streams.map(() => st.mlSurvey);
  demand.forEach((n, l) => workers.pool.demand(run, run.streams[l], n));

  // The survey, over the levels that have their samples (a prefix: level ℓ's check reads ℓ − 1).
  let complete = 0;
  while (complete < levels && run.streams[complete].acc.n >= st.mlSurvey) complete++;
  const S = survey(run.streams.slice(0, complete).map((s) => ({ fine: s.acc.values, coarse: s.acc.coarseValues })), st.mlSurvey);
  return { st, levels, scale, S, sweep, cost, dofs, neOf, run, M, spec, kl, klY };
}

export function renderMlmc(fig: Figure, st: State): ViewResult {
  const ctx = mlmcSession(st);
  const empty = () => fig.panels(1)[0].draw({ xlabel: "level ℓ", ylabel: "", series: [] });
  if (typeof ctx === "string") {
    empty();
    return { readout: ctx, animate: false };
  }
  const { run, S, M } = ctx;
  if (!S.length) {
    empty();
    return { readout: readout(ctx, M), animate: run.running };
  }
  const [variance, mean, consistency, kurtosis, samples, costs] = fig.panels(6, { cols: 2 });
  const scale = ctx.scale, levels = ctx.levels;
  drawVariance(variance, levels, S, scale);
  drawMean(mean, levels, S, scale);
  drawConsistency(consistency, ctx);
  drawKurtosis(kurtosis, ctx);
  drawSamples(samples, ctx);
  drawCost(costs, ctx);
  return { readout: readout(ctx, M), animate: run.running };
}

export interface MlmcSession {
  readonly st: State;
  /** Levels in use: the pane's, or fewer for a plate, whose finest stops at PLATE_MAX_NE. */
  readonly levels: number;
  readonly scale: number;
  readonly S: SurveyLevel[];
  readonly sweep: MlmcSweep | null;
  readonly cost: (l: number) => number;
  readonly dofs: (l: number) => number;
  readonly neOf: (l: number) => number;
  readonly run: McRun;
  /** KL terms, the spec every sample is solved from, and the expansions it reads. */
  readonly M: number;
  readonly spec: McSpec;
  readonly kl: KLData;
  readonly klY?: KLData;
}

type Ctx = MlmcSession;

/**
 * A level's moments as the level plots read them: the survey's fixed prefix
 * (MLMC view), or every sample so far (live view). NaN where a level has none.
 */
export interface LevelStats {
  readonly meanQ: number;
  readonly varQ: number;
  readonly meanY: number;
  readonly varY: number;
}

/**
 * Standard MC and what it sees (Q_ℓ) orange with square markers; MLMC and what
 * it sees (Y_ℓ) blue with round ones — the same in every panel, so "the orange
 * line is flat" and "the orange bar is tall" are statements about one estimator.
 */
export const SINGLE = SLOT[1];
export const MULTI = SLOT[0];
const DASH = [6, 4];
const MUTED = "rgba(207, 238, 255, 0.35)";
export const single = { color: SINGLE, marker: "square" } as const;
export const multi = { color: MULTI, marker: "circle" } as const;

/** The level axis, and the instrument-panel style every MLMC panel shares. */
export const levelAxis = (levels: number) => ({
  xlabel: "level ℓ", ylabel: "", legend: "top", xlim: [-0.3, levels - 0.7] as const,
  xticks: Array.from({ length: levels }, (_, l) => l), xfmt: String,
}) satisfies Partial<PlotSpec>;

/** The tolerance the pane shows, as an index into SWEEP. */
export const tolIndex = (st: State) => Math.max(0, Math.min(SWEEP.length - 1, st.cmpTol));

const levelsOf = (rows: readonly unknown[]) => rows.map((_, l) => l);
export const fmt2 = (x: number) => (Number.isFinite(x) ? x.toFixed(2) : "—");

/** A dotted vertical line at level `lv`: the live view's shown level. */
function markLevel(lv: number | undefined) {
  return (ctx: CanvasRenderingContext2D, a: { sx: (x: number) => number; box: { t: number; b: number } }) => {
    if (lv === undefined) return;
    ctx.strokeStyle = MUTED;
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 4]);
    ctx.beginPath(); ctx.moveTo(a.sx(lv), a.box.t); ctx.lineTo(a.sx(lv), a.box.b); ctx.stroke();
    ctx.setLineDash([]);
  };
}

/** (a) V[Q_ℓ] orange, V[Q_ℓ − Q_ℓ₋₁] blue, on a log₂ axis; the blue line's slope is −β. Relative to E[Q]². */
export function drawVariance(plot: Plot, levels: number, rows: readonly LevelStats[], scale: number, mark?: number): void {
  const l = levelsOf(rows), s2 = scale * scale;
  const beta = -slope(rows.map((r) => r.varY));
  plot.draw({
    ...levelAxis(levels), heading: "variance per level, ÷ E[Q]²", note: `β ≈ ${fmt2(beta)}`, ylog: true, ybase: 2,
    series: [
      { label: "V[Q_ℓ] — standard MC", ...single, x: l, y: rows.map((r) => r.varQ / s2), width: 2, markers: true },
      { label: "V[Q_ℓ − Q_ℓ₋₁] — MLMC", ...multi, x: l, y: rows.map((r, i) => (i === 0 ? NaN : r.varY / s2)), width: 2, markers: true },
    ],
    over: markLevel(mark),
    hover: (s, i) => `${s.label.split(" — ")[0]}, level ${i}\n${s.y[i].toExponential(2)} × E[Q]²  (2^${Math.log2(s.y[i]).toFixed(1)})`,
  });
}

/** (b) |E[Q_ℓ]| orange, |E[Q_ℓ − Q_ℓ₋₁]| blue, on a log₂ axis; the blue line's slope is −α. Relative to |E[Q]|. */
export function drawMean(plot: Plot, levels: number, rows: readonly LevelStats[], scale: number, mark?: number): void {
  const l = levelsOf(rows);
  const alpha = -slope(rows.map((r) => Math.abs(r.meanY)));
  plot.draw({
    ...levelAxis(levels), heading: "|mean| per level, ÷ |E[Q]|", note: `α ≈ ${fmt2(alpha)}`, ylog: true, ybase: 2,
    series: [
      { label: "|E[Q_ℓ]| — standard MC", ...single, x: l, y: rows.map((r) => Math.abs(r.meanQ) / scale), width: 2, markers: true },
      { label: "|E[Q_ℓ − Q_ℓ₋₁]| — MLMC", ...multi, x: l, y: rows.map((r, i) => (i === 0 ? NaN : Math.abs(r.meanY) / scale)), width: 2, markers: true },
    ],
    over: markLevel(mark),
    hover: (s, i) => `${s.label.split(" — ")[0]}, level ${i}\n${s.y[i].toExponential(2)} × |E[Q]|  (2^${Math.log2(s.y[i]).toFixed(1)})`,
  });
}

/** (c) |E[Y_ℓ] + E[Q_ℓ₋₁] − E[Q_ℓ]| over three standard errors: below 1 unless the coupling is broken. */
function drawConsistency(plot: Plot, { levels, S }: Ctx): void {
  const pts = S.filter((s) => s.level > 0);
  const top = Math.max(1.5, ...pts.map((s) => s.consistency).filter(Number.isFinite)) * 1.15;
  const worst = Math.max(...pts.map((s) => s.consistency).filter(Number.isFinite));
  plot.draw({
    ...levelAxis(levels), heading: "consistency check", ylim: [0, top],
    note: Number.isFinite(worst) ? `worst ${worst.toFixed(2)} · below 1 if coupled` : "",
    series: [
      { label: "|E[Y_ℓ] + E[Q_ℓ₋₁] − E[Q_ℓ]| / 3σ", ...multi, x: pts.map((s) => s.level), y: pts.map((s) => s.consistency), width: 2, markers: true },
      { label: "1", x: [-0.3, levels - 0.7], y: [1, 1], color: MUTED, width: 1, dash: [3, 4], inert: true, unlisted: true },
    ],
    hover: (_, i) => `level ${pts[i].level}\ncheck = ${pts[i].consistency.toFixed(2)}`,
  });
}

/** (d) Kurtosis of Y_ℓ: large means V_ℓ rests on a few rare samples. */
function drawKurtosis(plot: Plot, { levels, S }: Ctx): void {
  const pts = S.filter((s) => s.level > 0);
  const top = Math.max(3, ...pts.map((s) => s.kurtosis).filter(Number.isFinite)) * 1.15;
  const most = Math.max(...pts.map((s) => s.kurtosis).filter(Number.isFinite));
  plot.draw({
    ...levelAxis(levels), heading: "kurtosis of Y_ℓ", ylim: [0, top],
    note: Number.isFinite(most) ? `largest ${most.toFixed(1)} · 3 if Gaussian` : "",
    series: [{ label: "Q_ℓ − Q_ℓ₋₁", ...multi, x: pts.map((s) => s.level), y: pts.map((s) => s.kurtosis), width: 2, markers: true, unlisted: true }],
    hover: (_, i) => `level ${pts[i].level}\nkurtosis = ${pts[i].kurtosis.toFixed(1)}`,
  });
}

/** The N_ℓ a tolerance settled on, or asks for so far. */
const nOf = (a: Mlmc) => (a.status === "converged" || a.status === "failed" ? a.N : a.wants());

/** γ: the slope of log₂ C_ℓ, the cost of a sample, over the surveyed levels. */
const gammaOf = (c: Ctx) => slope(c.S.map((_, l) => c.cost(l)));

/**
 * (e) Samples per level for the tolerance the pane shows: MLMC's N_ℓ as blue
 * bars (hollow while settling), beside the orange bar of the samples standard
 * MC would need on level L (hollow: computed, never run). The other
 * tolerances' N_ℓ are faint lines behind.
 */
function drawSamples(plot: Plot, c: Ctx): void {
  const { st, levels, sweep } = c, k = tolIndex(st), runs = sweep?.runs ?? [];
  const eps = (i: number) => (SWEEP[i] * st.mlEps).toExponential(1);
  const series: Series[] = [];
  /** Each faint line's tolerance, for its hover: they share one legend key. */
  const other = new Map<Series, number>();
  runs.forEach((a, i) => {
    if (i === k) return;
    const N = nOf(a);
    const s: Series = {
      label: "other ε", x: N.map((_, l) => l), y: N.slice(), ...multi, width: 1, alpha: 0.3,
      unlisted: other.size > 0, hollow: !a.done,
    };
    other.set(s, i);
    series.push(s);
  });
  const shown = runs[k];
  if (shown) {
    const N = nOf(shown), settled = shown.status === "converged";
    series.push({
      label: settled ? "MLMC" : `MLMC (${shown.done ? shown.status : "settling"})`, ...multi, bars: true, hollow: !settled,
      x: N.map((_, l) => l), y: N.slice(),
    });
    const L = N.length - 1, mc = L < c.S.length ? mcCost(c, shown) / c.dofs(L) : NaN;
    if (Number.isFinite(mc)) series.push({ label: "standard MC (computed)", ...single, bars: true, hollow: true, x: [L], y: [Math.ceil(mc)] });
  }
  const gamma = gammaOf(c);
  plot.draw({
    ...levelAxis(levels), heading: shown ? `samples per level, ε = ${eps(k)}` : "samples per level",
    note: Number.isFinite(gamma) ? `γ ≈ ${fmt2(gamma)}` : "", ylog: true, series,
    title: shown ? undefined : "waiting for the survey…",
    hover: (s, i) => {
      const j = other.get(s);
      if (j !== undefined) return `MLMC, ε = ${eps(j)}\nlevel ${i}: ${fmtCount(s.y[i])} samples`;
      if (s.color === SINGLE) return `standard MC, ε = ${eps(k)}\nlevel ${s.x[i]}: ${fmtCount(s.y[i])} samples, computed from V[Q_L]`;
      return `MLMC, ε = ${eps(k)}\nlevel ${i}: ${fmtCount(s.y[i])} samples · cost ${fmtCount(c.cost(i))} each`;
    },
  });
}

/** Plain Monte Carlo to the same ε on level L: V[Q_L] / ((1 − θ)ε²) samples of one fine solve each. */
export function mcCost(c: Ctx, a: Mlmc): number {
  const { S } = c, L = Math.min(a.L, S.length - 1);
  return L < 0 ? NaN : (S[L].varQ * c.dofs(a.L)) / ((1 - a.theta) * a.opts.eps ** 2);
}

/**
 * Both estimators' costs at an absolute tolerance `eps`, predicted from the
 * survey as Giles' algorithm would settle them: the first level L ≥ 2 whose
 * bias estimate is within √θ ε, the optimal N_ℓ unrounded, and standard MC's
 * V[Q_L]/((1 − θ)ε²) solves on L. Null where even the finest level is too
 * coarse.
 */
export function predictCost(c: Ctx, eps: number): { L: number; mlmc: number; mc: number } | null {
  const { S } = c;
  const m = S.map((s) => Math.abs(s.meanY)), V = S.map((s) => s.varY);
  const alpha = Math.max(0.5, -slope(m)), budget = (1 - THETA) * eps * eps;
  for (let L = 2; L < S.length; L++) {
    const bias = Math.max(...[0, 1, 2].filter((i) => L - i >= 1).map((i) => m[L - i] * 2 ** (-alpha * i))) / (2 ** alpha - 1);
    if (!(bias <= Math.sqrt(THETA) * eps)) continue;
    let sum = 0;
    for (let l = 0; l <= L; l++) sum += Math.sqrt(V[l] * c.cost(l));
    return { L, mlmc: (sum * sum) / budget, mc: (S[L].varQ * c.dofs(L)) / budget };
  }
  return null;
}

/** Least-squares slope of log y against log x: a cost's power of ε. */
function logSlope(x: readonly number[], y: readonly number[]): number {
  if (x.length < 2) return NaN;
  const X = x.map(Math.log), Y = y.map(Math.log);
  const mx = X.reduce((s, v) => s + v, 0) / X.length, my = Y.reduce((s, v) => s + v, 0) / Y.length;
  let sxy = 0, sxx = 0;
  X.forEach((v, i) => { sxy += (v - mx) * (Y[i] - my); sxx += (v - mx) ** 2; });
  return sxy / sxx;
}

/**
 * (f) Cost against ε. Dashed: both estimators' costs predicted from the survey,
 * across the sweep and a factor 2 either side. Markers: each tolerance's run on
 * them — MLMC's Σ N_ℓ C_ℓ (hollow while it settles), and what standard MC
 * would cost for it (hollow: computed, not run). The note gives the slopes the
 * converged markers make, against ε⁻² for MLMC when β > γ.
 */
function drawCost(plot: Plot, c: Ctx): void {
  const { st, sweep, scale } = c;
  const rel = SWEEP.map((f) => f * st.mlEps);
  const lo = Math.min(...rel) / 2, hi = Math.max(...rel) * 2, n = Math.ceil(8 * Math.log10(hi / lo));
  const grid = Array.from({ length: n + 1 }, (_, k) => lo * (hi / lo) ** (k / n));
  const pr = c.S.length >= 3 ? grid.map((e) => predictCost(c, e * scale)) : [];
  const runs = (sweep?.runs ?? []).map((a, i) => ({ a, e: rel[i] })).filter(({ a }) => a.status === "converged" || !a.done);
  const done = runs.filter(({ a }) => a.status === "converged");
  const wanted = (a: Mlmc) => a.wants().reduce((s, n, l) => s + n * c.cost(l), 0);

  const mcLine: Series = { label: "standard MC", ...single, x: grid, y: pr.map((p) => p?.mc ?? NaN), width: 2, dash: DASH };
  const mlLine: Series = { label: "MLMC", ...multi, x: grid, y: pr.map((p) => p?.mlmc ?? NaN), width: 2, dash: DASH };
  const mcRuns: Series = {
    label: "standard MC", ...single, width: 0, markers: true, hollow: true, unlisted: true,
    x: done.map(({ e }) => e), y: done.map(({ a }) => mcCost(c, a)),
  };
  const mlRuns: Series = {
    label: "MLMC", ...multi, width: 0, markers: true, unlisted: true, hollow: runs.map(({ a }) => !a.done),
    x: runs.map(({ e }) => e), y: runs.map(({ a }) => (a.done ? a.cost : wanted(a))),
  };
  const slopes = [["MLMC", logSlope(done.map(({ e }) => e), done.map(({ a }) => a.cost))], ["MC", logSlope(mcRuns.x as number[], mcRuns.y as number[])]] as const;
  const fitted = slopes.filter(([, s]) => Number.isFinite(s)).map(([m, s]) => `${m} ${s.toFixed(1).replace("-", "−")}`);
  plot.draw({
    heading: `cost to reach ε  [${workUnit(st)}]`, legend: "top",
    note: fitted.length ? `measured slope ${fitted.join(" · ")}` : "dashed predicted · marks run",
    title: pr.length ? undefined : "waiting for the survey…",
    xlabel: "tolerance ε (relative)", ylabel: "", xlog: true, ylog: true, xlim: [lo, hi], xfmt: decade,
    series: [mcLine, mlLine, mcRuns, mlRuns],
    hover: (s, i) => {
      const eps = s.x[i].toExponential(1), cost = fmtCount(s.y[i]);
      if (s === mcLine) return `standard MC, predicted\nε = ${eps}: cost ${cost} · ${fmtCount(Math.ceil(s.y[i] / c.dofs(pr[i]!.L)))} samples on level ${pr[i]!.L}`;
      if (s === mlLine) return `MLMC, predicted\nε = ${eps}: cost ${cost} · levels 0–${pr[i]!.L}`;
      if (s === mcRuns) return `standard MC, computed\nε = ${eps}: cost ${cost} · level ${done[i].a.L}`;
      const { a } = runs[i];
      return a.done ? `MLMC, run\nε = ${eps}: cost ${cost} · levels 0–${a.L}` : `MLMC, settling\nε = ${eps}: asking ${cost} so far`;
    },
  });
}

/** Giles' `mlmc_test` printout: the survey table, the rates, and one row per tolerance. */
function readout(c: Ctx, M: number): string {
  const { st, scale, S, sweep, run } = c, d = displayOf(st);
  const v = (x: number) => valueOf(d, st.qoi, x, st.load);
  const levels = c.levels;
  const lines = [
    `${structureText(st)}, p = ${st.p} ${continuityName(continuityK(st.continuity, st.p))}; Q = ${QOI_NAME[st.qoi]}; ` +
      `levels 0 … ${levels - 1}, ne = ${st.ne0}·2^ℓ${isPlate(st) ? " per side" : ""}; ` +
      `${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L, σ = ${plain(st.sigma, 3)}, M = ${M}` +
      `${st.massFollows ? ", mass follows depth" : ""}${st.loadSigma > 0 ? `, load σ_q = ${plain(st.loadSigma, 3)}` : ""}; seed ${st.seed}`,
  ];
  if (levels < st.levels) lines.push(`(${st.levels - levels} of the levels asked for left out: a plate stops at ${PLATE_MAX_NE} × ${PLATE_MAX_NE} elements)`);
  const n = run.streams.map((s) => s.acc.n), total = n.reduce((a, b) => a + b, 0), wall = run.wallMs / 1000;
  lines.push(`${runStatus(run)}: ${total.toLocaleString()} samples (${n.map(fmtCount).join(" · ")} by level) on ${workers.size} workers` +
    (wall > 0 ? `, ${plain(wall, 3)} s` : ""));
  if (!S.length) return [...lines, "", "the survey is waiting for its first level…"].join("\n");
  if (!sweep) return lines.join("\n");

  const rel = (x: number) => x.toExponential(2);
  const rows = S.map((s, l) => [
    String(l), rel(s.meanY / scale), v(s.meanQ), rel(s.varY / scale ** 2), rel(s.varQ / scale ** 2),
    l === 0 ? "—" : Number.isFinite(s.kurtosis) ? s.kurtosis.toFixed(1) : "—",
    Number.isFinite(s.consistency) ? s.consistency.toFixed(2) : "—",
    fmtCount(c.cost(l)),
  ]);
  const alpha = -slope(S.map((s) => Math.abs(s.meanY))), beta = -slope(S.map((s) => s.varY));
  const gamma = slope(S.map((_, l) => c.cost(l)));
  lines.push(
    "",
    `survey, ${st.mlSurvey.toLocaleString()} samples per level; Y_ℓ = Q_ℓ − Q_ℓ₋₁ (Y₀ = Q₀); means relative to |E[Q₀]| = ${v(scale)}, variances to its square`,
    table(["ℓ", "E[Y_ℓ]", "E[Q_ℓ]", "V[Y_ℓ]", "V[Q_ℓ]", "kurtosis", "check", "cost C_ℓ"], rows),
    `α = ${fmt2(alpha)}  β = ${fmt2(beta)}  γ = ${fmt2(gamma)}  (cost in ${workUnit(st)})`,
  );

  const trows = sweep.runs.map((a, i) => {
    const e = SWEEP[i] * st.mlEps, ok = a.status === "converged";
    const mc = ok ? mcCost(c, a) : NaN;
    return [
      e.toExponential(1),
      ok ? `${v(a.estimate)} ± ${(Z95 * Math.sqrt(a.sampleVariance) / scale).toExponential(1)}` : a.status === "sampling" ? "running" : a.status,
      ok ? fmtCount(a.cost) : "—",
      Number.isFinite(mc) ? fmtCount(mc) : "—",
      Number.isFinite(mc) ? `${(mc / a.cost).toFixed(1)}×` : "—",
      a.N.length ? a.N.map(fmtCount).join(" ") : "—",
    ];
  });
  lines.push(
    "",
    `MLMC, θ = ${THETA}: bias² ≤ θε², variance ≤ (1 − θ)ε²; value ± 95% of the sampling error`,
    table(["ε", "value", "MLMC cost", "Std MC cost", "saving", "N_ℓ"], trows),
  );
  if (sweep.runs.some((a) => a.status === "failed"))
    lines.push(`failed: the bias test asked for a level finer than ℓ = ${levels - 1}; raise the number of levels or the tolerance`);
  const over = sweep.runs.filter((a) => a.status === "over budget");
  if (over.length) {
    const a = over[over.length - 1];
    lines.push(`over budget: ε = ${(SWEEP[sweep.runs.indexOf(a)] * st.mlEps).toExponential(1)} would need ${fmtCount(Math.max(...a.wants()))} samples on one level (limit ${fmtCount(MAX_N)})`);
  }
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}

