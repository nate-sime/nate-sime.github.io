// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 6 on screen: multilevel Monte Carlo, as Giles' `mlmc_test` lays it out.
 *
 * One experiment runs in the worker pool: a survey of every level of the
 * hierarchy at a fixed number of samples, and the adaptive algorithm at five
 * tolerances, ε_min·{16, 8, 4, 2, 1}. All of it reads one store of samples per
 * level — level ℓ's samples come from random stream ℓ, solved on ℓ and on its
 * parent with the same ω — so the five runs cost hardly more than the finest.
 * Tolerances are relative to |E[Q]| as the survey's level 0 estimates it: a
 * number fixed by the seed before any tolerance starts, so the runs stay
 * reproducible. (Q of the mean beam would not do: for the forced response it is
 * half of E[Q].) Four pictures, side by side:
 *
 * - mean vs level: |E[Q_ℓ]| and |E[Y_ℓ]|, Y_ℓ = Q_ℓ − Q_{ℓ−1}, with their 95%
 *   intervals: the corrections fall as 2^{−αℓ} — the bias the hierarchy leaves.
 * - variance vs level: V[Q_ℓ] stays put, V[Y_ℓ] falls as 2^{−βℓ}. The gap
 *   between them is what MLMC saves: a correction is far cheaper to estimate
 *   than Q itself.
 * - samples per level: the N_ℓ each tolerance settled on — most samples where
 *   they are cheapest, a handful on the finest levels, and a level more each
 *   time ε falls far enough that the bias asks for it.
 * - cost vs ε: ε² × cost for MLMC and for plain Monte Carlo on the finest level
 *   the same tolerance needed. With β > γ the MLMC line is flat — O(ε⁻²), as if
 *   the problem had no mesh at all — while plain Monte Carlo's climbs.
 */

import { SUPPORTS, admissible, freeDofs } from "../../beam/beam";
import { MlmcSweep, THETA, slope, survey, type Mlmc, type SurveyLevel } from "../../mc/mlmc";
import type { McRun, StreamSpec } from "../../mc/pool";
import type { McSpec } from "../../mc/sampler";
import { Z95 } from "../../mc/stats";
import { termsOf } from "../../random/kl";
import type { Figure } from "../figure";
import { SLOT, decade, type Plot, type Series } from "../plot";
import { beamCaseOf, continuityK, continuityName, displayOf, fieldSpecOf, type State } from "../state";
import { plain, referenceNote } from "../units";
import { QOI_NAME, valueOf } from "./convergence";
import { KERNEL_NAME, klOf } from "./field";
import { fmtCount, runStatus } from "./montecarlo";
import { memo, table, type ViewResult } from "./view";
import { workers } from "./workers";

/** The tolerances of a sweep, as multiples of the finest. */
export const SWEEP = [16, 8, 4, 2, 1];
/** Warm-up samples on levels 0, 1, 2 before the first allocation. */
const N0 = 200;
/** Most samples a tolerance may ask of one level: 2·10⁶ is 32 MB of Q and Q_c, and minutes of work. */
const MAX_N = 2e6;
const INK = "rgba(207, 238, 255, 0.85)";

const sweeps = memo<MlmcSweep>();

export function renderMlmc(fig: Figure, st: State): ViewResult {
  const p = st.p, k = continuityK(st.continuity, p), ne0 = st.ne0, levels = st.levels;
  const supports = SUPPORTS[st.supports];
  const why = admissible({ p, k, ne: ne0, supports });
  const empty = () => fig.panels(1)[0].draw({ xlabel: "level ℓ", ylabel: "", series: [] });
  if (why) {
    empty();
    return { readout: `${continuityName(Math.max(k, 0))}, degree ${p}, ${ne0} elements: ${why}`, animate: false };
  }
  if (levels < 3) {
    empty();
    return { readout: "Multilevel Monte Carlo fits its rates on levels 1 … L and starts from levels 0, 1, 2: raise the number of levels to at least 3.", animate: false };
  }

  const kl = klOf(st), M = termsOf(kl);
  const spec: McSpec = { p, k, beam: beamCaseOf(st), qoi: st.qoi, field: fieldSpecOf(st), seed: st.seed };
  const neOf = (l: number) => ne0 * 2 ** l;
  const dofs = (l: number) => (l < 0 ? 0 : freeDofs({ p, k, ne: neOf(l), supports }));
  /** Cost of one sample of level ℓ: its fine and its coarse solve, in free coefficients solved for. */
  const cost = (l: number) => dofs(l) + dofs(l - 1);
  const key = JSON.stringify(["mlmc", spec, M, ne0, levels, st.mlSurvey, st.mlEps]);
  const run = workers.pool.open(key, () => ({ spec, kl }));
  const streamOf = (l: number): StreamSpec => ({ ne: neOf(l), coarse: l > 0, wantW: false, stream: l });
  for (let l = 0; l < levels; l++) run.stream(streamOf(l));
  // The scale ε is relative to: |E[Q₀]| over the survey's samples of level 0, once they are in.
  const level0 = run.streams[0].acc;
  const scale = level0.n >= st.mlSurvey ? Math.abs(level0.prefix(st.mlSurvey).q.mean) || 1 : NaN;
  const sweep = Number.isFinite(scale)
    ? sweeps(`${run.id}`, () => new MlmcSweep(levels, st.mlSurvey, SWEEP.map((f) => f * st.mlEps * scale),
      { N0: Math.min(N0, st.mlSurvey), Lmin: 2, cost }, MAX_N))
    : null;
  const demand = sweep ? sweep.step((l) => run.streams[l].acc) : run.streams.map(() => st.mlSurvey);
  demand.forEach((n, l) => workers.pool.demand(run, run.streams[l], n));

  // The survey, over the levels that have their samples (a prefix: level ℓ's check reads ℓ − 1).
  let complete = 0;
  while (complete < levels && run.streams[complete].acc.n >= st.mlSurvey) complete++;
  const S = survey(run.streams.slice(0, complete).map((s) => ({ fine: s.acc.values, coarse: s.acc.coarseValues })), st.mlSurvey);

  const ctx: Ctx = { st, scale, S, sweep, cost, dofs, neOf, run };
  if (!S.length) {
    empty();
    return { readout: readout(ctx, M), animate: run.running };
  }
  const [mean, variance, samples, costs] = fig.panels(4, { cols: 2 });
  drawMean(mean, ctx);
  drawVariance(variance, ctx);
  drawSamples(samples, ctx);
  drawCost(costs, ctx);
  return { readout: readout(ctx, M), animate: run.running };
}

interface Ctx {
  readonly st: State;
  readonly scale: number;
  readonly S: SurveyLevel[];
  readonly sweep: MlmcSweep | null;
  readonly cost: (l: number) => number;
  readonly dofs: (l: number) => number;
  readonly neOf: (l: number) => number;
  readonly run: McRun;
}

const levelAxis = (levels: number) => ({
  xlabel: "level ℓ  (ne = ne₀·2^ℓ)", xlim: [-0.3, levels - 0.7] as const,
  xfmt: (v: number) => (Number.isInteger(Math.round(v * 1e6) / 1e6) ? String(Math.round(v)) : ""),
});

/** A line through the fitted points at the fitted rate, 2^{−rate·ℓ}, over levels 1 … L. */
function rateLine(y: number[], rate: number, label: string, color: string): Series | null {
  const pts = y.map((v, l) => [l, v] as const).filter(([l, v]) => l >= 1 && v > 0 && Number.isFinite(v));
  if (pts.length < 2 || !Number.isFinite(rate)) return null;
  const mid = pts.reduce((s, [l, v]) => s + Math.log2(v) + rate * l, 0) / pts.length;
  const x = [pts[0][0], pts[pts.length - 1][0]];
  return { label, x, y: x.map((l) => 2 ** (mid - rate * l)), color, width: 1.25, dash: [2, 4], inert: true };
}

function bounds(series: Series[]): [number, number] {
  const ys = series.flatMap((s) => Array.from(s.y)).filter((v) => v > 0 && Number.isFinite(v));
  return ys.length ? [Math.min(...ys) / 4, Math.max(...ys) * 4] : [1e-6, 1];
}

function drawMean(plot: Plot, { st, scale, S }: Ctx): void {
  const l = S.map((v) => v.level);
  const q = S.map((v) => Math.abs(v.meanQ) / scale), y = S.map((v, i) => (i === 0 ? NaN : Math.abs(v.meanY) / scale));
  const half = S.map((v) => (Z95 * Math.sqrt(v.varY / v.N)) / scale);
  const alpha = -slope(S.map((v) => Math.abs(v.meanY)));
  const series: Series[] = [
    { label: "|E[Q_ℓ]|", x: l, y: q, color: SLOT[0], width: 2, markers: true },
    { label: `|E[Y_ℓ]| = |E[Q_ℓ − Q_ℓ₋₁]|`, x: l, y, color: SLOT[1], width: 2, markers: true },
    { label: "  ± 95% (sampling)", x: [], y: [], color: SLOT[1], width: 1, alpha: 0.6, inert: true },
  ];
  const fit = rateLine(y, alpha, `2^(−αℓ), α ≈ ${Number.isFinite(alpha) ? alpha.toFixed(2) : "—"}`, INK);
  if (fit) series.push(fit);
  const [lo, hi] = bounds(series);
  plot.draw({
    title: `mean against level — ${S.length ? S[0].N.toLocaleString() : 0} samples per level (survey)`,
    ...levelAxis(st.levels), ylabel: "relative to |E[Q]|", ylog: true, ylim: [lo, hi], legend: "bl", series,
    under: (ctx, a) => {
      ctx.strokeStyle = SLOT[1];
      ctx.globalAlpha = 0.6;
      ctx.lineWidth = 1.5;
      for (let i = 1; i < S.length; i++) {
        const top = y[i] + half[i], bot = y[i] - half[i];
        ctx.beginPath();
        ctx.moveTo(a.sx(i), a.sy(top));
        ctx.lineTo(a.sx(i), a.sy(bot > 0 ? bot : lo));
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
    },
    hover: (s, i) => `${s.label}, level ${l[i]}\n${s.y[i].toExponential(3)} (relative)` +
      (s.color === SLOT[1] ? `\n± ${half[i].toExponential(2)} (95%)` : ""),
  });
}

function drawVariance(plot: Plot, { st, scale, S }: Ctx): void {
  const l = S.map((v) => v.level), s2 = scale * scale;
  const q = S.map((v) => v.varQ / s2), y = S.map((v, i) => (i === 0 ? NaN : v.varY / s2));
  const beta = -slope(S.map((v) => v.varY));
  const series: Series[] = [
    { label: "V[Q_ℓ]", x: l, y: q, color: SLOT[0], width: 2, markers: true },
    { label: "V[Y_ℓ] = V[Q_ℓ − Q_ℓ₋₁]", x: l, y, color: SLOT[1], width: 2, markers: true },
  ];
  const fit = rateLine(y, beta, `2^(−βℓ), β ≈ ${Number.isFinite(beta) ? beta.toFixed(2) : "—"}`, INK);
  if (fit) series.push(fit);
  const [lo, hi] = bounds(series);
  plot.draw({
    title: "variance against level: what each level's estimate costs in samples",
    ...levelAxis(st.levels), ylabel: "relative to E[Q]²", ylog: true, ylim: [lo, hi], legend: "bl", series,
    hover: (s, i) => `${s.label}, level ${l[i]}\n${s.y[i].toExponential(3)} (relative)` +
      (i > 0 && s.color === SLOT[1] ? `\nV[Y]/V[Q] = ${(S[i].varY / S[i].varQ).toExponential(2)}\nkurtosis ${S[i].kurtosis.toFixed(1)}` : ""),
  });
}

function drawSamples(plot: Plot, { st, sweep }: Ctx): void {
  const series: Series[] = (sweep?.runs ?? []).map((a, i) => {
    const N = a.status === "converged" || a.status === "failed" ? a.N : a.wants();
    return {
      label: `ε = ${(SWEEP[i] * st.mlEps).toExponential(1)}${a.done ? (a.status === "converged" ? "" : ` (${a.status})`) : " (running)"}`,
      x: N.map((_, l) => l), y: N.slice(), color: SLOT[i % SLOT.length], width: 2, markers: true,
      dash: a.done ? undefined : [5, 4], alpha: a.done ? 1 : 0.6,
    };
  });
  const [lo, hi] = bounds(series);
  plot.draw({
    title: "samples per level, N_ℓ ∝ √(V_ℓ / C_ℓ), for each tolerance",
    ...levelAxis(st.levels), ylabel: "N_ℓ", ylog: true, ylim: [Math.max(1, lo), hi], legend: "tr", series,
    hover: (s, i) => `${s.label}\nlevel ${i}: N = ${s.y[i].toLocaleString()}`,
  });
}

/** Plain Monte Carlo to the same ε on level L: V[Q_L] / ((1 − θ)ε²) samples of one fine solve each. */
function mcCost(c: Ctx, a: Mlmc): number {
  const { S } = c, L = Math.min(a.L, S.length - 1);
  return L < 0 ? NaN : (S[L].varQ * c.dofs(a.L)) / ((1 - a.theta) * a.opts.eps ** 2);
}

function drawCost(plot: Plot, c: Ctx): void {
  const { st, sweep } = c, done = (sweep?.runs ?? []).map((a, i) => ({ a, e: SWEEP[i] * st.mlEps })).filter(({ a }) => a.status === "converged");
  const e = done.map(({ e }) => e);
  const ml = done.map(({ a, e }) => e * e * a.cost), mc = done.map(({ a, e }) => e * e * mcCost(c, a));
  const series: Series[] = [
    { label: "multilevel Monte Carlo", x: e, y: ml, color: SLOT[1], width: 2.5, markers: true },
    { label: "plain Monte Carlo on level L", x: e, y: mc, color: SLOT[0], width: 2.5, markers: true },
  ];
  const all = SWEEP.map((f) => f * st.mlEps), [lo, hi] = bounds(series);
  plot.draw({
    title: done.length ? "cost to reach a root-mean-square error ε, times ε²" : "cost against ε — waiting for the first tolerance to converge",
    xlabel: "ε relative to |E[Q]|", ylabel: "ε² × cost  [free coefficients solved for]", xlog: true, ylog: true,
    xlim: [Math.min(...all) / 1.5, Math.max(...all) * 1.5], ylim: [lo, hi], legend: "tr", series,
    xfmt: decade,
    hover: (s, i) => {
      const { a } = done[i];
      return `${s.label}, ε = ${e[i].toExponential(1)}\nL = ${a.L}\ncost ${fmtCount(s.y[i] / (e[i] * e[i]))}\n` +
        `ε²·cost = ${s.y[i].toPrecision(3)}`;
    },
  });
}

function readout(c: Ctx, M: number): string {
  const { st, scale, S, sweep, run } = c, d = displayOf(st);
  const v = (x: number) => valueOf(d, st.qoi, x, st.load);
  const levels = st.levels;
  const lines = [
    `levels ℓ = 0 … ${levels - 1}: ne = ${st.ne0}·2^ℓ (to ${c.neOf(levels - 1)}), p = ${st.p} ${continuityName(continuityK(st.continuity, st.p))}, ` +
      `${st.supports}; Q = ${QOI_NAME[st.qoi]}`,
    `input: ${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L, σ = ${plain(st.sigma, 3)}, M = ${M}` +
      `${st.massFollows ? ", mass follows depth" : ""}${st.loadSigma > 0 ? `, load σ_q = ${plain(st.loadSigma, 3)}` : ""}; seed ${st.seed}; ` +
      `level ℓ draws from random stream ℓ`,
    "",
  ];
  const n = run.streams.map((s) => s.acc.n), total = n.reduce((a, b) => a + b, 0), wall = run.wallMs / 1000;
  lines.push(`${runStatus(run)}: ${total.toLocaleString()} samples (${n.map(fmtCount).join(" · ")} by level) on ${workers.size} workers` +
    (wall > 0 ? `, ${plain(wall, 3)} s` : ""));
  if (!S.length) return [...lines, "", "the survey is waiting for its first level…"].join("\n");
  if (!sweep) return lines.join("\n");

  // ---- the survey ----
  const ms = run.streams.map((s) => s.msPerSample);
  const rows = S.map((s, l) => [
    String(l), String(c.neOf(l)), v(s.meanQ),
    l === 0 ? "—" : (Math.abs(s.meanY) / scale).toExponential(2),
    l === 0 ? "—" : (Z95 * Math.sqrt(s.varY / s.N) / scale).toExponential(1),
    (s.varY / scale ** 2).toExponential(2),
    l === 0 ? "—" : (s.varY / s.varQ).toExponential(2),
    Number.isFinite(s.kurtosis) ? s.kurtosis.toFixed(1) : "—",
    Number.isFinite(s.consistency) ? s.consistency.toFixed(2) : "—",
    String(c.cost(l)),
    Number.isFinite(ms[l]) ? ms[l].toFixed(3) : "—",
  ]);
  const alpha = -slope(S.map((s) => Math.abs(s.meanY))), beta = -slope(S.map((s) => s.varY));
  const gamma = slope(S.map((_, l) => c.cost(l))), gammaMs = slope(ms.slice(0, S.length));
  lines.push(
    "",
    `survey: ${st.mlSurvey.toLocaleString()} samples on each level; Y₀ = Q₀, Y_ℓ = Q_ℓ − Q_ℓ₋₁ (same ω); relative to |E[Q]| ≈ |E[Q₀]| = ${v(scale)}`,
    table(["ℓ", "ne", "E[Q_ℓ]", "|E[Y_ℓ]|", "± 95%", "V[Y_ℓ]", "V[Y]/V[Q]", "kurt", "check", "C_ℓ", "ms"], rows),
    `α ≈ ${fmt2(alpha)} (bias),  β ≈ ${fmt2(beta)} (variance),  γ = ${fmt2(gamma)} (cost model: free coefficients of both solves; ` +
      `${Number.isFinite(gammaMs) ? `measured ${fmt2(gammaMs)} in worker time` : "measured — "})`,
    regime(beta, gamma),
  );
  const bad = S.filter((s) => s.consistency > 1);
  if (bad.length) lines.push(`check > 1 on level${bad.length > 1 ? "s" : ""} ${bad.map((s) => s.level).join(", ")}: E[Y_ℓ] ≠ E[Q_ℓ] − E[Q_ℓ₋₁] beyond three standard errors (expected now and then by chance; always, if the coupling were broken)`);
  const kurt = S.filter((s) => s.level > 0 && s.kurtosis > 100);
  if (kurt.length) lines.push(`kurtosis > 100 on level${kurt.length > 1 ? "s" : ""} ${kurt.map((s) => s.level).join(", ")}: V[Y] there rests on a few rare samples, and the allocation trusts it less than it seems`);
  if (S.length > 1 && S[1].varY / S[1].varQ > 0.1)
    lines.push(`V[Y₁]/V[Q₁] = ${(S[1].varY / S[1].varQ).toFixed(2)}: the coarsest mesh (${st.ne0} elements, ${plain(st.ell * st.ne0, 2)} per correlation length) barely sees the field, so its corrections are hardly smaller than Q`);

  // ---- the tolerances ----
  const trows = sweep.runs.map((a, i) => {
    const e = SWEEP[i] * st.mlEps;
    const est = a.N.length ? `${v(a.estimate)} ± ${(Z95 * Math.sqrt(a.sampleVariance) / scale).toExponential(1)}` : "—";
    const mc = a.status === "converged" ? mcCost(c, a) : NaN;
    return [
      e.toExponential(1), a.status, a.N.length ? String(a.L) : "—", est,
      Number.isFinite(a.bias) ? (a.bias / scale).toExponential(1) : "—",
      a.N.length ? a.N.map(fmtCount).join(" ") : "—",
      a.N.length ? fmtCount(a.cost) : "—",
      Number.isFinite(mc) ? fmtCount(mc) : "—",
      Number.isFinite(mc) ? `${(mc / a.cost).toFixed(1)}×` : "—",
    ];
  });
  lines.push(
    "",
    `adaptive MLMC (Giles): variance (1 − θ)ε², bias² θε², θ = ${THETA}; ε relative to |E[Q]|; ± is 95% of the sampling error alone`,
    table(["ε", "status", "L", "estimate ± 95%", "bias est.", "N_ℓ", "cost", "MC cost", "saving"], trows),
  );
  if (sweep.runs.some((a) => a.status === "failed"))
    lines.push(`failed: the bias test still asked for a finer level than ℓ = ${levels - 1}; raise the number of levels, or the tolerance`);
  const over = sweep.runs.filter((a) => a.status === "over budget");
  if (over.length) {
    const a = over[over.length - 1];
    lines.push(`over budget: ε = ${(SWEEP[sweep.runs.indexOf(a)] * st.mlEps).toExponential(1)} would need ${fmtCount(Math.max(...a.wants()))} samples on one level ` +
      `(the limit is ${fmtCount(MAX_N)}) — V[Q₀]/E[Q]² = ${(S[0].varQ / scale ** 2).toPrecision(3)} sets that, not the mesh`);
  }
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}

const fmt2 = (x: number) => (Number.isFinite(x) ? x.toFixed(2) : "—");

/** Giles' complexity theorem, read off the measured rates. */
function regime(beta: number, gamma: number): string {
  if (!Number.isFinite(beta) || !Number.isFinite(gamma)) return "";
  if (beta > gamma * 1.05)
    return "β > γ: the variance falls faster than the cost rises, so the coarsest levels carry the work and MLMC costs O(ε⁻²) — as if the mesh were free.";
  if (beta > gamma * 0.95)
    return "β ≈ γ: every level costs about the same, and MLMC costs O(ε⁻² (log ε)²).";
  return "β < γ: the finest levels dominate, and MLMC costs O(ε^(−2−(γ−β)/α)) — still below plain Monte Carlo's ε^(−2−γ/α).";
}
