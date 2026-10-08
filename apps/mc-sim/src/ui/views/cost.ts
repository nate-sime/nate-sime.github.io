// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Plain Monte Carlo against multilevel Monte Carlo, at the same tolerance:
 * how many samples each takes, on which meshes, and what that costs — the
 * comparison the live view (`live.ts`) shows when "compare cost with MC" is
 * pressed. Not a view of its own: these are its pictures and its table.
 *
 * It reads the MLMC run — the same survey, the same sweep of five
 * tolerances, the same samples, kept in the pool — for one tolerance ε
 * chosen in the pane. Both estimators are held to the same budget: bias from
 * the finest level L the sweep needed for ε, and sampling variance
 * (1 − θ)ε². Multilevel Monte Carlo's N_ℓ are the ones its adaptive
 * algorithm settled on. Plain Monte Carlo's count is what it would need on
 * level L for the same variance,
 *
 *   N = V[Q_L] / ((1 − θ) ε²),
 *
 * from the survey's measured V[Q_L] — computed, not run: on a fine level at a
 * fine tolerance, running it is the cost being compared. Costs are in the
 * cost model's unit (free coefficients for a beam, dofs × bandwidth² for a
 * plate), a plain-MC sample costing one fine solve, C_L^MC = dofs_L, and a
 * multilevel one its fine and coarse solve, C_ℓ = dofs_ℓ + dofs_ℓ₋₁.
 *
 * Four pictures:
 * - samples per level at ε: MLMC's N_ℓ on every level beside plain MC's N, all
 *   on level L — most of MLMC's samples sit where a sample is cheapest;
 * - cost to reach ε, on a linear axis so the ratio reads as length: MLMC's bar
 *   split into its levels' shares, plain MC's one solid bar;
 * - samples against ε for every tolerance: MLMC may take more samples in all
 *   than plain MC, and still cost a fraction — the samples are cheap;
 * - ε² × cost against ε, flat for MLMC when β > γ (`drawCost`, from the
 *   MLMC view).
 */

import type { Mlmc } from "../../mc/mlmc";
import { SLOT, type Axes, type Plot, type Series } from "../plot";
import { plain } from "../units";
import { SWEEP, mcCost, type MlmcSession } from "./mlmc";
import { fmtCount } from "./montecarlo";
import { workUnit } from "./structure";
import { table } from "./view";

export const MLMC_INK = SLOT[1];
const MC_INK = SLOT[0];
const INK = "rgba(207, 238, 255, 0.92)";
const SURFACE_GAP = "#05050c";

/** One tolerance of the sweep, both ways. */
export interface Pair {
  readonly eps: number;
  readonly alg: Mlmc;
  readonly converged: boolean;
  /** Said beside a tolerance that has not converged. */
  readonly note: string;
  /** MLMC's samples per level 0 … L, and their cost per sample. */
  readonly N: readonly number[];
  readonly C: readonly number[];
  readonly mlCost: number;
  /** Plain MC on level L: samples, cost per sample, total. */
  readonly L: number;
  readonly nMC: number;
  readonly cMC: number;
  readonly mcCost: number;
}

/** Tolerance `k` of the sweep (0 the coarsest), MLMC's numbers and plain MC's beside them; null until it has a level structure the survey covers. */
export function pairOf(c: MlmcSession, k: number): Pair | null {
  const alg = c.sweep?.runs[k];
  if (!alg) return null;
  const converged = alg.status === "converged";
  const N = (converged ? alg.N : alg.wants()).slice();
  const L = N.length - 1;
  if (L >= c.S.length) return null;
  const C = N.map((_, l) => c.cost(l));
  const mlCost = N.reduce((s, n, l) => s + n * C[l], 0);
  const cMC = c.dofs(L);
  const total = mcCost(c, alg);
  const note = converged ? "" : alg.status === "sampling" ? " (MLMC still settling)"
    : alg.status === "over budget" ? " (over budget: what MLMC would need)" : " (MLMC wanted a finer level than the run has)";
  return { eps: SWEEP[k] * c.st.mlEps, alg, converged, note, N, C, mlCost, L, nMC: Math.ceil(total / cMC), cMC, mcCost: total };
}

/** A bar from the bottom of the axes up to `top`, rounded at its data end. */
export function bar(ctx: CanvasRenderingContext2D, a: Axes, x: number, halfWidth: number, top: number, color: string): void {
  const X0 = a.sx(x - halfWidth), X1 = a.sx(x + halfWidth), Y = Math.max(a.box.t, a.sy(top)), B = a.box.b;
  if (!(B - Y > 0)) return;
  ctx.fillStyle = color;
  ctx.beginPath();
  ctx.roundRect(X0, Y, X1 - X0, B - Y, [4, 4, 0, 0]);
  ctx.fill();
}

export function label(ctx: CanvasRenderingContext2D, X: number, Y: number, text: string, align: CanvasTextAlign = "center"): void {
  ctx.fillStyle = INK;
  ctx.font = "11px ui-monospace, monospace";
  ctx.textAlign = align;
  ctx.textBaseline = "bottom";
  ctx.fillText(text, X, Y);
}

/**
 * MLMC's samples on each level beside plain MC's on level L. The live view
 * passes the level whose pair it is showing, which is outlined — or, past the
 * L this tolerance needs, marked where its bar would be — and its own title.
 */
export function drawPerLevel(plot: Plot, p: Pair, opts: { title?: string; highlight?: number } = {}): void {
  const hl = opts.highlight;
  const levels = Math.max(p.L + 1, hl === undefined ? 0 : hl + 1, 1), w = 0.17;
  // Headroom on the log axis, sized in pixels: the legend's two rows and the tallest bar's count need
  // ~62 px at the top, whatever the panel's height — and the panel is short when it shares the live view.
  const tallest = Math.max(p.nMC, ...p.N), lo = Math.max(1, Math.min(p.nMC, ...p.N) / 4);
  const boxH = Math.max(60, plot.canvas.clientHeight - 30 - 46), f = Math.min(0.7, 62 / boxH);
  const top = 10 ** ((Math.log10(tallest) - f * Math.log10(lo)) / (1 - f));
  const hidden = "rgba(0, 0, 0, 0)";
  plot.draw({
    title: opts.title ?? `samples per level to reach ε = ${p.eps.toExponential(1)}${p.note}`,
    xlabel: "level ℓ  (ne = ne₀·2^ℓ)", ylabel: "samples", ylog: true,
    xlim: [-0.6, levels - 0.4], ylim: [lo, top],
    xfmt: (v) => (Number.isInteger(Math.round(v * 1e6) / 1e6) && v >= 0 && v < levels ? String(Math.round(v)) : ""),
    legend: "tr",
    // The bars are drawn below; the series carry the legend swatches and the hover points at the bars' tops.
    series: [
      { label: "multilevel MC: N_ℓ on each level", x: [], y: [], color: MLMC_INK, width: 8, inert: true },
      { label: `plain MC: all N on level ${p.L}`, x: [], y: [], color: MC_INK, width: 8, inert: true },
      { label: "multilevel MC", x: p.N.map((_, l) => l - w), y: p.N.slice(), color: hidden, width: 0.001, unlisted: true },
      { label: "plain MC", x: [p.L + w], y: [p.nMC], color: hidden, width: 0.001, unlisted: true },
    ],
    under: (ctx, a) => {
      p.N.forEach((n, l) => bar(ctx, a, l - w, w * 0.92, n, MLMC_INK));
      bar(ctx, a, p.L + w, w * 0.92, p.nMC, MC_INK);
    },
    over: (ctx, a) => {
      p.N.forEach((n, l) => label(ctx, a.sx(l - w), a.sy(n) - 3, fmtCount(n)));
      label(ctx, a.sx(p.L + w), a.sy(p.nMC) - 3, fmtCount(p.nMC));
      if (hl === undefined) return;
      ctx.strokeStyle = INK;
      ctx.lineWidth = 2;
      if (hl <= p.L) {
        // An outline 3px outside the bar (centre hl − w, half-width 0.92w), leaving the bar's own colour untouched.
        const X0 = a.sx(hl - w - 0.92 * w) - 3, X1 = a.sx(hl - w + 0.92 * w) + 3;
        const Y = Math.max(a.box.t, a.sy(p.N[hl])) - 3;
        ctx.beginPath();
        ctx.roundRect(X0, Y, X1 - X0, a.box.b - Y, [6, 6, 0, 0]);
        ctx.stroke();
        label(ctx, a.sx(hl - w), a.box.b - 4, "shown");
      } else {
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(a.sx(hl), a.box.t);
        ctx.lineTo(a.sx(hl), a.box.b);
        ctx.stroke();
        ctx.setLineDash([]);
        // Up the line itself, rotated: beside it the bars and the plot's edge leave no room.
        ctx.save();
        ctx.translate(a.sx(hl) - 5, (a.box.t + a.box.b) / 2 + 20);
        ctx.rotate(-Math.PI / 2);
        label(ctx, 0, 0, `shown ℓ = ${hl}: past L = ${p.L}`);
        ctx.restore();
      }
    },
    hover: (s, i) => s.label === "plain MC"
      ? `plain MC on level ${p.L}\nN = ${fmtCount(p.nMC)}\n${fmtCount(p.cMC)} per sample (one fine solve)\ncost ${fmtCount(p.mcCost)}`
      : `MLMC, level ${i}\nN_ℓ = ${fmtCount(p.N[i])}\n${fmtCount(p.C[i])} per sample (fine${i ? " + coarse" : ""} solve)\ncost ${fmtCount(p.N[i] * p.C[i])} (${(100 * p.N[i] * p.C[i] / p.mlCost).toFixed(1)}%)`,
  });
}

export function drawTotal(plot: Plot, c: MlmcSession, p: Pair): void {
  const unit = workUnit(c.st), max = Math.max(p.mcCost, p.mlCost);
  const ratio = p.mcCost / p.mlCost;
  plot.draw({
    title: `cost to reach ε = ${p.eps.toExponential(1)} — plain MC ${ratio >= 1 ? `${plain(ratio, 3)}× multilevel MC's` : `${plain(1 / ratio, 3)}× cheaper`}`,
    xlabel: `cost  [${unit}]`, ylabel: "", xlim: [0, max * 1.12], ylim: [0.3, 2.7],
    yfmt: (v) => (Math.abs(v - 1) < 1e-9 ? "plain MC" : Math.abs(v - 2) < 1e-9 ? "MLMC" : ""),
    xfmt: (v) => fmtCount(v),
    series: [
      { label: "hover", x: [p.mcCost], y: [1], color: "rgba(0,0,0,0)", width: 0.001, unlisted: true },
      { label: "hover", x: p.N.map((_, l) => p.N.slice(0, l + 1).reduce((s, n, j) => s + n * p.C[j], 0)), y: p.N.map(() => 2), color: "rgba(0,0,0,0)", width: 0.001, unlisted: true },
    ],
    under: (ctx, a) => {
      const h = 0.32;
      const rect = (x0: number, x1: number, y: number, color: string, round: boolean) => {
        const X0 = a.sx(x0), X1 = a.sx(x1), Y0 = a.sy(y + h), Y1 = a.sy(y - h);
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.roundRect(X0, Y0, Math.max(1, X1 - X0), Y1 - Y0, round ? [0, 4, 4, 0] : 0);
        ctx.fill();
      };
      rect(0, p.mcCost, 1, MC_INK, true);
      // MLMC's bar in its levels' shares, a 2px surface gap between them.
      let x = 0;
      p.N.forEach((n, l) => {
        const share = n * p.C[l];
        rect(x, x + share, 2, MLMC_INK, l === p.L);
        if (l > 0) {
          ctx.fillStyle = SURFACE_GAP;
          ctx.fillRect(a.sx(x) - 1, a.sy(2 + h), 2, a.sy(2 - h) - a.sy(2 + h));
        }
        x += share;
      });
    },
    over: (ctx, a) => {
      label(ctx, a.sx(p.mcCost) + 6, a.sy(1) + 5, fmtCount(p.mcCost), "left");
      label(ctx, a.sx(p.mlCost) + 6, a.sy(2) + 5, fmtCount(p.mlCost), "left");
      // Each level's share named inside it, where it fits.
      let x = 0;
      p.N.forEach((n, l) => {
        const share = n * p.C[l], X0 = a.sx(x), X1 = a.sx(x + share);
        if (X1 - X0 > 34) {
          ctx.fillStyle = "#05050c";
          ctx.font = "11px ui-monospace, monospace";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillText(`ℓ=${l}`, (X0 + X1) / 2, a.sy(2));
        }
        x += share;
      });
    },
    hover: (s, i) => s.y[i] === 1
      ? `plain MC on level ${p.L}\n${fmtCount(p.nMC)} samples × ${fmtCount(p.cMC)}\n= ${fmtCount(p.mcCost)}`
      : `MLMC through level ${i}\nlevel ${i}: ${fmtCount(p.N[i])} × ${fmtCount(p.C[i])} = ${fmtCount(p.N[i] * p.C[i])}\nrunning total ${fmtCount(s.x[i])}`,
  });
}

export function drawSamplesVsEps(plot: Plot, c: MlmcSession, pairs: readonly (Pair | null)[]): void {
  const done = pairs.filter((p): p is Pair => p !== null && p.converged);
  const e = done.map((p) => p.eps);
  const series: Series[] = [
    { label: "plain MC: N, all on level L", x: e, y: done.map((p) => p.nMC), color: MC_INK, width: 2.5, markers: true },
    { label: "multilevel MC: Σ N_ℓ, all levels", x: e, y: done.map((p) => p.N.reduce((s, n) => s + n, 0)), color: MLMC_INK, width: 2.5, markers: true },
    { label: "multilevel MC: N_L, finest level only", x: e, y: done.map((p) => p.N[p.L]), color: MLMC_INK, width: 1.5, dash: [5, 4], markers: true },
  ];
  const all = SWEEP.map((f) => f * c.st.mlEps);
  const ys = series.flatMap((s) => Array.from(s.y)).filter((v) => v > 0);
  plot.draw({
    title: done.length ? "samples against ε: more in all for MLMC, few on the fine mesh" : "samples against ε — waiting for the first tolerance to converge",
    xlabel: "ε relative to |E[Q]|", ylabel: "samples", xlog: true, ylog: true, legend: "tr",
    xlim: [Math.min(...all) / 1.5, Math.max(...all) * 1.5],
    ylim: ys.length ? [Math.min(...ys) / 4, Math.max(...ys) * 4] : [1, 1e6],
    series,
    hover: (s, i) => `${s.label}\nε = ${e[i].toExponential(1)}, L = ${done[i].L}\n${fmtCount(s.y[i])} samples`,
  });
}

/** The comparison as readout lines: the chosen tolerance level by level, then every tolerance that has converged. */
export function costLines(c: MlmcSession, k: number): string[] {
  const pairs = SWEEP.map((_, i) => pairOf(c, i));
  const lines: string[] = [];
  const p = pairs[k];
  if (p) {
    const sumN = p.N.reduce((s, n) => s + n, 0), ratio = p.mcCost / p.mlCost;
    const share0 = (p.N[0] * p.C[0]) / p.mlCost;
    lines.push(
      `cost against plain MC at ε = ${p.eps.toExponential(1)} (relative to |E[Q]|), both with bias from level L = ${p.L} and sampling variance (1 − θ)ε²:`,
      table(["", "samples", "on", "per sample", "cost", "share"], [
        ...p.N.map((n, l) => [`MLMC, level ${l}`, fmtCount(n), `${c.neOf(l)}${l ? ` + ${c.neOf(l - 1)}` : ""}`, fmtCount(p.C[l]), fmtCount(n * p.C[l]), `${(100 * n * p.C[l] / p.mlCost).toFixed(1)}%`]),
        ["MLMC, in all", fmtCount(sumN), "", "", fmtCount(p.mlCost), "100%"],
        [`plain MC, level ${p.L}`, fmtCount(p.nMC), String(c.neOf(p.L)), fmtCount(p.cMC), fmtCount(p.mcCost), ""],
      ]),
      `plain Monte Carlo costs ${plain(ratio, 3)}× as much${p.note}. ` +
        `MLMC takes ${sumN >= p.nMC ? `${plain(sumN / p.nMC, 3)}× as many samples in all` : `${plain(p.nMC / sumN, 3)}× fewer samples`}, ` +
        `but ${(100 * p.N[0] / sumN).toFixed(1)}% of them on level 0, where one costs ${plain(p.cMC / p.C[0], 3)}× less than a level-${p.L} solve; ` +
        `level 0 carries ${(100 * share0).toFixed(0)}% of its cost.`,
      `plain MC's N = V[Q_L]/((1 − θ)ε²) uses the survey's V[Q_${p.L}]: the count it would need, not a run — on a fine level at a fine tolerance, running it is the cost in question.`,
    );
  } else lines.push(c.sweep ? "cost against plain MC: this tolerance has not yet settled on its levels." : "cost against plain MC: waiting for the survey…");

  const done = pairs.filter((q): q is Pair => q !== null && q.converged);
  if (done.length) {
    lines.push("", `every tolerance (costs in ${workUnit(c.st)}):`, table(
      ["ε", "L", "MLMC N_ℓ", "Σ N_ℓ", "MLMC cost", "MC N", "MC cost", "MC / MLMC"],
      done.map((q) => [
        q.eps.toExponential(1), String(q.L), q.N.map(fmtCount).join(" "), fmtCount(q.N.reduce((s, n) => s + n, 0)),
        fmtCount(q.mlCost), fmtCount(q.nMC), fmtCount(q.mcCost), `${plain(q.mcCost / q.mlCost, 3)}×`,
      ]),
    ));
    const r = done.map((q) => q.mcCost / q.mlCost);
    const deeper = done[done.length - 1].L > done[0].L;
    if (r.length > 1 && r[r.length - 1] > r[0])
      lines.push(`the gap widens as ε falls, ${plain(r[0], 2)}× to ${plain(r[r.length - 1], 2)}× here: plain MC pays ε⁻² samples at the finest level's price` +
        `${deeper ? ", and that level moves finer as ε falls" : ""}; MLMC pays its ε⁻² mostly at the coarsest level's, with a fixed warm-up of samples on the fine levels that matters less and less.`);
  }
  return lines;
}
