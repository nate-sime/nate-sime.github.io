// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 5 on screen: plain Monte Carlo on one level of the hierarchy.
 *
 * Samples run in the worker pool while the view redraws; every sample is also
 * solved on the parent level with the same ω. Four pictures of one run, side
 * by side:
 *
 * - histogram: the distribution of Q_ℓ, a normal of the same mean and spread,
 *   and two lines — the estimate of E[Q], and Q(E[inputs]), the answer for the
 *   mean beam. They differ (Jensen: Q is not linear in the stiffness), which is
 *   why a random input cannot be replaced by its mean.
 * - running mean: the estimate against N with its 95% CLT interval, closing
 *   as σ/√N.
 * - error vs N: the root-mean-square error split as MSE = bias² + σ²/N. The
 *   sampling part falls with N; the bias — indicated by |E[Q_ℓ − Q_{ℓ−1}]|, which
 *   bounds it for any convergence rate α ≥ 1 — does not. Past their crossing
 *   N* = σ²/bias², more samples on this mesh buy nothing: refine the mesh, or
 *   share the work across meshes, which is multilevel Monte Carlo.
 * - deflection bands: the mean field and its ±σ, ±2σ bands, the first few
 *   sample paths, and the deflection of the mean beam.
 */

import { SUPPORTS, admissible } from "../../beam/beam";
import { qoiPoint } from "../../hierarchy";
import type { McRun } from "../../mc/pool";
import { Sampler, fieldX, type McSpec } from "../../mc/sampler";
import { Z95, histogram } from "../../mc/stats";
import { termsOf, type KLData } from "../../random/kl";
import type { Figure } from "../figure";
import { SLOT, type Axes, type Plot, type Series } from "../plot";
import { beamCaseOf, continuityK, continuityName, displayOf, fieldSpecOf, type State } from "../state";
import { axisUnit, deflectionScale, plain, qoiScale, referenceNote } from "../units";
import { QOI_NAME, valueOf } from "./convergence";
import { KERNEL_NAME, fill, klOf } from "./field";
import { memo, type ViewResult } from "./view";
import { workers } from "./workers";

const INK = "rgba(207, 238, 255, 0.85)";
const BAR = "rgba(57, 135, 229, 0.55)";

const meanField = memo<{ Q: number; w: Float64Array }>();

export function renderMonteCarlo(fig: Figure, st: State): ViewResult {
  const p = st.p, k = continuityK(st.continuity, p);
  const ne = st.ne0 * 2 ** st.mcLevel, coarse = st.mcLevel > 0;
  const supports = SUPPORTS[st.supports];
  const why = admissible({ p, k, ne: coarse ? ne / 2 : ne, supports });
  if (why) {
    fig.panels(1)[0].draw({ xlabel: "Q", ylabel: "density", series: [] });
    return { readout: `${continuityName(Math.max(k, 0))}, degree ${p}, ${coarse ? ne / 2 : ne} elements: ${why}`, animate: false };
  }

  const kl = klOf(st), M = termsOf(kl);
  const bc = beamCaseOf(st);
  const spec: McSpec = { p, k, beam: bc, qoi: st.qoi, field: fieldSpecOf(st), seed: st.seed };
  const key = JSON.stringify(["mc", spec, M, ne, coarse]);
  const run = workers.pool.ensure(key, () => ({ spec, kl }), { ne, coarse, wantW: true, stream: 0 }, st.mcSamples);
  const mf = meanField(key, () =>
    new Sampler({ ...spec, field: { ...spec.field, sigma: 0, loadSigma: 0 } }, kl).solve(0, ne, true) as { Q: number; w: Float64Array });

  const d = displayOf(st), acc = run.acc;
  const status = runStatus(run);
  const head = `level ℓ = ${st.mcLevel}: ${ne} elements${coarse ? ` (parent ${ne / 2}, same ω)` : ""}, p = ${p} ${continuityName(k)}, ` +
    `${st.supports}; Q = ${QOI_NAME[st.qoi]}\n` +
    `input: ${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L, σ = ${plain(st.sigma, 3)}, M = ${M}` +
    `${st.massFollows ? ", mass follows depth" : ""}${st.loadSigma > 0 ? `, load σ_q = ${plain(st.loadSigma, 3)}` : ""}; seed ${st.seed}`;
  if (acc.n < 2) {
    fig.panels(1)[0].draw({ xlabel: "Q", ylabel: "density", series: [] });
    return { readout: `${head}\n\n${status}: waiting for the first samples…`, animate: run.running };
  }

  const sc = qoiScale(d, st.qoi, st.load);
  const u = d.dimensional ? axisUnit(Math.max(Math.abs(acc.min), Math.abs(acc.max)) * sc.factor, sc.unit) : { factor: 1, label: "" };
  const kq = sc.factor * u.factor;
  const qlabel = `Q${u.label ? ` [${u.label}]` : ""}`;

  const [hist, running, error, bands] = fig.panels(4, { cols: 2 });
  drawHistogram(hist, run, mf.Q, kq, qlabel);
  drawRunning(running, run, mf.Q, kq, qlabel);
  drawError(error, run, coarse);
  drawBands(bands, st, run, mf.w);

  return { readout: readout(st, run, mf.Q, head, status, coarse, kl), animate: run.running };
}

function drawHistogram(plot: Plot, run: McRun, Qmf: number, k: number, xlabel: string): void {
  const acc = run.acc, { mean, sd } = acc.q;
  const h = histogram(acc.values, acc.min, acc.max, sd);
  const edges = Array.from(h.edges, (e) => e * k), dens = Array.from(h.density, (v) => v / k);
  const pad = 0.04 * (edges[edges.length - 1] - edges[0] || Math.abs(edges[0]) || 1);
  const lo = Math.min(edges[0], Qmf * k) - pad, hi = Math.max(edges[edges.length - 1], Qmf * k) + pad;
  const xs = Array.from({ length: 161 }, (_, i) => lo + ((hi - lo) * i) / 160);
  const s = sd * k, m = mean * k;
  const pdf = xs.map((x) => Math.exp(-0.5 * ((x - m) / s) ** 2) / (s * Math.sqrt(2 * Math.PI)));
  const top = Math.max(...dens, ...pdf.filter(Number.isFinite)) * 1.15;
  const centers = dens.map((_, i) => 0.5 * (edges[i] + edges[i + 1]));
  plot.draw({
    title: `distribution of Q_ℓ — ${acc.n.toLocaleString()} samples`,
    xlabel, ylabel: "probability density", xlim: [lo, hi], ylim: [0, top],
    series: [
      { label: "samples (histogram)", x: [], y: [], color: BAR, width: 8, inert: true },
      { label: "bin", x: centers, y: dens, color: "rgba(0, 0, 0, 0)", width: 0.001, unlisted: true },
      { label: "normal, same mean and σ", x: xs, y: pdf, color: SLOT[1], width: 1.5, dash: [6, 4], inert: true },
      { label: "E[Q] estimate", x: [], y: [], color: INK, width: 2, inert: true },
      { label: "Q(E[inputs]), the mean beam", x: [], y: [], color: SLOT[2], width: 2, dash: [3, 3], inert: true },
    ],
    under: (ctx, a) => {
      ctx.fillStyle = BAR;
      for (let i = 0; i < dens.length; i++) {
        const x0 = a.sx(edges[i]), x1 = a.sx(edges[i + 1]), y = a.sy(dens[i]);
        ctx.fillRect(x0 + 0.5, y, Math.max(1, x1 - x0 - 1), a.sy(0) - y);
      }
    },
    over: (ctx, a) => {
      vline(ctx, a, m, INK, []);
      vline(ctx, a, Qmf * k, SLOT[2], [3, 3]);
    },
    hover: (_, i) => `bin [${edges[i].toPrecision(4)}, ${edges[i + 1].toPrecision(4)}]\n` +
      `${Math.round(dens[i] * (edges[i + 1] - edges[i]) * acc.n)} samples`,
  });
}

function drawRunning(plot: Plot, run: McRun, Qmf: number, k: number, ylabel: string): void {
  const t = run.acc.trajectory;
  const n = t.map((c) => c.n), m = t.map((c) => c.mean * k);
  const half = t.map((c) => Z95 * (c.sd / Math.sqrt(c.n)) * k);
  const lo = m.map((v, i) => v - half[i]), hi = m.map((v, i) => v + half[i]);
  // Scale to the interval from N ≈ 20 on; the first few are as wide as Q itself.
  const from = Math.max(0, n.findIndex((v) => v >= 20));
  const ys = [...lo.slice(from), ...hi.slice(from), Qmf * k].filter(Number.isFinite);
  const span = Math.max(...ys) - Math.min(...ys) || Math.abs(m[m.length - 1]) * 1e-3;
  const last = m[m.length - 1];
  plot.draw({
    title: "running mean with its 95% interval, Q̄_N ± 1.96 σ̂/√N",
    xlabel: "samples N", ylabel, xlog: true,
    xlim: [n[0], Math.max(run.target, n[n.length - 1]) * 1.2],
    ylim: [Math.min(...ys) - 0.1 * span, Math.max(...ys) + 0.1 * span],
    series: [
      { label: "running mean Q̄_N", x: n, y: m, color: SLOT[0], width: 2 },
      { label: "95% interval", x: [], y: [], color: "rgba(57, 135, 229, 0.5)", width: 8, inert: true },
      { label: "current estimate", x: [n[0], run.target * 1.2], y: [last, last], color: INK, width: 1, dash: [6, 5], inert: true },
      { label: "Q(E[inputs]), the mean beam", x: [n[0], run.target * 1.2], y: [Qmf * k, Qmf * k], color: SLOT[2], width: 1.5, dash: [3, 3], inert: true },
    ],
    under: (ctx, a) => fill(ctx, a, n, lo, hi, "rgba(57, 135, 229, 0.18)"),
    hover: (_, i) => `N = ${n[i]}\nQ̄ = ${m[i].toPrecision(6)}\n± ${half[i].toPrecision(3)} (95%)`,
  });
}

function drawError(plot: Plot, run: McRun, coarse: boolean): void {
  const acc = run.acc, t = acc.trajectory, scale = Math.abs(acc.q.mean) || 1;
  const n = t.map((c) => c.n), se = t.map((c) => c.sd / Math.sqrt(c.n) / scale);
  const sd = acc.q.sd / scale;
  const bias = coarse && acc.dq.n > 1 ? Math.abs(acc.dq.mean) / scale : NaN;
  const nStar = Number.isFinite(bias) && bias > 0 ? (sd / bias) ** 2 : NaN;
  const nMax = Math.max(run.target, n[n.length - 1], Number.isFinite(nStar) ? 3 * nStar : 0) * 1.5;
  const ext = geometric(Math.max(1, n[0]), nMax);
  const series: Series[] = [
    { label: "sampling error σ̂/√N (measured)", x: n, y: se, color: SLOT[0], width: 2.5 },
    { label: "σ̂/√N extrapolated", x: ext, y: ext.map((N) => sd / Math.sqrt(N)), color: SLOT[0], width: 1.25, dash: [2, 3], inert: true },
  ];
  if (Number.isFinite(bias)) {
    const band = (Z95 * acc.dq.se) / scale;
    series.push(
      { label: "bias indicator |E[Q_ℓ − Q_ℓ₋₁]|", x: [ext[0], nMax], y: [bias, bias], color: SLOT[1], width: 2, inert: true },
      { label: "  ± 95%", x: [ext[0], nMax], y: [bias + band, bias + band], color: SLOT[1], width: 1, alpha: 0.5, dash: [2, 3], inert: true, unlisted: true },
      { label: "  ± 95%", x: [ext[0], nMax], y: [Math.max(bias - band, 1e-17), Math.max(bias - band, 1e-17)], color: SLOT[1], width: 1, alpha: 0.5, dash: [2, 3], inert: true, unlisted: true },
      { label: "root MSE √(bias² + σ²/N)", x: ext, y: ext.map((N) => Math.hypot(bias, sd / Math.sqrt(N))), color: INK, width: 1.5, dash: [6, 4], inert: true },
    );
  }
  const ys = series.flatMap((s) => Array.from(s.y)).filter((v) => v > 0 && Number.isFinite(v));
  plot.draw({
    title: coarse ? "error against samples: MSE = bias² + σ²/N" : "sampling error against samples (level 0 has no parent: no bias indicator)",
    xlabel: "samples N", ylabel: "error relative to |E[Q]|", xlog: true, ylog: true, legend: "tr",
    xlim: [ext[0], nMax], ylim: [Math.min(...ys) / 3, Math.max(...ys) * 3],
    series,
    over: (ctx, a) => {
      if (!Number.isFinite(nStar)) return;
      vline(ctx, a, nStar, SLOT[1], [4, 4]);
      ctx.fillStyle = SLOT[1];
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "left";
      ctx.fillText(` N* = σ²/bias² ≈ ${fmtCount(nStar)}`, a.sx(nStar), a.box.b - 10);
    },
    hover: (_, i) => `N = ${n[i]}\nσ̂/√N = ${se[i].toExponential(2)} (relative)`,
  });
}

function drawBands(plot: Plot, st: State, run: McRun, wMean: Float64Array): void {
  const acc = run.acc, d = displayOf(st), L = d.dimensional ? d.ref.L : 1;
  const W = d.dimensional ? deflectionScale(d.ref, st.load) : 1;
  const mean = acc.field.mean, sd = acc.field.sd();
  const top = Math.max(...Array.from(mean, (m, i) => Math.abs(m) + 2 * sd[i]), ...Array.from(wMean, Math.abs)) * W;
  const u = d.dimensional ? axisUnit(top, "m") : { factor: 1, label: "" };
  const k = W * u.factor, X = Array.from(fieldX, (x) => x * L);
  const at = (z: number) => Array.from(mean, (m, i) => (m + z * sd[i]) * k);
  const series: Series[] = acc.paths.slice(0, 8).map((w, i) => ({
    label: `sample ${i + 1}`, x: X, y: Array.from(w, (v) => v * k), color: SLOT[i % SLOT.length], width: 1, alpha: 0.5, unlisted: true,
  }));
  series.push(
    { label: "mean deflection E[w](x)", x: X, y: at(0), color: SLOT[0], width: 2.5 },
    { label: "± σ, ± 2σ bands", x: [], y: [], color: "rgba(57, 135, 229, 0.45)", width: 8, inert: true },
    { label: "w of the mean beam", x: X, y: Array.from(wMean, (v) => v * k), color: SLOT[2], width: 1.5, dash: [3, 3] },
  );
  const ys = [...at(-2), ...at(2)];
  const hi = Math.max(...ys, 0), lo = Math.min(...ys, 0), pad = 0.08 * (hi - lo || 1);
  plot.draw({
    title: `deflection field over ${acc.n.toLocaleString()} samples — the 'function' quantity of interest`,
    xlabel: d.dimensional ? "x [m]" : "x / L",
    ylabel: d.dimensional ? `w [${u.label}]  (downward)` : `w EI₀ / ${st.load === "uniform" ? "q₀L⁴" : "P₀L³"}  (downward)`,
    // Drawn sagging, so the root's corner at the bottom left is the empty one.
    xlim: [0, L], ylim: [lo - pad, hi + pad], yflip: true, legend: "bl",
    series,
    under: (ctx, a) => {
      fill(ctx, a, X, at(-2), at(2), "rgba(57, 135, 229, 0.13)");
      fill(ctx, a, X, at(-1), at(1), "rgba(57, 135, 229, 0.22)");
    },
    hover: (s, i) => `${s.label}\nx = ${s.x[i].toPrecision(3)}\nw = ${s.y[i].toPrecision(4)}` +
      (s.label.startsWith("mean deflection") ? `\nσ = ${(sd[i] * k).toPrecision(3)}` : ""),
  });
}

function readout(st: State, run: McRun, Qmf: number, head: string, status: string, coarse: boolean, kl: KLData): string {
  const acc = run.acc, d = displayOf(st), q = acc.q;
  const v = (x: number) => valueOf(d, st.qoi, x, st.load);
  const rel = (x: number) => (x / Math.abs(q.mean)).toExponential(2);
  const nw = workers.size, wall = run.wallMs / 1000;
  const lines = [
    head,
    "",
    `${status}: N = ${acc.n.toLocaleString()} of ${run.target.toLocaleString()} on ${nw} worker${nw === 1 ? "" : "s"}` +
      (wall > 0 ? ` — ${fmtCount(acc.n / wall)} samples/s, ${(acc.cpuMs / acc.n).toFixed(2)} ms per sample per worker` : "") +
      (acc.waiting ? ` (${acc.waiting} waiting on an earlier batch)` : ""),
    "",
    `E[Q_ℓ] ≈ ${v(q.mean)} ± ${v(Z95 * q.se)}  (95%: ± 1.96 σ̂/√N, relative ${rel(Z95 * q.se)})`,
    `σ̂ = ${v(q.sd)}  (coefficient of variation ${(100 * q.sd / Math.abs(q.mean)).toFixed(2)}%)`,
    `Q(E[inputs]) = ${v(Qmf)} for the mean beam — the Jensen gap E[Q] − Q(E[·]) is ${rel(q.mean - Qmf)} of E[Q]` +
      (Math.abs(q.mean - Qmf) > Z95 * q.se ? "" : " (not yet resolved by the sampling error)"),
  ];
  if (coarse && acc.dq.n > 1) {
    const dq = acc.dq, bias = Math.abs(dq.mean), band = Z95 * dq.se;
    const nStar = (q.sd / bias) ** 2, resolved = bias > band;
    lines.push(
      "",
      `correction Y = Q_ℓ − Q_ℓ₋₁ (same ω):  E[Y] ≈ ${rel(dq.mean)} ± ${rel(band)},  σ̂_Y = ${rel(dq.sd)}  (relative to E[Q])`,
      `V[Y] / V[Q] = ${(dq.variance / q.variance).toExponential(2)} — the coupled correction is that much less variable than Q itself;` +
        " multilevel Monte Carlo spends its samples on corrections for that reason",
      `MSE split now:  bias² ≲ E[Y]² = ${(bias / Math.abs(q.mean)) ** 2 > 0 ? ((bias / q.mean) ** 2).toExponential(2) : "0"},  ` +
        `σ̂²/N = ${((q.se / q.mean) ** 2).toExponential(2)}  (both relative)`,
      resolved
        ? `N* = σ̂²/E[Y]² ≈ ${fmtCount(nStar)}: ${acc.n < nStar ? "sampling still dominates; more samples help" : "past it, the mesh dominates — more samples on this level buy nothing"}`
        : `E[Y] is not yet resolved from zero (|E[Y]| < its 95% band): the bias is below ~${rel(bias + band)}, and the sampling error still dominates`,
      "|E[Y_ℓ]| bounds the bias of level ℓ when |Q_ℓ − Q| falls geometrically at any rate α ≥ 1 (the bias is then |E[Y]|/(2^α − 1))",
    );
  } else if (!coarse) lines.push("", "level 0 has no parent level: raise ℓ to measure the correction Y = Q_ℓ − Q_ℓ₋₁ and with it the bias.");
  const perEll = st.ell * st.ne0 * (coarse ? 2 ** (st.mcLevel - 1) : 1);
  if (perEll < 1) lines.push(`the ${coarse ? "parent " : ""}mesh has ${plain(perEll, 2)} elements per correlation length: too coarse to see the field it is given`);
  if (termsOf(kl) < st.terms) lines.push(`only ${termsOf(kl)} KL terms are above round-off for this kernel and length`);
  if (st.qoi !== "omega1" && st.load === "point" && st.loadSigma > 0) lines.push(`point load: random magnitude 1 + σ_q ξ′₀ at x = ${qoiPoint(st.supports)}`);
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}

export function runStatus(run: McRun): string {
  return run.error ? `stopped: ${run.error}` : run.running ? "running" : run.done ? "done" : run.paused ? "paused" : "waiting";
}

export function vline(ctx: CanvasRenderingContext2D, a: Axes, x: number, color: string, dash: number[]): void {
  if (!Number.isFinite(x)) return;
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.setLineDash(dash);
  ctx.beginPath();
  ctx.moveTo(a.sx(x), a.box.t);
  ctx.lineTo(a.sx(x), a.box.b);
  ctx.stroke();
  ctx.setLineDash([]);
}

export function geometric(a: number, b: number, n = 60): number[] {
  return Array.from({ length: n }, (_, i) => a * (b / a) ** (i / (n - 1)));
}

export function fmtCount(n: number): string {
  if (!Number.isFinite(n)) return "—";
  return n >= 1e6 ? n.toExponential(1) : Math.round(n).toLocaleString();
}
