// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 8 on screen: multilevel Monte Carlo at work, the solver beside its
 * statistics.
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
 * Below, the statistics of the same run: the distribution of Q₀ with the
 * shown sample and the MLMC estimate marked, the variance of Q_ℓ and Y_ℓ
 * against level, and ε² × cost against ε. It is the MLMC view's run — the
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
import { SWEEP, drawCost, drawVariance, fmt2, mlmcSession, type MlmcSession } from "./mlmc";
import { fmtCount, runStatus, vline } from "./montecarlo";
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
  const [solver, hist, variance, cost] = fig.panels(4, { cols: 3, wideFirst: true });
  const lv = Math.max(1, Math.min(st.liveLevel, c.levels - 1));
  const acc = c.run.streams[lv].acc;
  if (acc.n < 1) {
    solver.draw({ xlabel: "", ylabel: "", series: [], title: `waiting for the first samples of level ${lv}…` });
    for (const pl of [hist, variance, cost]) pl.draw({ xlabel: "", ylabel: "", series: [] });
    return { readout: summary(c, lv, null), animate: true };
  }

  const i = Math.floor(t / DWELL_MS) % Math.min(acc.n, CYCLE);
  const sampler = samplers(JSON.stringify([c.spec, c.M]), () => new Sampler(c.spec, c.kl, c.klY));
  const pair = pairs(JSON.stringify([c.spec, c.M, c.neOf(0), lv, i]), () => solvePair(c, sampler, lv, i));
  const stored = { Q: acc.values[i], Qc: acc.coarseValues[i] };
  const th = (2 * Math.PI * t) / PERIOD_MS;

  if (isPlate(st)) drawPlates(solver, c, lv, pair, th);
  else drawBeams(solver, c, lv, pair, th);
  drawHistogram(hist, c, pair);
  if (c.S.length) {
    drawVariance(variance, c);
    drawCost(cost, c);
  } else for (const pl of [variance, cost]) pl.draw({ xlabel: "", ylabel: "", series: [], title: "waiting for the survey…" });
  return { readout: summary(c, lv, { pair, stored }), animate: true };
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
  // Below the curves, in screen space: the field as a strip, then one row of mesh ticks per level.
  const extra = 0.62 * (hi - lo);
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

function drawHistogram(plot: Plot, c: MlmcSession, p: Pair): void {
  const acc = c.run.streams[0].acc, { mean, sd } = acc.q;
  if (acc.n < 2) { plot.draw({ xlabel: "Q", ylabel: "", series: [] }); return; }
  const h = histogram(acc.values, acc.min, acc.max, sd), edges = h.edges, dens = h.density;
  const best = (c.sweep?.runs ?? []).filter((r) => r.status === "converged").pop();
  const est = best ? best.estimate : NaN, half = best ? Z95 * Math.sqrt(best.sampleVariance) : NaN;
  const pad = 0.04 * (edges[edges.length - 1] - edges[0] || 1);
  const lo = Math.min(edges[0], p.fine.Q, p.coarse.Q) - pad, hi = Math.max(edges[edges.length - 1], p.fine.Q, p.coarse.Q) + pad;
  const topD = Math.max(...dens) * 1.2;
  const series: Series[] = [
    { label: `Q₀, ${acc.n.toLocaleString()} samples`, x: [], y: [], color: "rgba(57, 135, 229, 0.55)", width: 8, inert: true },
    { label: "bin", x: Array.from(dens, (_, k) => 0.5 * (edges[k] + edges[k + 1])), y: Array.from(dens), color: "rgba(0, 0, 0, 0)", width: 0.001, unlisted: true },
    { label: "E[Q₀]", x: [], y: [], color: MUTED, width: 2, inert: true },
    { label: "the shown sample, Q_ℓ", x: [], y: [], color: SLOT[0], width: 2, inert: true },
  ];
  if (best) series.push({ label: "MLMC estimate ± 95%", x: [], y: [], color: INK, width: 2, inert: true });
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

function summary(c: MlmcSession, lv: number, shown: { pair: Pair; stored: { Q: number; Qc: number } } | null): string {
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
  if (c.S.length > 1) {
    const S = c.S, alpha = -slope(S.map((s) => Math.abs(s.meanY))), beta = -slope(S.map((s) => s.varY)), gamma = slope(S.map((_, l) => c.cost(l)));
    lines.push("", `survey of ${S.length} level${S.length > 1 ? "s" : ""}: α ≈ ${fmt2(alpha)}, β ≈ ${fmt2(beta)}, γ = ${fmt2(gamma)}; ` +
      `V[Y_ℓ]/V[Q_ℓ] = ${S.slice(1).map((s) => (s.varY / s.varQ).toExponential(1)).join(", ")}`);
  }
  const runs = c.sweep?.runs ?? [];
  if (runs.length) {
    lines.push("", `tolerances (θ = ${THETA}): ` + runs.map((r, k) => `ε = ${(SWEEP[k] * st.mlEps).toExponential(1)} ${r.status}` +
      (r.status === "converged" ? ` (L = ${r.L}, ${v(r.estimate)} ± ${(Z95 * Math.sqrt(r.sampleVariance) / c.scale).toExponential(1)})` : "")).join("; "));
  }
  if (d.dimensional) lines.push(referenceNote(d.ref));
  return lines.join("\n");
}
