// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 4 on screen: the random input, every picture of it at once.
 *
 * - stiffness (full width, on top): a handful of samples of e(x, ω) = EI/EI₀,
 *   the shaded 5–95% band (exact: the field is lognormal pointwise) and the
 *   mean, which is the deterministic section whatever the truncation. Sample 1
 *   is also marked at the quadrature points of a mesh of `ne` elements — all
 *   that a solve on that mesh ever sees of it. With fewer points than wiggles
 *   per correlation length, a coarse level is solving a different beam from a
 *   fine one, and its corrections will not shrink: the lesson of MC_PLAN.md's
 *   last pitfall.
 * - load: the same for q(x, ω), when the load is random.
 * - KL spectrum: λ_j against j, the decay the kernel's smoothness sets, the
 *   exact eigenvalues where they are known, and where the truncation cuts.
 *
 * With the plate as the structure, the field is two-dimensional and separable
 * (`random/field2d.ts`): two samples of log(D/D₀) on the diverging scale, with
 * the mesh of `plateNe` elements drawn over the first; the share of the
 * variance the M kept products hold at each point; and the product spectrum
 * λ_i λ_j, sorted, with the truncation.
 */

import { gaussLegendre } from "../../quad";
import { FieldOnGrid, productTerms } from "../../random/field2d";
import { FieldAt, draw, lognormalQuantile, realise, type Draw } from "../../random/field";
import {
  MAX_TERMS, captured, exponentialEigenvalues, karhunenLoeve, termsOf, truncate, type KLData, type Kernel,
} from "../../random/kl";
import { sectionOf } from "../../hierarchy";
import type { Figure } from "../figure";
import { colourBar, linspace, paint, plateOutline, shapeLimits, type Grid } from "../heatmap";
import { SLOT, type Axes, type Plot, type Series } from "../plot";
import { fieldSpecOf, referenceOf, type State } from "../state";
import { axisUnit, flexuralRigidity, plain, si } from "../units";
import { isPlate, klsOf } from "./structure";
import { ro, type Block, type Tile } from "../readout";
import { memo, type ViewResult } from "./view";

export const KERNEL_NAME: Record<Kernel, string> = {
  exponential: "exponential",
  matern32: "Matérn ν = 3/2",
  gaussian: "squared exponential",
};

const SAMPLES = 8;
const GRID = 257;
const BAND = "rgba(57, 135, 229, 0.16)";
const MEAN_INK = "rgba(207, 238, 255, 0.85)";

const fullKL = memo<KLData>();

/** The expansion for the pane's kernel and length, truncated to its terms; solved once per kernel and length. */
export function klOf(st: State): KLData {
  const full = fullKL(JSON.stringify([st.kernel, st.ell]), () => karhunenLoeve(st.kernel, st.ell, MAX_TERMS));
  return truncate(full, st.terms);
}

const grid = Float64Array.from({ length: GRID }, (_, i) => i / (GRID - 1));
const gridField = memo<FieldAt>();

/** What every panel shares: the expansion on the plot grid, the first samples' normals, the x axis. */
interface Ctx {
  readonly st: State;
  readonly kl: KLData;
  readonly f: FieldAt;
  readonly draws: Draw[];
  readonly X: number[];
  readonly L: number;
  readonly xlabel: string;
}

export function renderField(fig: Figure, st: State): ViewResult {
  if (isPlate(st)) return renderPlateField(fig, st);
  const kl = klOf(st), M = termsOf(kl), spec = fieldSpecOf(st);
  const f = gridField(JSON.stringify([st.kernel, st.ell, M]), () => new FieldAt(kl, grid));
  const d = st.dimensional, L = d ? st.L : 1;
  const ellText = d ? si(st.ell * st.L, "m", 3) : `${plain(st.ell, 3)} L`;
  const header = `${KERNEL_NAME[st.kernel]} kernel, correlation length ℓ = ${ellText}, σ = ${plain(st.sigma, 3)}, ` +
    `M = ${M} term${M === 1 ? "" : "s"}${M < st.terms ? ` (of ${st.terms} asked: the rest are below round-off)` : ""}`;
  const c: Ctx = {
    st, kl, f, L, xlabel: d ? "x [m]" : "x / L",
    draws: Array.from({ length: SAMPLES }, (_, i) => draw(spec, M, st.seed, i)),
    X: Array.from(grid, (x) => x * L),
  };

  // The stiffness across the top; below it the load (when it is random) and the spectrum.
  const randomLoad = st.loadSigma > 0;
  const plots = randomLoad ? fig.panels(3, { cols: 2, wideFirst: true }) : fig.panels(2, { cols: 1 });
  const input: Block[] = [ro.lead(header), ...stiffness(plots[0], c)];
  if (randomLoad) input.push(...load(plots[1], c));
  else input.push(ro.note(`the load is deterministic: raise "load σ_q" to give it a random part, and a panel of its samples joins these`));
  const terms = spectrum(plots[plots.length - 1], c);
  return { readout: { sections: [{ title: "RANDOM INPUT", blocks: input }, { title: "KARHUNEN–LOÈVE SPECTRUM", blocks: terms }] }, animate: false };
}

const seriesOf = (c: Ctx, ys: Float64Array[], scale: number, name: string): Series[] => ys.map((y, i) => ({
  label: `${name}, sample ${i + 1}`, x: c.X, y: Array.from(y, (v) => v * scale), color: SLOT[i % SLOT.length],
  width: i === 0 ? 2 : 1.25, alpha: i === 0 ? 1 : 0.55, unlisted: true,
}));

function stiffness(plot: Plot, c: Ctx): Block[] {
  const { st, kl, f, X, L, xlabel } = c, d = st.dimensional, spec = fieldSpecOf(st);
  const sec = sectionOf(st.section), e0 = Float64Array.from(grid, sec.stiffness ?? (() => 1));
  const ones = new Float64Array(GRID).fill(1);
  const scale = d ? flexuralRigidity(referenceOf(st)) : 1;
  const ys = c.draws.map((dr) => realise(f, spec, dr, e0, ones).e);
  const top = Math.max(...ys.flatMap((y) => Array.from(y))) * scale;
  const u = d ? axisUnit(top, "N·m²") : { factor: 1, label: "" };
  const k = scale * u.factor;
  const lo = lognormalQuantile(f, st.sigma, -1.6449).map((v, i) => v * e0[i] * k);
  const hi = lognormalQuantile(f, st.sigma, 1.6449).map((v, i) => v * e0[i] * k);

  // What a mesh of ne elements samples of sample 1: its Gauss points.
  const rule = gaussLegendre(st.p + 1), pts: number[] = [];
  for (let e = 0; e < st.ne; e++) rule.x.forEach((xi) => pts.push((e + xi) / st.ne));
  const fm = new FieldAt(kl, pts);
  const e0m = Float64Array.from(pts, sec.stiffness ?? (() => 1));
  const seen = realise(fm, spec, c.draws[0], e0m, new Float64Array(pts.length).fill(1)).e;

  plot.draw({
    title: `stiffness e(x, ω) = EI / EI₀ — ${SAMPLES} samples, and what ${st.ne} elements see of sample 1`,
    xlabel,
    ylabel: d ? `EI [${u.label}]` : "EI / EI₀",
    xlim: [0, L],
    series: [
      ...seriesOf(c, ys, k, "e"),
      { label: "mean (the deterministic section)", x: X, y: Array.from(e0, (v) => v * k), color: MEAN_INK, width: 1.25, dash: [6, 5] },
      { label: "5–95% band", x: [], y: [], color: "rgba(57, 135, 229, 0.5)", width: 8, inert: true },
      {
        label: `sample 1 at the ${pts.length} Gauss points of ${st.ne} elements`, x: pts.map((x) => x * L), y: Array.from(seen, (v) => v * k),
        color: SLOT[0], width: 0.001, markers: true,
      },
    ],
    under: (ctx, a) => fill(ctx, a, X, lo, hi, BAND),
    hover: (s, i) => `${s.label}\n${xlabel} = ${s.x[i].toPrecision(3)}\n${s.y[i].toPrecision(4)}`,
  });

  const mid = (GRID - 1) / 2, sM = f.s[mid];
  const perEll = st.ell * st.ne;
  return [
    ro.tiles([
      {
        label: "variance captured", value: `${(100 * captured(kl)).toFixed(2)}%`, detail: "Σ_{j<M} λ_j, of σ²",
        help: "The share of the field's variance the M kept terms carry; the rest is truncated away.",
      },
      {
        label: "coefficient of variation of EI", value: `${(100 * Math.sqrt(Math.expm1(st.sigma ** 2 * sM))).toFixed(1)}%`,
        detail: `√(e^{σ² s_M} − 1), s_M = ${plain(sM, 4)} at midspan`,
        help: "The stiffness's spread relative to its mean at midspan. The mean is exactly the section's at every x, whatever M.",
      },
      perTile(perEll, `mesh of ${st.ne} elements`),
    ]),
    ro.note(st.massFollows ? "mass follows depth: ρA ∝ (EI/EI₀)^{1/3}, a random depth d with I ∝ d³, A ∝ d" : "mass deterministic: the field is a random modulus"),
    ro.note(`seed ${st.seed}: sample i is drawn from (seed, i) alone, so these are samples 1–${SAMPLES} of every Monte Carlo run with this seed`),
  ];
}

function load(plot: Plot, c: Ctx): Block[] {
  const { st, f, X, L, xlabel } = c, spec = fieldSpecOf(st);
  const ones = new Float64Array(GRID).fill(1);
  const ys = c.draws.map((dr) => realise(f, spec, dr, ones, ones).q!);
  const lo = f.s.map((s) => 1 - 1.6449 * st.loadSigma * Math.sqrt(s)), hi = f.s.map((s) => 1 + 1.6449 * st.loadSigma * Math.sqrt(s));
  plot.draw({
    title: `distributed load q(x, ω) = q₀ (1 + σ_q g′) — ${SAMPLES} samples`,
    xlabel, ylabel: "q / q₀", xlim: [0, L],
    series: [...seriesOf(c, ys, 1, "q"), { label: "mean", x: X, y: X.map(() => 1), color: MEAN_INK, width: 1.25, dash: [6, 5], inert: true }],
    under: (ctx, a) => fill(ctx, a, X, lo, hi, BAND),
  });
  const blocks = [
    ro.lead(`load σ_q = ${plain(st.loadSigma, 3)}, independent of the stiffness (its own normals, the same kernel)`),
    ro.note("shaded: the 5–95% band of q, Gaussian pointwise; a point load gets the random magnitude 1 + σ_q ξ′₀ instead"),
  ];
  if (st.load === "point") blocks.push(ro.note("the beam's load is a point load, so the Monte Carlo views use that magnitude, not this field"));
  return blocks;
}

function spectrum(plot: Plot, c: Ctx): Block[] {
  const { st, kl } = c;
  const full = fullKL(JSON.stringify([st.kernel, st.ell]), () => karhunenLoeve(st.kernel, st.ell, MAX_TERMS));
  const shown = Math.min(128, full.spectrum.length);
  const lam = Array.from(full.spectrum.subarray(0, shown)), j = lam.map((_, i) => i + 1);
  const M = termsOf(kl);
  const series: Series[] = [
    { label: "retained λ_j (Nyström)", x: j.slice(0, M), y: lam.slice(0, M), color: SLOT[0], width: 2, markers: true },
    { label: "dropped λ_j", x: j.slice(M - 1), y: lam.slice(M - 1), color: SLOT[0], width: 1.25, alpha: 0.45 },
  ];
  let exact: Float64Array | null = null;
  if (st.kernel === "exponential") {
    exact = exponentialEigenvalues(st.ell, shown);
    series.push({ label: "exact (exponential kernel)", x: j, y: Array.from(exact), color: MEAN_INK, width: 1.25, dash: [6, 5] });
  }
  const pos = lam.filter((v) => v > 0);
  plot.draw({
    title: `Karhunen–Loève eigenvalues — ${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L`,
    xlabel: "mode j", ylabel: "λ_j (share of the variance)", xlog: true, ylog: true,
    xlim: [0.8, shown * 1.2], ylim: [Math.max(Math.min(...pos) / 3, 1e-17), 2],
    series,
    over: (ctx, a) => {
      ctx.strokeStyle = SLOT[1];
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      const X = a.sx(M + 0.5);
      ctx.beginPath(); ctx.moveTo(X, a.box.t); ctx.lineTo(X, a.box.b); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = SLOT[1];
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "left";
      ctx.fillText(` truncation, M = ${M}`, X, a.box.t + 12);
    },
    hover: (s, i) => `${s.label}\nj = ${s.x[i]}\nλ = ${s.y[i].toExponential(3)}`,
  });
  const rows = lam.slice(0, 10).map((v, i) => [
    String(i + 1), v.toExponential(4), exact ? exact[i].toExponential(4) : "—",
    exact ? Math.abs(v / exact[i] - 1).toExponential(1) : "—",
    (100 * lam.slice(0, i + 1).reduce((a, b) => a + b, 0)).toFixed(2) + "%",
  ]);
  const decay = { exponential: "j⁻² (paths continuous, nowhere differentiable)", matern32: "j⁻⁴ (paths once differentiable)", gaussian: "faster than any power (paths analytic)" }[st.kernel];
  return [
    ro.table(["j", "λ_j", "exact", "rel. error", "cumulative"], rows, { caption: `the first ten eigenvalues — decay ${decay}` }),
    ro.note(`retained: ${(100 * captured(kl)).toFixed(3)}% of the variance — the rest is a truncation error shared by every level, ` +
      "so it biases the model, not the hierarchy"),
    ro.note("Nyström on 384 Gauss nodes; against the exponential kernel's exact values the j-th is good to about (j/384)²"),
  ];
}

/** Elements per correlation length, flagged when the mesh is coarser than the field. */
function perTile(perEll: number, mesh: string): Tile {
  return {
    label: "elements per correlation length", value: plain(perEll, 3), detail: mesh,
    verdict: perEll < 1 ? { text: "coarser than the field", tone: "warn" } : { text: "resolves the field", tone: "good" },
    help: "How many elements of the mesh a correlation length spans. Below one, the mesh is coarser than the field: a solve on it cannot see the wiggles it averages over.",
  };
}

/** Shade between two curves. */
export function fill(
  ctx: CanvasRenderingContext2D, a: { sx: (x: number) => number; sy: (y: number) => number },
  x: ArrayLike<number>, lo: ArrayLike<number>, hi: ArrayLike<number>, color: string,
): void {
  ctx.fillStyle = color;
  ctx.beginPath();
  for (let i = 0; i < x.length; i++) (i ? ctx.lineTo : ctx.moveTo).call(ctx, a.sx(x[i]), a.sy(hi[i]));
  for (let i = x.length - 1; i >= 0; i--) ctx.lineTo(a.sx(x[i]), a.sy(lo[i]));
  ctx.closePath();
  ctx.fill();
}

const PLATE_RES = 97;
const plateGrids = memo<{ f: FieldOnGrid; samples: Grid[] }>();

/** The plate's field: two samples, the variance the truncation keeps, and the product spectrum. */
function renderPlateField(fig: Figure, st: State): ViewResult {
  const { kl, klY, M } = klsOf(st), spec = fieldSpecOf(st), a = st.aspect, kly = klY ?? kl;
  const nx = PLATE_RES, ny = Math.max(9, Math.round(PLATE_RES / a) | 1);
  const { f, samples } = plateGrids(JSON.stringify([st.kernel, st.ell, a, M, st.sigma, st.seed]), () => {
    const f = new FieldOnGrid(kl, kly, st.terms, linspace(0, a, nx), linspace(0, 1, ny), a);
    const ones = new Float64Array(f.n).fill(1);
    const samples = [0, 1].map((i): Grid => {
      const e = realise(f, spec, draw(spec, f.M, st.seed, i), ones, ones).e;
      return { nx, ny, x0: 0, x1: a, y0: 0, y1: 1, v: e.map(Math.log) };
    });
    return { f, samples };
  });
  const [s1, s2, share, spectrumPlot] = fig.panels(4, { cols: 2 });
  const top = 3 * st.sigma || 1;
  const d = st.dimensional, L = d ? st.L : 1;
  const axis = d ? ["x [m]", "y [m]"] : ["x / L", "y / L"];
  const inPlate = (ax: Axes): Axes => ({ ...ax, sx: (x: number) => ax.sx(x * L), sy: (y: number) => ax.sy(y * L) });
  samples.forEach((g, i) => {
    const plot = i === 0 ? s1 : s2;
    plot.draw({
      title: `log(D / D₀), sample ${i + 1}${i === 0 ? ` — and the ${st.plateNe} × ${st.plateNe} mesh` : ""}`,
      xlabel: axis[0], ylabel: axis[1], ...shapeLimits(plot, a * L, L), series: [],
      under: (ctx, ax) => {
        const pa = inPlate(ax);
        paint(ctx, pa, g, "diverging", top);
        if (i > 0) return;
        ctx.strokeStyle = "rgba(5, 5, 12, 0.6)";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let e = 1; e < st.plateNe; e++) {
          const X = pa.sx((a * e) / st.plateNe), Y = pa.sy(e / st.plateNe);
          ctx.moveTo(X, pa.sy(0)); ctx.lineTo(X, pa.sy(1));
          ctx.moveTo(pa.sx(0), Y); ctx.lineTo(pa.sx(a), Y);
        }
        ctx.stroke();
      },
      over: (ctx, ax) => {
        plateOutline(ctx, inPlate(ax), a, 1, "FFFF");
        colourBar(ctx, ax, "diverging", top, "log D/D₀", (v) => plain(v, 2));
      },
    });
  });
  share.draw({
    title: `pointwise variance kept by ${f.M} products, s_M(x, y)`,
    xlabel: axis[0], ylabel: axis[1], ...shapeLimits(share, a * L, L), series: [],
    under: (ctx, ax) => paint(ctx, inPlate(ax), { nx, ny, x0: 0, x1: a, y0: 0, y1: 1, v: f.s }, "sequential", 1),
    over: (ctx, ax) => {
      plateOutline(ctx, inPlate(ax), a, 1, "FFFF");
      colourBar(ctx, ax, "sequential", 1, "s_M", (v) => plain(v, 2));
    },
  });

  // Every product the two 1D expansions hold, sorted, and where the truncation cuts.
  const all = productTerms(kl, kly, Infinity).values, shown = Math.min(all.length, 1024);
  const j = Array.from({ length: shown }, (_, i) => i + 1), lam = Array.from(all.subarray(0, shown));
  spectrumPlot.draw({
    title: `product spectrum λ_i^x λ_j^y — ${KERNEL_NAME[st.kernel]}, ℓ = ${plain(st.ell, 3)} L`,
    xlabel: "term t (products, sorted)", ylabel: "λ_t (share of the variance)", xlog: true, ylog: true,
    xlim: [0.8, shown * 1.2], ylim: [Math.max(Math.min(...lam.filter((v) => v > 0)) / 3, 1e-17), 2],
    series: [
      { label: "kept", x: j.slice(0, f.M), y: lam.slice(0, f.M), color: SLOT[0], width: 2 },
      { label: "dropped", x: j.slice(f.M - 1), y: lam.slice(f.M - 1), color: SLOT[0], width: 1.25, alpha: 0.45 },
    ],
    over: (ctx, ax) => {
      ctx.strokeStyle = SLOT[1];
      ctx.lineWidth = 1.5;
      ctx.setLineDash([4, 4]);
      const X = ax.sx(f.M + 0.5);
      ctx.beginPath(); ctx.moveTo(X, ax.box.t); ctx.lineTo(X, ax.box.b); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = SLOT[1];
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "left";
      ctx.fillText(` truncation, M = ${f.M}`, X, ax.box.t + 12);
    },
    hover: (se, i) => `${se.label}\nt = ${se.x[i]}\nλ = ${se.y[i].toExponential(3)}`,
  });

  const kept = f.terms.values.reduce((x, y) => x + y, 0), total = all.reduce((x, y) => x + y, 0);
  const sM = f.s[((ny - 1) / 2) * nx + (nx - 1) / 2];
  const perEll = st.ell * st.plateNe;
  const blocks: Block[] = [
    ro.lead(`${KERNEL_NAME[st.kernel]} kernel, separable: C₁(|Δx|) C₁(|Δy|) with ℓ = ${d ? si(st.ell * st.L, "m", 3) : `${plain(st.ell, 3)} L`} along both sides, ` +
      `σ = ${plain(st.sigma, 3)}, M = ${f.M} products${f.M < st.terms ? ` (of ${st.terms} asked)` : ""}`),
    ro.tiles([
      {
        label: "variance kept", value: `${(100 * kept).toFixed(2)}%`, detail: `the two 1D expansions hold ${(100 * total).toFixed(2)}%`,
        help: "A separable field needs more terms than a 1D one for the same share: two decaying spectra multiplied and sorted decay more slowly than either.",
      },
      {
        label: "coefficient of variation of D", value: `${(100 * Math.sqrt(Math.expm1(st.sigma ** 2 * sM))).toFixed(1)}%`,
        detail: `√(e^{σ² s_M} − 1), s_M = ${plain(sM, 4)} at the centre`,
      },
      perTile(perEll, `plate mesh of ${st.plateNe} × ${st.plateNe}`),
    ]),
    ro.note(st.kernel === "gaussian"
      ? "the squared exponential is the one kernel whose product is also isotropic, e^{−r²/2ℓ²}."
      : "for this kernel the product is not isotropic: correlation falls faster along a diagonal than along an axis."),
    ro.note(st.massFollows ? "mass follows depth: ρt ∝ (D/D₀)^{1/3}, a random thickness t with D ∝ t³" : "mass deterministic: the field is a random modulus"),
    ro.note(`seed ${st.seed}: these are samples 1 and 2 of every Monte Carlo run on this plate with this seed`),
  ];
  if (st.loadSigma > 0) blocks.push(ro.note(`load σ_q = ${plain(st.loadSigma, 3)}: the pressure is 1 + σ_q g′(x, y), g′ an independent field from the same expansion`));
  return { readout: { sections: [{ title: "RANDOM INPUT ON THE PLATE", blocks }] }, animate: false };
}
