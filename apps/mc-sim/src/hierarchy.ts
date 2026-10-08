// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The hierarchy of discretisations, and how fast a quantity of interest settles
 * along it.
 *
 * Level ℓ is the same spline space (degree p, continuity k) on ne₀·2^ℓ
 * elements. Multilevel Monte Carlo rests on three rates measured on exactly this
 * ladder:
 *
 *   |Q_ℓ − Q|        ~ h_ℓ^α     the bias a level leaves,
 *   |Q_ℓ − Q_{ℓ−1}|  ~ h_ℓ^α     what each correction adds — computable without Q,
 *   cost_ℓ           ~ h_ℓ^{−γ}  what a sample costs.
 *
 * (The variance rate β of the corrections joins them once the input is random.)
 * The successive differences are the point: they are what MLMC actually
 * estimates, and they track the true error wherever the true error is known —
 * which is why the view plots both.
 *
 * Quantities of interest, all nondimensional here:
 *
 *   ω₁          first natural frequency, √λ₁,
 *   deflection  w at the tip of a cantilever, at midspan otherwise,
 *   compliance  ℓ(w) = a(w, w), the work done by the load,
 *   field       ‖w‖ = (∫₀¹ w²)^½, an RMS deflection — a function's error, not a
 *               functional's: its difference is ‖w_ℓ − w_{ℓ−1}‖, not a
 *               difference of norms;
 *   response    |w(x_q)| under the load applied harmonically at Ω, with Rayleigh
 *               damping: the steady amplitude of the forced vibration.
 *
 * The forced response solves
 *
 *   (K − Ω²M + iΩC) u = F,   C = a M + b K,
 *
 * and reads |u| at the same point as the deflection. Ω is set as a fraction of
 * ω₁ of the uniform beam on the same supports, and a, b so that the first two
 * modes of that beam have the damping ratio ζ (the usual way Rayleigh damping is
 * fitted: ζ_n = a/(2ω_n) + bω_n/2). Damping is not optional: undamped, the
 * amplitude is unbounded wherever a sample's ω₁ meets Ω, and so is its variance.
 */

import { complexSolve } from "./band";
import { Beam, SUPPORTS, type Load, type SupportName } from "./beam/beam";
import { exactEigenvalues, exactStatic, pinnedResponse, type StaticSolution } from "./beam/exact";
import { gaussLegendre } from "./quad";
import { evalField } from "./spline";

export type QoI = "omega1" | "deflection" | "compliance" | "field" | "response";
export type LoadCase = "uniform" | "point";
export type Section = "uniform" | "tapered";

export interface BeamCase {
  readonly supports: SupportName;
  readonly load: LoadCase;
  readonly section: Section;
  /** How the load is applied harmonically, for the forced response; ignored by every other quantity. */
  readonly forcing?: Forcing;
}

export interface Forcing {
  /** Ω / ω₁ of the uniform beam on the same supports. */
  readonly ratio: number;
  /** Rayleigh damping ratio of that beam's first two modes. */
  readonly zeta: number;
}

export const DEFAULT_FORCING: Forcing = { ratio: 0.8, zeta: 0.02 };

/** The forcing frequency and the Rayleigh coefficients, C = a M + b K. */
export interface Harmonic {
  readonly Omega: number;
  readonly a: number;
  readonly b: number;
}

export function harmonicOf(c: BeamCase): Harmonic {
  const { ratio, zeta } = c.forcing ?? DEFAULT_FORCING;
  const [l1, l2] = exactEigenvalues(SUPPORTS[c.supports], 2)!;
  const w1 = Math.sqrt(l1), w2 = Math.sqrt(l2);
  return { Omega: ratio * w1, a: (2 * zeta * w1 * w2) / (w1 + w2), b: (2 * zeta) / (w1 + w2) };
}

/** Where deflection is read, and where the point load acts: the tip, or midspan. */
export const qoiPoint = (s: SupportName): number => (s === "cantilever" ? 1 : 0.5);

/**
 * A tapered section: depth falling linearly to half at x = 1, width constant,
 * so I ∝ depth³ and A ∝ depth. No closed form — the hierarchy then measures
 * itself against itself, as it must once the input is random.
 */
const taper = (x: number) => 1 - 0.5 * x;
export function sectionOf(section: Section): { stiffness?: (x: number) => number; mass?: (x: number) => number } {
  return section === "tapered" ? { stiffness: (x) => taper(x) ** 3, mass: taper } : {};
}

export function loadOf(c: BeamCase): Load {
  return c.load === "uniform" ? { q: () => 1 } : { forces: [{ x: qoiPoint(c.supports), P: 1 }] };
}

export function exactStaticOf(c: BeamCase): StaticSolution | null {
  if (c.section !== "uniform") return null;
  return exactStatic(SUPPORTS[c.supports], c.load === "uniform"
    ? { q: 1, forces: [] }
    : { q: 0, forces: [{ x: qoiPoint(c.supports), P: 1 }] });
}

export interface Level {
  readonly level: number;
  readonly ne: number;
  readonly h: number;
  readonly dofs: number;
  readonly Q: number;
  /** |Q_ℓ − Q_{ℓ−1}| (‖w_ℓ − w_{ℓ−1}‖ for the field); NaN on the coarsest level. */
  readonly dQ: number;
  /** |Q_ℓ − Q| against the closed form; NaN where there is none. */
  readonly err: number;
  /**
   * Round-off in Q: the change when K, M and F are perturbed entrywise by
   * (p + 1)ε. It grows as h falls (cond K ~ h⁻⁴); below it, refining buys nothing.
   */
  readonly noise: number;
  /** Wall time to assemble and solve, ms. */
  readonly ms: number;
}

export interface Hierarchy {
  readonly levels: Level[];
  readonly exact: number | null;
  /** Rate fitted to the successive differences, and to the true error where known. */
  readonly alpha: number | null;
  readonly alphaExact: number | null;
  /** Rate at which dofs grow as h falls — 1 for any 1D space. */
  readonly gamma: number | null;
}

export interface HierarchySpec {
  readonly p: number;
  readonly k: number;
  readonly beam: BeamCase;
  readonly qoi: QoI;
  readonly ne0: number;
  /** Number of levels, ℓ = 0 … levels − 1. */
  readonly levels: number;
}

/** What theory predicts for α, given smooth data; null where it says nothing simple. */
export function theoryRate(qoi: QoI, p: number, beam: BeamCase): number | null {
  if (qoi === "omega1") return 2 * (p - 1);
  if (beam.load === "point") return null; // w‴ jumps under the load: the smooth-data rates do not apply
  if (qoi === "compliance") return 2 * (p - 1);
  if (qoi === "field") return Math.min(p + 1, 2 * (p - 1));
  return null; // a point value: no single rate covers it
}

/**
 * Q for one solved beam — the one definition the hierarchy and every Monte
 * Carlo sample share. `c` is the static solution (empty for ω₁, which needs
 * none); for the forced response it is the real part of the amplitude and `ci`
 * the imaginary.
 */
export function evaluateQoI(
  beam: Beam, qoi: QoI, load: Load, xq: number, h?: Harmonic,
): { Q: number; c: Float64Array; ci?: Float64Array } {
  if (qoi === "omega1") return { Q: Math.sqrt(beam.modes(1).values[0]), c: new Float64Array(0) };
  const F = beam.load(load);
  if (qoi === "response") {
    if (!h) throw new Error("the forced response needs a forcing frequency and damping");
    const { Omega: W, a, b } = h, n = beam.space.n;
    // (1 + iΩb) K + (−Ω² + iΩa) M.
    const u = complexSolve(beam.K, beam.M, [1, W * b], [-W * W, W * a], F.subarray(beam.lo, beam.hi));
    const c = new Float64Array(n), ci = new Float64Array(n);
    c.set(u.re, beam.lo);
    ci.set(u.im, beam.lo);
    return { Q: Math.hypot(beam.evaluate(c, xq)[0], beam.evaluate(ci, xq)[0]), c, ci };
  }
  const c = beam.solve(F);
  const Q = qoi === "deflection" ? beam.evaluate(c, xq)[0]
    : qoi === "compliance" ? F.reduce((s, f, i) => s + f * c[i], 0)
    : fieldNorm(beam, c, null);
  return { Q, c };
}

const now = () => globalThis.performance?.now() ?? Date.now();

export function runHierarchy(spec: HierarchySpec): Hierarchy {
  const { p, k, beam: bc, qoi, ne0 } = spec;
  const section = sectionOf(bc.section), load = loadOf(bc);
  const xq = qoiPoint(bc.supports);
  const exactS = exactStaticOf(bc);
  const harmonic = qoi === "response" ? harmonicOf(bc) : undefined;
  let exact: number | null = null;
  if (qoi === "omega1") {
    const ev = bc.section === "uniform" ? exactEigenvalues(SUPPORTS[bc.supports], 1) : null;
    exact = ev ? Math.sqrt(ev[0]) : null;
  } else if (qoi === "response") {
    exact = bc.section === "uniform" && bc.supports === "pinned–pinned" ? pinnedResponse(bc.load, harmonic!) : null;
  } else if (exactS) {
    exact = qoi === "deflection" ? exactS.w(xq) : qoi === "compliance" ? exactS.compliance : rms((x) => exactS.w(x));
  }

  const levels: Level[] = [];
  let prev: { beam: Beam; c: Float64Array; Q: number; noise: number } | null = null;
  const measure = (beam: Beam) => evaluateQoI(beam, qoi, load, xq, harmonic);
  for (let l = 0; l < spec.levels; l++) {
    const ne = ne0 * 2 ** l;
    const base = { p, k, ne, supports: SUPPORTS[bc.supports], ...section };
    const t0 = now();
    const beam = new Beam(base);
    const { Q, c } = measure(beam);
    const ms = now() - t0;
    // Round-off, measured rather than assumed: the same level re-solved with
    // every assembled entry perturbed at the size f64 perturbs it anyway.
    let noise = 4 * EPS * Math.abs(Q);
    for (const seed of [1, 2]) {
      const pb = new Beam({ ...base, perturb: { size: (p + 1) * EPS, seed } });
      const r = measure(pb);
      noise = Math.max(noise, qoi === "field" ? fieldNorm(beam, c, { beam: pb, c: r.c }) : Math.abs(r.Q - Q));
    }
    let dQ = NaN, err = NaN;
    if (prev) dQ = qoi === "field" ? fieldNorm(beam, c, prev) : Math.abs(Q - prev.Q);
    if (exact !== null) err = qoi === "field" ? fieldError(beam, c, exactS!) : Math.abs(Q - exact);
    levels.push({ level: l, ne, h: 1 / ne, dofs: beam.dofs, Q, dQ, err, noise, ms });
    prev = { beam, c, Q, noise };
  }

  const h = levels.map((v) => v.h), noise = levels.map((v) => v.noise);
  // A difference carries the round-off of both its levels.
  const dNoise = noise.map((n, i) => Math.max(n, noise[i - 1] ?? n));
  return {
    levels, exact,
    alpha: fitRate(h, levels.map((v) => v.dQ), dNoise),
    alphaExact: exact === null ? null : fitRate(h, levels.map((v) => v.err), noise),
    gamma: fitRate(h, levels.map((v) => v.dofs), h.map(() => 0), -1),
  };
}

const EPS = 2 ** -52;

/**
 * Least-squares slope of log e against log h over the finest three points that
 * stand clear of round-off — at least `CLEAR` times the measured floor — so the
 * asymptotic regime is fitted, not the pre-asymptotic coarse levels nor the
 * fine ones where round-off has taken over. `sign` = −1 fits a growth rate.
 *
 * The measured floor is a typical perturbation's effect, not the worst one's,
 * and near a resonance it can sit an order of magnitude under the round-off a
 * solve actually makes. So an error is also taken to have reached round-off
 * where it stops falling: nothing past its smallest value is fitted.
 */
export const CLEAR = 10;
export function fitRate(h: number[], e: number[], floor: number[], sign = 1): number | null {
  let last = e.length - 1;
  if (sign > 0) {
    let best = Infinity;
    e.forEach((v, i) => { if (Number.isFinite(v) && v > 0 && v < best) { best = v; last = i; } });
  }
  const pts = h.map((hi, i) => [Math.log(hi), Math.log(e[i])] as const)
    .filter(([, le], i) => i <= last && Number.isFinite(le) && e[i] > CLEAR * floor[i])
    .slice(-3);
  if (pts.length < 2) return null;
  const mx = pts.reduce((s, [x]) => s + x, 0) / pts.length, my = pts.reduce((s, [, y]) => s + y, 0) / pts.length;
  let sxy = 0, sxx = 0;
  for (const [x, y] of pts) { sxy += (x - mx) * (y - my); sxx += (x - mx) ** 2; }
  return (sign * sxy) / sxx;
}

/** (∫₀¹ f²)^½ by Gauss on a fixed fine partition. */
function rms(f: (x: number) => number, parts = 256): number {
  const { x, w } = gaussLegendre(8);
  let s = 0;
  for (let e = 0; e < parts; e++) x.forEach((xi, q) => (s += (w[q] / parts) * f((e + xi) / parts) ** 2));
  return Math.sqrt(s);
}

/** ‖w_ℓ‖, or ‖w_ℓ − w_{ℓ−1}‖ — the coarse mesh nests in the fine, so fine Gauss is exact. */
function fieldNorm(beam: Beam, c: Float64Array, coarse: { beam: Beam; c: Float64Array } | null): number {
  return quadFine(beam, (x) => beam.evaluate(c, x)[0] - (coarse ? evalField(coarse.beam.space, coarse.c, x)[0] : 0));
}

function fieldError(beam: Beam, c: Float64Array, ex: StaticSolution): number {
  return quadFine(beam, (x) => beam.evaluate(c, x)[0] - ex.w(x));
}

function quadFine(beam: Beam, f: (x: number) => number): number {
  const { x, w } = gaussLegendre(beam.spec.p + 3), { breaks, ne } = beam.space;
  let s = 0;
  for (let e = 0; e < ne; e++) {
    const a = breaks[e], h = breaks[e + 1] - a;
    x.forEach((xi, q) => (s += h * w[q] * f(a + h * xi) ** 2));
  }
  return Math.sqrt(s);
}
