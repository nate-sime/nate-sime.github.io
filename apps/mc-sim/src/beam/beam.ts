// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The Euler–Bernoulli beam, nondimensional, in a spline space.
 *
 *   (e(x) w″)″ + μ(x) ẅ = q(x) + Σ P_j δ(x − a_j)       on (0, 1),
 *
 * lengths in units of the span L, stiffness e = EI/EI₀ and mass μ = ρA/(ρA)₀
 * relative to their reference values, so a uniform beam has e = μ = 1 and its
 * natural frequencies are ω̂_n = (β_n L)² — `ui/units.ts` puts the metres and
 * hertz back. The weak form is
 *
 *   a(w, v) = ∫ e w″ v″,   m(w, v) = ∫ μ w v,
 *   ℓ(v)    = ∫ q v + Σ P_j v(a_j) + Σ M_j v′(a_j),
 *
 * which asks for w ∈ H², hence a C¹ space: degree p ≥ 2 and continuity k ≥ 1.
 * Natural conditions (zero moment, zero shear at a free end; zero moment at a
 * pin) are carried by the weak form; essential ones are imposed strongly. With
 * an open knot vector w(0) is c₀ alone and w′(0) involves only c₀, c₁, so
 *
 *   clamped  ⇒ c₀ = c₁ = 0,   pinned ⇒ c₀ = 0,   free ⇒ nothing
 *
 * (mirrored at x = 1). The constrained functions are therefore always a prefix
 * and a suffix, and the free system is a contiguous principal block [lo, hi) of
 * the full banded matrices — still banded, no renumbering.
 *
 * Coefficients arrive sampled at the quadrature points (`sampleAt`), never as
 * callbacks inside the element loop: that is the form a random-field sample
 * takes, and the form a GPU assembly kernel would read. A Monte Carlo run
 * builds one tabulation per level and hands it to every sample's beam, with the
 * sample's stiffness, mass and load already evaluated at its points.
 */

import { SymBand, cholSolve, cholesky } from "../band";
import { bandSpectrum, lowestModes, type Modes } from "../eig";
import { gaussLegendre } from "../quad";
import {
  basisRow, evalField, firstDof, sampleAt, splineSpace, tabulate,
  type SplineSpace, type Tabulation,
} from "../spline";

/** A coefficient: a function of x, or its values already at the tabulation's quadrature points. */
export type Coefficient = ((x: number) => number) | Float64Array;

const sampled = (t: Tabulation, c: Coefficient): Float64Array => (typeof c === "function" ? sampleAt(t, c) : c);

export type End = "clamped" | "pinned" | "free";

export interface Supports {
  readonly left: End;
  readonly right: End;
}

export const SUPPORTS = {
  cantilever: { left: "clamped", right: "free" },
  "clamped–clamped": { left: "clamped", right: "clamped" },
  "pinned–pinned": { left: "pinned", right: "pinned" },
  "clamped–pinned": { left: "clamped", right: "pinned" },
} as const satisfies Record<string, Supports>;

export type SupportName = keyof typeof SUPPORTS;

export interface PointForce { readonly x: number; readonly P: number; }
export interface PointMoment { readonly x: number; readonly M: number; }

export interface Load {
  readonly q?: Coefficient;
  readonly forces?: readonly PointForce[];
  readonly moments?: readonly PointMoment[];
}

export interface BeamSpec {
  readonly p: number;
  readonly ne: number;
  /** Continuity; defaults to the maximal C^{p−1}. */
  readonly k?: number;
  readonly supports: Supports;
  readonly stiffness?: Coefficient;
  readonly mass?: Coefficient;
  /** Gauss points per element; defaults to p + 1, exact for a uniform mass matrix. */
  readonly nq?: number;
  /**
   * Round-off probe: multiply every assembled entry of K, M and F by (1 + δ),
   * δ uniform in ±size, seeded. A backward-error model of what f64 does to the
   * solve anyway, used to measure how much of a result is round-off.
   */
  readonly perturb?: { readonly size: number; readonly seed: number };
}

const CONSTRAINED: Record<End, number> = { clamped: 2, pinned: 1, free: 0 };

/** Why a spec cannot be solved, or null if it can. */
export function admissible(spec: BeamSpec): string | null {
  const k = spec.k ?? spec.p - 1;
  if (spec.p < 2 || k < 1)
    return "The beam's weak form needs w″ to be square-integrable, so the space must be at least C¹ (p ≥ 2, k ≥ 1).";
  const { left, right } = spec.supports;
  // Two constraints in all is exactly the line between a mechanism and a beam:
  // pinned–free and free–free are short of it, every other pairing reaches it.
  if (CONSTRAINED[left] + CONSTRAINED[right] < 2)
    return "These supports leave a rigid-body motion free: the beam is a mechanism.";
  if (freeDofs(spec) < 1) return "The mesh is too coarse to leave any free coefficient.";
  return null;
}

/** Free coefficients of a spec, without building it: the size of every solve, and the unit Monte Carlo costs are counted in. */
export function freeDofs(spec: BeamSpec): number {
  const n = splineSpace(spec.p, spec.ne, spec.k ?? spec.p - 1).n;
  return n - CONSTRAINED[spec.supports.left] - CONSTRAINED[spec.supports.right];
}

export class Beam {
  readonly space: SplineSpace;
  readonly tab: Tabulation;
  /** Free coefficients are [lo, hi). */
  readonly lo: number;
  readonly hi: number;
  /** Full matrices, every coefficient. */
  readonly Kfull: SymBand;
  readonly Mfull: SymBand;
  /** Free blocks. */
  readonly K: SymBand;
  readonly M: SymBand;
  private LK?: SymBand;
  private readonly jitter?: () => number;

  /**
   * `tab`, if given, is reused rather than rebuilt — it must tabulate this
   * spec's space to second derivatives, and it then fixes the quadrature.
   */
  constructor(readonly spec: BeamSpec, tab?: Tabulation) {
    const why = admissible(spec);
    if (why) throw new Error(why);
    const k = spec.k ?? spec.p - 1;
    if (tab) {
      const s = tab.space;
      if (s.p !== spec.p || s.k !== k || s.ne !== spec.ne || tab.d < 2)
        throw new Error("the tabulation given is not of this beam's space");
    }
    this.space = tab?.space ?? splineSpace(spec.p, spec.ne, k);
    this.tab = tab ?? tabulate(this.space, gaussLegendre(spec.nq ?? spec.p + 1), 2);
    const e = sampled(this.tab, spec.stiffness ?? (() => 1));
    const mu = sampled(this.tab, spec.mass ?? (() => 1));
    this.Kfull = assemble(this.tab, e, 2);
    this.Mfull = assemble(this.tab, mu, 0);
    if (spec.perturb) {
      this.jitter = jitter(spec.perturb.size, spec.perturb.seed);
      this.Kfull.data.forEach((v, i, a) => (a[i] = v * this.jitter!()));
      this.Mfull.data.forEach((v, i, a) => (a[i] = v * this.jitter!()));
    }
    this.lo = CONSTRAINED[spec.supports.left];
    this.hi = this.space.n - CONSTRAINED[spec.supports.right];
    this.K = this.Kfull.slice(this.lo, this.hi);
    this.M = this.Mfull.slice(this.lo, this.hi);
  }

  /** Number of free coefficients — the size of every solve. */
  get dofs(): number { return this.hi - this.lo; }

  /** Load vector over every coefficient. */
  load(load: Load): Float64Array {
    const { tab, space } = this, P1 = space.p + 1, stride = (tab.d + 1) * P1;
    const F = new Float64Array(space.n);
    if (load.q) {
      const q = sampled(tab, load.q);
      for (let e = 0; e < space.ne; e++) {
        const f0 = firstDof(space, e);
        for (let iq = 0; iq < tab.nq; iq++) {
          const i = e * tab.nq + iq, wq = tab.w[i] * q[i];
          for (let a = 0; a < P1; a++) F[f0 + a] += wq * tab.B[i * stride + a];
        }
      }
    }
    for (const { x, P } of load.forces ?? []) basisRow(space, x, 0).forEach((v, i) => (F[i] += P * v));
    for (const { x, M } of load.moments ?? []) basisRow(space, x, 1).forEach((v, i) => (F[i] += M * v));
    if (this.jitter) F.forEach((v, i) => (F[i] = v * this.jitter!()));
    return F;
  }

  /** Static deflection: full coefficient vector, zeros on the constrained ones. */
  solve(load: Load | Float64Array): Float64Array {
    const F = load instanceof Float64Array ? load : this.load(load);
    this.LK ??= cholesky(this.K);
    const c = new Float64Array(this.space.n);
    c.set(cholSolve(this.LK, F.subarray(this.lo, this.hi)), this.lo);
    return c;
  }

  /** The `count` lowest modes; vectors are full coefficient vectors, M-orthonormal. */
  modes(count: number): Modes {
    const r = lowestModes(this.K, this.M, count);
    const vectors = r.vectors.map((v) => {
      const c = new Float64Array(this.space.n);
      c.set(v, this.lo);
      return c;
    });
    return { values: r.values, vectors, iterations: r.iterations };
  }

  /** Every discrete eigenvalue λ_h = ω̂_h², ascending. Dense: keep the mesh modest. */
  spectrum(): Float64Array {
    return bandSpectrum(this.K, this.M);
  }

  /** w and its derivatives up to `d` at x. */
  evaluate(c: ArrayLike<number>, x: number, d = 0): Float64Array {
    return evalField(this.space, c, x, d);
  }
}

/**
 * ∫ coef · B^(r)_i B^(r)_j over every element, into a symmetric band of half-
 * width p. One loop serves both matrices: r = 2 is the stiffness, r = 0 the mass.
 */
function assemble(t: Tabulation, coef: Float64Array, r: number): SymBand {
  const s = t.space, P1 = s.p + 1, stride = (t.d + 1) * P1;
  const A = new SymBand(s.n, s.p);
  for (let e = 0; e < s.ne; e++) {
    const f0 = firstDof(s, e);
    for (let iq = 0; iq < t.nq; iq++) {
      const i = e * t.nq + iq, wq = t.w[i] * coef[i], o = i * stride + r * P1;
      for (let a = 0; a < P1; a++) {
        const Ba = wq * t.B[o + a];
        for (let b = 0; b <= a; b++) A.add(f0 + a, f0 + b, Ba * t.B[o + b]);
      }
    }
  }
  return A;
}

/** 1 + δ, δ uniform in ±size, from a seeded xorshift — reproducible run to run. */
export function jitter(size: number, seed: number): () => number {
  let x = (seed * 0x9e3779b1) >>> 0 || 1;
  return () => {
    x ^= x << 13; x ^= x >>> 17; x ^= x << 5; x >>>= 0;
    return 1 + size * (2 * (x / 0x100000000) - 1);
  };
}
