// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Multilevel Monte Carlo (Giles, Operations Research 2008; Acta Numerica 2015).
 *
 * The finest level's expectation is a telescoping sum over the hierarchy,
 *
 *   E[Q_L] = E[Q_0] + Σ_{ℓ=1}^{L} E[Q_ℓ − Q_{ℓ−1}] = Σ_ℓ E[Y_ℓ],
 *
 * and each term is estimated by its own Monte Carlo average over N_ℓ samples —
 * independent across levels (each level reads its own random stream), but with
 * Q_ℓ and Q_{ℓ−1} computed from the same ω within a level, so Y_ℓ is small and
 * V_ℓ = V[Y_ℓ] smaller still. The estimator's mean-square error splits as
 *
 *   MSE = (E[Q_L] − E[Q])² + Σ_ℓ V_ℓ / N_ℓ,
 *
 * bias² from the finest mesh plus the sampling variance of every level. Given
 * a per-sample cost C_ℓ, the N_ℓ minimising the cost Σ N_ℓ C_ℓ for a variance
 * budget (1 − θ)ε² are (Lagrange)
 *
 *   N_ℓ = ⌈ (1 − θ)⁻¹ ε⁻² √(V_ℓ / C_ℓ) Σ_k √(V_k C_k) ⌉,
 *
 * and the bias, θε² of the budget, sets L. With |E[Y_ℓ]| ~ 2^{−αℓ},
 * V_ℓ ~ 2^{−βℓ} and C_ℓ ~ 2^{γℓ}, the total cost to reach ε is ε⁻² when β > γ
 * — the cost of Monte Carlo on a single mesh of fixed cost, however fine the
 * mesh must be — against ε^{−2−γ/α} for plain Monte Carlo on the mesh fine
 * enough for the bias.
 *
 * `Mlmc` is Giles' adaptive algorithm (the 2015 `mlmc.m`), as a state machine
 * in rounds: `wants()` names the cumulative samples each level needs, the
 * caller makes sure they exist, and `update()` reads the statistics of exactly
 * those prefixes and decides the next round. Nothing depends on how or when the
 * samples were produced, so a run is reproducible from its seed, and one store
 * of samples per level serves every tolerance of a sweep at once.
 *
 * Each round:
 *  1. estimate m_ℓ = |E[Y_ℓ]|, V_ℓ (floored on fine levels, below, while their
 *     few samples could make either look zero);
 *  2. fit α, β by regression on levels 1…L (γ is known: the cost model);
 *  3. set N_ℓ by the formula above;
 *  4. once no level wants more than 1% more, estimate the remaining bias from
 *     the last three corrections as max_i m_{L−i} 2^{−αi} / (2^α − 1); if it is
 *     above √θ ε, add level L + 1 (its variance extrapolated, V_L / 2^β) and go
 *     on; otherwise stop.
 *
 * One departure from Giles' code: step 4 reads corrections only, levels ≥ 1.
 * His reads m_{L−2} even at L = 2, where it is |E[Q₀]| — the quantity itself,
 * not a correction — and so asks for a third level whatever the bias.
 */

export interface MlmcOptions {
  /** Target root-mean-square error, absolute. */
  readonly eps: number;
  /** Share of ε² given to bias²; Giles' default 1/4. */
  readonly theta?: number;
  /** Warm-up samples on each of the first levels. */
  readonly N0: number;
  /** Levels the run starts with are 0 … Lmin; it never goes past Lmax. */
  readonly Lmin: number;
  readonly Lmax: number;
  /** Cost of one sample of level ℓ (fine and coarse solve), in any fixed unit. */
  readonly cost: (level: number) => number;
}

/** What a level's first N samples say: E and V of the correction Y_ℓ. */
export interface LevelMoments {
  readonly mean: number;
  readonly variance: number;
}

/** `failed`: the bias asked for a level past Lmax; `over budget`: some level asked for more samples than the caller allows. */
export type MlmcStatus = "sampling" | "converged" | "failed" | "over budget";

/** One round of the algorithm: what it knew and what it asked for next. */
export interface Round {
  readonly L: number;
  readonly N: readonly number[];
  readonly estimate: number;
  readonly alpha: number;
  readonly beta: number;
}

export const THETA = 0.25;

export class Mlmc {
  readonly theta: number;
  /** Finest level in use. */
  L: number;
  /** Cumulative samples wanted on levels 0 … L. */
  private target: number[];
  /** Samples whose statistics the last update read. */
  N: number[] = [];
  /** |E[Y_ℓ]| and V[Y_ℓ] as last read, floored as in step 1, and the signed means. */
  m: number[] = [];
  V: number[] = [];
  means: number[] = [];
  /** Fitted rates; α ≥ ½ and β ≥ ½ as in Giles' code, so a level is never assumed free. */
  alpha = NaN;
  beta = NaN;
  status: MlmcStatus = "sampling";
  /** The remaining bias estimate when it stopped. */
  bias = NaN;
  readonly rounds: Round[] = [];

  constructor(readonly opts: MlmcOptions) {
    if (opts.Lmin < 2) throw new Error("the rates are fitted on levels 1 … L, so start with at least levels 0, 1, 2");
    this.theta = opts.theta ?? THETA;
    this.L = opts.Lmin;
    this.target = Array.from({ length: this.L + 1 }, () => opts.N0);
  }

  /** Samples level ℓ must have (cumulative) before `update` may be called. */
  wants(): readonly number[] {
    return this.target;
  }

  get done(): boolean {
    return this.status !== "sampling";
  }

  /** Stop where it stands: the caller will not produce what it wants. */
  stop(status: Exclude<MlmcStatus, "sampling" | "converged">): void {
    if (!this.done) this.status = status;
  }

  /** Σ_ℓ E[Y_ℓ]: the estimate of E[Q_L]. */
  get estimate(): number {
    return this.means.reduce((s, v) => s + v, 0);
  }

  /** Σ_ℓ V_ℓ / N_ℓ: the estimator's sampling variance, from the unfloored V_ℓ. */
  sampleVariance = NaN;

  /** Σ_ℓ N_ℓ C_ℓ. */
  get cost(): number {
    return this.N.reduce((s, n, l) => s + n * this.opts.cost(l), 0);
  }

  /**
   * Read the moments of each level's first `wants()[ℓ]` samples and decide the
   * next round. `stats(ℓ, n)` must be the moments of Y_ℓ over samples 0 … n − 1
   * of level ℓ's stream.
   */
  update(stats: (level: number, n: number) => LevelMoments): void {
    if (this.done) return;
    const { opts, theta } = this, eps2 = opts.eps ** 2, L = this.L;
    this.N = this.target.slice();
    const read = this.N.map((n, l) => stats(l, n));
    this.means = read.map((r) => r.mean);
    this.sampleVariance = read.reduce((s, r, l) => s + (r.variance || 0) / this.N[l], 0);
    const m = read.map((r) => Math.abs(r.mean));
    const V = read.map((r) => (Number.isFinite(r.variance) ? Math.max(0, r.variance) : 0));
    const C = this.N.map((_, l) => opts.cost(l));
    this.alpha = Math.max(0.5, -slope(m));
    this.beta = Math.max(0.5, -slope(V));
    // A fine level with few samples can show a mean or a variance of zero (or
    // merely too small): hold each above half the extrapolation of the level
    // below it, so neither the bias test nor the allocation trusts it.
    for (let l = 2; l <= L; l++) {
      m[l] = Math.max(m[l], (0.5 * m[l - 1]) / 2 ** this.alpha);
      V[l] = Math.max(V[l], (0.5 * V[l - 1]) / 2 ** this.beta);
    }
    this.m = m;
    this.V = V;
    this.rounds.push({ L, N: this.N.slice(), estimate: this.estimate, alpha: this.alpha, beta: this.beta });

    let want = allocate(V, C, (1 - theta) * eps2);
    // Nearly there (no level short by more than 1%): is the mesh fine enough?
    if (!want.some((n, l) => n > 1.01 * this.N[l])) {
      const a = this.alpha;
      this.bias = Math.max(...[0, 1, 2].filter((i) => L - i >= 1).map((i) => m[L - i] * 2 ** (-a * i))) / (2 ** a - 1);
      if (this.bias > Math.sqrt(theta) * opts.eps) {
        if (L === opts.Lmax) {
          this.status = "failed";
          return;
        }
        // A new level: no samples yet, its variance extrapolated, its count from the formula.
        this.L = L + 1;
        V.push(V[L] / 2 ** this.beta);
        C.push(opts.cost(L + 1));
        want = allocate(V, C, (1 - theta) * eps2);
      }
    }
    // Never fewer than already read, and at least two on a level, so it has a variance.
    this.target = want.map((n, l) => Math.max(n, this.N[l] ?? 0, 2));
    if (this.target.every((n, l) => n === this.N[l])) this.status = "converged";
  }
}

/** The optimal N_ℓ for a variance budget `budget`: ⌈√(V_ℓ/C_ℓ) Σ√(V_k C_k) / budget⌉. */
export function allocate(V: readonly number[], C: readonly number[], budget: number): number[] {
  const sum = V.reduce((s, v, l) => s + Math.sqrt(v * C[l]), 0);
  return V.map((v, l) => Math.ceil((Math.sqrt(v / C[l]) * sum) / budget));
}

/**
 * Least-squares slope of log₂ x_ℓ against ℓ over levels 1 … L — the rates
 * are about corrections, and level 0 holds Q₀ itself, not a correction.
 */
export function slope(x: readonly number[]): number {
  const pts = x.map((v, l) => [l, Math.log2(v)] as const).filter(([l, y]) => l >= 1 && Number.isFinite(y));
  if (pts.length < 2) return NaN;
  const mx = pts.reduce((s, [l]) => s + l, 0) / pts.length, my = pts.reduce((s, [, y]) => s + y, 0) / pts.length;
  let sxy = 0, sxx = 0;
  for (const [l, y] of pts) { sxy += (l - mx) * (y - my); sxx += (l - mx) ** 2; }
  return sxy / sxx;
}

/**
 * The survey Giles' `mlmc_test` starts with: a fixed number of samples on every
 * level, and what they say about the hierarchy before any tolerance is chosen.
 */
export interface SurveyLevel {
  readonly level: number;
  readonly N: number;
  /** E and V of Q_ℓ (the fine solve of the pair). */
  readonly meanQ: number;
  readonly varQ: number;
  /** E and V of Y_ℓ = Q_ℓ − Q_{ℓ−1} (Q₀ on level 0). */
  readonly meanY: number;
  readonly varY: number;
  /** Kurtosis of Y_ℓ: large means V_ℓ is estimated from a few rare samples. */
  readonly kurtosis: number;
  /**
   * The telescoping identity E[Y_ℓ] = E[Q_ℓ] − E[Q_{ℓ−1}], with E[Q_{ℓ−1}]
   * from level ℓ − 1's own (independent) samples: the discrepancy over three
   * standard errors, as Giles' `mlmc_test` checks it. Well above 1, the coarse
   * solve within a level is not the same function of ω as the level below's
   * fine solve — a coupling bug.
   */
  readonly consistency: number;
}

/** Mean, variance and kurtosis by two passes — exact enough for a fourth moment. */
export function moments4(x: ArrayLike<number>): { mean: number; variance: number; kurtosis: number } {
  const n = x.length;
  let s = 0;
  for (let i = 0; i < n; i++) s += x[i];
  const mean = s / n;
  let m2 = 0, m4 = 0;
  for (let i = 0; i < n; i++) { const d = (x[i] - mean) ** 2; m2 += d; m4 += d * d; }
  return { mean, variance: n > 1 ? m2 / (n - 1) : NaN, kurtosis: m2 > 0 ? (n * m4) / (m2 * m2) : NaN };
}

/** The survey of levels 0 … fine.length − 1, from each level's first N samples: Q_ℓ and Q_{ℓ−1} in sample order. */
export function survey(levels: readonly { readonly fine: ArrayLike<number>; readonly coarse: ArrayLike<number> }[], N: number): SurveyLevel[] {
  const out: SurveyLevel[] = [];
  levels.forEach(({ fine, coarse }, l) => {
    const Qf = Array.prototype.slice.call(fine, 0, N) as number[];
    const Y = l === 0 ? Qf : Qf.map((q, i) => q - coarse[i]);
    const q = moments4(Qf), y = moments4(Y);
    let consistency = NaN;
    if (l > 0) {
      const prev = out[l - 1];
      const se = (Math.sqrt(y.variance) + Math.sqrt(prev.varQ) + Math.sqrt(q.variance)) / Math.sqrt(N);
      consistency = Math.abs(y.mean + prev.meanQ - q.mean) / (3 * se);
    }
    out.push({ level: l, N, meanQ: q.mean, varQ: q.variance, meanY: y.mean, varY: y.variance, kurtosis: y.kurtosis, consistency });
  });
  return out;
}

/** What a sweep reads from a level: how many samples it has, and the moments of Y over any prefix. */
export interface LevelSource {
  readonly n: number;
  prefix(n: number): { readonly y: LevelMoments };
}

/**
 * Giles' `mlmc_test` as one experiment: a survey of `levels` levels at `surveyN`
 * samples each, and the adaptive algorithm at every tolerance in `eps` — all
 * reading one store of samples per level. `step` advances whatever the samples
 * in hand allow and returns what each level should have next; the caller makes
 * the samples (in workers, in any order) and calls it again. Each tolerance
 * reads only its own prefixes, so it lands exactly where it would have alone.
 * A tolerance that would need more than `maxN` samples on some level is
 * stopped "over budget" rather than left to run for hours.
 */
export class MlmcSweep {
  readonly runs: Mlmc[];

  constructor(
    readonly levels: number, readonly surveyN: number, readonly eps: readonly number[],
    base: Omit<MlmcOptions, "eps" | "Lmax">, readonly maxN = Infinity,
  ) {
    this.runs = eps.map((e) => new Mlmc({ ...base, eps: e, Lmax: levels - 1 }));
  }

  step(source: (level: number) => LevelSource): number[] {
    for (const alg of this.runs)
      while (!alg.done) {
        if (alg.wants().some((n) => n > this.maxN)) alg.stop("over budget");
        else if (alg.wants().every((n, l) => source(l).n >= n)) alg.update((l, n) => source(l).prefix(n).y);
        else break;
      }
    const demand = Array.from({ length: this.levels }, () => this.surveyN);
    for (const alg of this.runs)
      if (!alg.done) alg.wants().forEach((n, l) => (demand[l] = Math.max(demand[l], n)));
    return demand;
  }

  get done(): boolean {
    return this.runs.every((a) => a.done);
  }
}
