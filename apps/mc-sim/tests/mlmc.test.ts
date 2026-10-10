// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { SUPPORTS, freeDofs } from "../src/beam/beam";
import { exactStaticOf, runHierarchy, type BeamCase } from "../src/hierarchy";
import { Mlmc, MlmcSweep, allocate, slope, survey, type MlmcOptions } from "../src/mc/mlmc";
import { Sampler, type McSpec } from "../src/mc/sampler";
import { Accumulator } from "../src/mc/stats";
import type { FieldSpec } from "../src/random/field";
import { karhunenLoeve, type KLData } from "../src/random/kl";

/** A perfectly correlated field, as in montecarlo.test.ts: e = exp(σξ − σ²/2), one number per sample. */
const constantKL: KLData = {
  kernel: "gaussian", ell: 1e9,
  x: Float64Array.of(0.5), w: Float64Array.of(1),
  spectrum: Float64Array.of(1), values: Float64Array.of(1), vectors: Float64Array.of(1),
};

const field = (o: Partial<FieldSpec> = {}): FieldSpec =>
  ({ kernel: "gaussian", ell: 1e9, sigma: 0.4, terms: 1, massFollows: false, loadSigma: 0, ...o });

const pinned: BeamCase = { supports: "pinned–pinned", load: "uniform", section: "uniform" };
const spec = (o: Partial<McSpec> = {}): McSpec =>
  ({ p: 2, k: 1, beam: pinned, qoi: "compliance", field: field(), seed: 1, ...o });

/**
 * Level stores filled on demand, synchronously: what the worker pool does in
 * the app, here in one thread. Level ℓ is ne₀·2^ℓ elements, its parent ne/2,
 * and its samples come from stream ℓ.
 */
class Levels {
  readonly acc: Accumulator[] = [];
  constructor(readonly sampler: Sampler, readonly ne0: number) {}

  ensure(l: number, n: number): Accumulator {
    while (this.acc.length <= l) this.acc.push(new Accumulator(0));
    const a = this.acc[l], from = a.n, count = n - from;
    if (count <= 0) return a;
    const Q = new Float64Array(count), Qc = new Float64Array(count);
    for (let i = 0; i < count; i++) {
      const r = this.sampler.sample(from + i, this.ne0 * 2 ** l, l > 0, false, l);
      Q[i] = r.Q; Qc[i] = r.Qc;
    }
    a.push({ from, count, Q, Qc, W: null, ms: 0 });
    return a;
  }

  run(opts: MlmcOptions): Mlmc {
    const alg = new Mlmc(opts);
    while (!alg.done) {
      alg.wants().forEach((n, l) => this.ensure(l, n));
      alg.update((l, n) => {
        const y = this.acc[l].prefix(n).y;
        return { mean: y.mean, variance: y.variance };
      });
    }
    return alg;
  }
}

const dofCost = (s: McSpec, ne0: number) => (l: number) => {
  const d = (ne: number) => freeDofs({ p: s.p, k: s.k, ne, supports: SUPPORTS[s.beam.supports] });
  return d(ne0 * 2 ** l) + (l > 0 ? d(ne0 * 2 ** (l - 1)) : 0);
};

describe("the optimal allocation", () => {
  it("meets the variance budget at the least cost", () => {
    const V = [1, 0.1, 0.01, 1e-3], C = [1, 2, 4, 8], budget = 1e-4;
    const N = allocate(V, C, budget);
    expect(V.reduce((s, v, l) => s + v / N[l], 0)).toBeLessThanOrEqual(budget * (1 + 1e-12));
    // Any other split meeting the same budget costs more: move variance between two levels.
    const cost = (n: number[]) => n.reduce((s, v, l) => s + v * C[l], 0);
    for (const [a, b] of [[0, 1], [1, 3], [0, 3]]) {
      const n = N.slice();
      n[a] = Math.ceil(n[a] * 1.3);
      // Shrink level b only as far as the budget allows.
      const room = budget - V.reduce((s, v, l) => s + (l === b ? 0 : v / n[l]), 0);
      n[b] = Math.ceil(V[b] / room);
      expect(cost(n)).toBeGreaterThan(cost(N) * 0.999);
    }
  });

  it("fits a slope on levels 1 … L only", () => {
    expect(slope([123, 2 ** -3, 2 ** -7, 2 ** -11])).toBeCloseTo(-4, 12);
  });
});

describe("the survey against closed forms", () => {
  it("measures β = 2α of the deterministic hierarchy, as a perfectly correlated field predicts", () => {
    // Compliance ∝ 1/e, so Y_ℓ = ΔQ_ℓ · e^{−σξ + σ²/2} with ΔQ_ℓ the stage 3
    // successive difference: V[Y_ℓ] = ΔQ_ℓ² e^{2σ²}(e^{σ²} − 1), exactly.
    const sigma = 0.4, N = 4000, ne0 = 4;
    for (const p of [2, 3]) {
      const s = spec({ p, k: p - 1 });
      const lv = new Levels(new Sampler(s, constantKL), ne0);
      const S = survey([0, 1, 2, 3, 4].map((l) => {
        const a = lv.ensure(l, N);
        return { fine: a.values, coarse: a.coarseValues };
      }), N);
      const H = runHierarchy({ p, k: p - 1, beam: pinned, qoi: "compliance", ne0, levels: 5 });
      const factor = Math.exp(2 * sigma ** 2) * Math.expm1(sigma ** 2);
      for (let l = 1; l <= 4; l++) {
        const ratio = S[l].varY / (H.levels[l].dQ ** 2 * factor);
        // A sample variance over N lognormal draws: within ~4 standard errors.
        expect(Math.abs(ratio - 1)).toBeLessThan(4 * Math.sqrt((S[l].kurtosis - 1) / N));
        expect(S[l].consistency).toBeLessThan(1.5);
      }
      const beta = -slope(S.map((v) => v.varY)), alpha = -slope(S.map((v) => Math.abs(v.meanY)));
      expect(beta / 2).toBeCloseTo(H.alpha!, 0);
      expect(Math.abs(beta - 2 * H.alpha!)).toBeLessThan(0.15);
      expect(Math.abs(alpha - H.alpha!)).toBeLessThan(0.1);
    }
  });

  it("finds V[Y_ℓ] far below V[Q_ℓ] and falling, once the mesh resolves the field", () => {
    const f = field({ kernel: "matern32", ell: 0.3, sigma: 0.5, terms: 32 });
    const s = spec({ p: 3, k: 2, qoi: "deflection", beam: { ...pinned, supports: "cantilever" }, field: f });
    const lv = new Levels(new Sampler(s, karhunenLoeve(f.kernel, f.ell, f.terms)), 4);
    const N = 400;
    const S = survey([0, 1, 2, 3, 4].map((l) => { const a = lv.ensure(l, N); return { fine: a.values, coarse: a.coarseValues }; }), N);
    expect(S[1].varY / S[1].varQ).toBeLessThan(1e-2);
    for (let l = 2; l <= 4; l++) expect(S[l].varY).toBeLessThan(S[l - 1].varY / 4);
    for (let l = 1; l <= 4; l++) expect(S[l].consistency).toBeLessThan(1.5);
  });
});

describe("multilevel Monte Carlo", () => {
  const sigma = 0.4, ne0 = 4;
  // Pinned–pinned, uniform load: compliance 1/120 for the uniform beam, and
  // E[1/e] = e^{σ²} for the lognormal factor.
  const truth = exactStaticOf(pinned)!.compliance * Math.exp(sigma ** 2);
  const opts = (s: McSpec, eps: number): MlmcOptions =>
    ({ eps, N0: 200, Lmin: 2, Lmax: 8, cost: dofCost(s, ne0) });

  it("is unbiased for the finest level: the telescoping sum adds up", () => {
    const s = spec({ p: 2, k: 1 });
    const lv = new Levels(new Sampler(s, constantKL), ne0);
    const alg = lv.run(opts(s, 1e-3 * truth));
    expect(alg.status).toBe("converged");
    expect(alg.L).toBeGreaterThanOrEqual(4); // measured: L = 4, N = 231620, 5164, 913, 161, 29
    // E[Q_L] in closed form on the discrete mesh: the deterministic Q_L times e^{σ²}.
    const QL = runHierarchy({ p: 2, k: 1, beam: pinned, qoi: "compliance", ne0, levels: alg.L + 1 }).levels[alg.L].Q;
    expect(Math.abs(alg.estimate - QL * Math.exp(sigma ** 2))).toBeLessThan(4 * Math.sqrt(alg.sampleVariance));
    // And the sampling variance is within the (1 − θ)ε² budget it was allocated.
    expect(alg.sampleVariance).toBeLessThan(1.05 * (1 - alg.theta) * (1e-3 * truth) ** 2);
  });

  it("delivers a root-mean-square error within ε, over independent runs", () => {
    const eps = 4e-3 * truth, runs = 32;
    let se = 0;
    const Ls: number[] = [];
    for (let seed = 1; seed <= runs; seed++) {
      const s = spec({ p: 2, k: 1, seed });
      const alg = new Levels(new Sampler(s, constantKL), ne0).run(opts(s, eps));
      se += (alg.estimate - truth) ** 2;
      Ls.push(alg.L);
    }
    const rmse = Math.sqrt(se / runs);
    // χ² with 32 degrees of freedom: the measured RMSE scatters ±13% about the
    // true one. Measured: 1.08 ε, every run settling on L = 3.
    expect(rmse).toBeLessThan(1.3 * eps);
    expect(Math.max(...Ls)).toBeGreaterThan(2); // the bias test did add levels
  });

  it("gives the same answer from a shared store as from a fresh one, whatever else read it", () => {
    const s = spec({ p: 2, k: 1, seed: 5 });
    const shared = new Levels(new Sampler(s, constantKL), ne0);
    const coarse = shared.run(opts(s, 4e-3 * truth));
    const fine = shared.run(opts(s, 1e-3 * truth));
    const alone = new Levels(new Sampler(s, constantKL), ne0).run(opts(s, 1e-3 * truth));
    expect(fine.N).toEqual(alone.N);
    expect(fine.estimate).toBe(alone.estimate);
    expect(coarse.cost).toBeLessThan(fine.cost);
  });

  it("costs O(ε⁻²) when the variance decays faster than the cost grows", () => {
    const s = spec({ p: 2, k: 1, seed: 3 });
    const lv = new Levels(new Sampler(s, constantKL), ne0);
    const work = [4e-3, 2e-3, 1e-3].map((e) => { const a = lv.run(opts(s, e * truth)); return a.cost * (e * truth) ** 2; });
    // β = 4 > γ = 1: ε²·cost settles to a constant (up to the integer N_ℓ and the
    // levels added). Measured: 1.06, 1.01, 1.02, with L = 3, 3, 4.
    expect(Math.max(...work) / Math.min(...work)).toBeLessThan(2);
  });

  it("runs a sweep of tolerances from one store, and stops any that would overrun its budget", () => {
    const s = spec({ p: 2, k: 1, seed: 5 });
    const lv = new Levels(new Sampler(s, constantKL), ne0);
    const eps = [4e-3, 1e-3, 1e-4].map((e) => e * truth);
    const sweep = new MlmcSweep(6, 300, eps, { N0: 200, Lmin: 2, cost: dofCost(s, ne0) }, 5e5);
    // As the app drives it: make what it asks for, in no particular order, until it asks for nothing new.
    for (let guard = 0; guard < 100 && !sweep.done; guard++) {
      const demand = sweep.step((l) => lv.ensure(l, 0));
      [...demand.keys()].reverse().forEach((l) => lv.ensure(l, demand[l]));
    }
    expect(sweep.runs.map((a) => a.status)).toEqual(["converged", "converged", "over budget"]);
    // ε = 10⁻⁴ wants ~2·10⁷ samples on level 0; the sweep stopped it before asking for them.
    expect(lv.acc[0].n).toBeLessThanOrEqual(5e5);
    // Each tolerance landed exactly where it would have alone.
    const alone = new Levels(new Sampler(s, constantKL), ne0).run(opts(s, eps[1]));
    expect(sweep.runs[1].N).toEqual(alone.N);
    expect(sweep.runs[1].estimate).toBe(alone.estimate);
  });
});
