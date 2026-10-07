// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { FIELD_POINTS, Sampler, type McSpec } from "../src/mc/sampler";
import { Accumulator, Z95, histogram, type Batch } from "../src/mc/stats";
import type { FieldSpec } from "../src/random/field";
import { karhunenLoeve, type KLData } from "../src/random/kl";

/**
 * A perfectly correlated field: one mode, φ ≡ 1, λ = 1 (a squared-exponential
 * kernel of enormous length, which the Nyström interpolant turns into a
 * constant). Then e = exp(σξ − σ²/2) is one number per sample, K scales by it
 * exactly — in the discrete problem as in the continuous one — and every
 * moment of Q is known in closed form on any mesh.
 */
const constantKL: KLData = {
  kernel: "gaussian", ell: 1e9,
  x: Float64Array.of(0.5), w: Float64Array.of(1),
  spectrum: Float64Array.of(1), values: Float64Array.of(1), vectors: Float64Array.of(1),
};

const field = (o: Partial<FieldSpec> = {}): FieldSpec =>
  ({ kernel: "gaussian", ell: 1e9, sigma: 0.4, terms: 1, massFollows: false, loadSigma: 0, ...o });

const spec = (o: Partial<McSpec> = {}): McSpec => ({
  p: 4, k: 3, beam: { supports: "pinned–pinned", load: "uniform", section: "uniform" },
  qoi: "deflection", field: field(), seed: 1, ...o,
});

/** Mean and standard error of Q over samples [0, n). */
function run(s: Sampler, n: number, ne: number) {
  const acc = new Accumulator(FIELD_POINTS);
  const Q = new Float64Array(n), Qc = new Float64Array(n).fill(NaN);
  for (let i = 0; i < n; i++) Q[i] = s.solve(i, ne).Q;
  acc.push({ from: 0, count: n, Q, Qc, W: null, ms: 0 });
  return acc;
}

describe("plain Monte Carlo against closed forms", () => {
  const sigma = 0.4;

  it("deflection: E[w] = w₀ E[1/e] = w₀ e^{σ²}, with the lognormal variance", () => {
    // Quartic splines hold the uniform-load solution exactly: no bias at all.
    const acc = run(new Sampler(spec(), constantKL), 4000, 4);
    const w0 = 5 / 384, mean = w0 * Math.exp(sigma ** 2);
    const sd = w0 * Math.sqrt(Math.expm1(sigma ** 2)) * Math.exp(sigma ** 2);
    expect(Math.abs(acc.q.mean - mean)).toBeLessThan(4 * acc.q.se);
    expect(acc.q.sd / sd).toBeCloseTo(1, 1);
  });

  it("ω₁: E = ω₁ₕ e^{−σ²/8}, and e^{−σ²/9} when the mass follows the depth", () => {
    const det = new Sampler(spec({ qoi: "omega1", field: field({ sigma: 0 }) }), constantKL).solve(0, 8).Q;
    const a = run(new Sampler(spec({ qoi: "omega1" }), constantKL), 3000, 8);
    expect(Math.abs(a.q.mean - det * Math.exp(-(sigma ** 2) / 8))).toBeLessThan(4 * a.q.se);
    const b = run(new Sampler(spec({ qoi: "omega1", field: field({ massFollows: true }) }), constantKL), 3000, 8);
    expect(Math.abs(b.q.mean - det * Math.exp(-(sigma ** 2) / 9))).toBeLessThan(4 * b.q.se);
    // Same ω, same ξ, sample by sample: ω ∝ e^{1/2} with a random modulus,
    // e^{1/3} with a random depth, so both recover the same e.
    for (const i of [0, 5, 99]) expect((b.values[i] / det) ** 3).toBeCloseTo((a.values[i] / det) ** 2, 10);
  });

  it("a random load adds its own variance and leaves the mean alone", () => {
    const s = new Sampler(spec({ field: field({ sigma: 0, loadSigma: 0.5 }) }), constantKL);
    const acc = run(s, 4000, 4);
    // q = 1 + 0.5 ξ′ uniformly, so w = (1 + 0.5 ξ′) w₀.
    expect(Math.abs(acc.q.mean - 5 / 384)).toBeLessThan(4 * acc.q.se);
    expect(acc.q.sd / (0.5 * 5 / 384)).toBeCloseTo(1, 1);
  });

  it("covers the truth with its 95% interval about 95% of the time", () => {
    const w0 = 5 / 384, mean = w0 * Math.exp(sigma ** 2);
    let hits = 0;
    const runs = 300;
    for (let seed = 1; seed <= runs; seed++) {
      const acc = run(new Sampler(spec({ seed, p: 4 }), constantKL), 60, 2);
      if (Math.abs(acc.q.mean - mean) < Z95 * acc.q.se) hits++;
    }
    // Binomial(300, 0.95) has sd ≈ 3.8 hits; the lognormal's skew costs a
    // little coverage at N = 60, so the band is lopsided toward under-coverage.
    expect(hits / runs).toBeGreaterThan(0.88);
    expect(hits / runs).toBeLessThan(0.99);
  });
});

describe("statistics in sample order", () => {
  it("are the same bits however the batches were cut and in whatever order they came", () => {
    const kl = karhunenLoeve("matern32", 0.3, 16);
    const s = new Sampler(spec({ p: 3, k: 2, field: field({ kernel: "matern32", ell: 0.3, terms: 16 }) }), kl);
    const n = 40, Q = new Float64Array(n), Qc = new Float64Array(n), W = new Float64Array(n * FIELD_POINTS);
    for (let i = 0; i < n; i++) {
      const r = s.sample(i, 8, true, true);
      Q[i] = r.Q; Qc[i] = r.Qc; W.set(r.w!, i * FIELD_POINTS);
    }
    const cut = (bounds: number[]): Batch[] => bounds.slice(0, -1).map((from, j) => {
      const to = bounds[j + 1];
      return { from, count: to - from, Q: Q.slice(from, to), Qc: Qc.slice(from, to), W: W.slice(from * FIELD_POINTS, to * FIELD_POINTS), ms: 1 };
    });
    const a = new Accumulator(FIELD_POINTS), b = new Accumulator(FIELD_POINTS);
    for (const x of cut([0, 40])) a.push(x);
    const pieces = cut([0, 3, 4, 17, 30, 40]);
    for (const x of [pieces[3], pieces[1], pieces[4], pieces[0], pieces[2]]) b.push(x);
    expect(b.n).toBe(40);
    expect(b.q.mean).toBe(a.q.mean);
    expect(b.dq.variance).toBe(a.dq.variance);
    expect(Array.from(b.field.mean)).toEqual(Array.from(a.field.mean));
    expect(b.trajectory).toEqual(a.trajectory);
    // A batch that arrives early waits: nothing past the gap is folded.
    const c = new Accumulator(FIELD_POINTS);
    c.push(pieces[2]);
    expect(c.n).toBe(0);
    expect(c.waiting).toBe(13);
  });

  it("draws a histogram that is a density", () => {
    const v = Float64Array.from({ length: 5000 }, (_, i) => Math.sin(i) ** 3);
    const h = histogram(v, -1, 1, 0.5);
    let area = 0;
    h.density.forEach((d, i) => (area += d * (h.edges[i + 1] - h.edges[i])));
    expect(area).toBeCloseTo(1, 12);
  });
});

describe("coupled levels", () => {
  it("make the correction far less variable than Q, and less so with every level", () => {
    const f = field({ kernel: "matern32", ell: 0.3, sigma: 0.5, terms: 32 });
    const kl = karhunenLoeve(f.kernel, f.ell, f.terms);
    const s = new Sampler(spec({ p: 3, k: 2, beam: { supports: "cantilever", load: "uniform", section: "uniform" }, field: f }), kl);
    const n = 200, V: number[] = [], VQ: number[] = [];
    for (const ne of [8, 16, 32, 64]) {
      const acc = new Accumulator(FIELD_POINTS);
      const Q = new Float64Array(n), Qc = new Float64Array(n);
      for (let i = 0; i < n; i++) { const r = s.sample(i, ne, true); Q[i] = r.Q; Qc[i] = r.Qc; }
      acc.push({ from: 0, count: n, Q, Qc, W: null, ms: 0 });
      V.push(acc.dq.variance);
      VQ.push(acc.q.variance);
    }
    // Measured: V[Y]/V[Q] = 2.7e-4, 7.5e-6, 1.2e-7, 3.9e-10 for ne = 8 … 64.
    expect(V[0] / VQ[0]).toBeLessThan(1e-3);
    for (let l = 1; l < V.length; l++) expect(V[l - 1] / V[l]).toBeGreaterThan(4);
  });
});
