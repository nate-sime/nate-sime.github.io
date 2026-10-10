// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { FieldAt, draw, lognormalQuantile, realise, type FieldSpec } from "../src/random/field";
import {
  captured, correlation, exponentialEigenvalues, karhunenLoeve, modesAt, termsOf,
} from "../src/random/kl";
import { CHANNEL, normals, philox4x32, uniforms } from "../src/random/philox";

const hex = (a: Uint32Array) => Array.from(a, (v) => v.toString(16).padStart(8, "0"));

describe("Philox4x32-10", () => {
  it("reproduces the Random123 known-answer vectors", () => {
    const o = new Uint32Array(4), F = 0xffffffff;
    expect(hex(philox4x32(0, 0, 0, 0, 0, 0, o))).toEqual(["6627e8d5", "e169c58d", "bc57ac4c", "9b00dbd8"]);
    expect(hex(philox4x32(F, F, F, F, F, F, o))).toEqual(["408f276d", "41c83b0e", "a20bc7c6", "6d5451fd"]);
    expect(hex(philox4x32(0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344, 0xa4093822, 0x299f31d0, o)))
      .toEqual(["d16cfe09", "94fdcceb", "5001e420", "24126ea1"]);
  });

  it("addresses normals by (seed, sample, channel, j), whatever else is drawn", () => {
    const a = normals(7, 123, CHANNEL.stiffness, 10), b = normals(7, 123, CHANNEL.stiffness, 3);
    expect(Array.from(b)).toEqual(Array.from(a.subarray(0, 3)));
    expect(normals(7, 124, 0, 1)[0]).not.toBe(a[0]);
    expect(normals(8, 123, 0, 1)[0]).not.toBe(a[0]);
    expect(normals(7, 123, CHANNEL.load, 1)[0]).not.toBe(a[0]);
  });

  it("gives standard normals and uniforms, uncorrelated across samples and components", () => {
    const n = 40000;
    let s = 0, s2 = 0, s4 = 0, cross = 0, lag = 0, us = 0;
    let prev = 0;
    for (let i = 0; i < n; i++) {
      const z = normals(1, i, 0, 2), u = uniforms(1, i, 1, 1)[0];
      s += z[0]; s2 += z[0] ** 2; s4 += z[0] ** 4; cross += z[0] * z[1]; lag += z[0] * prev; us += u;
      prev = z[0];
    }
    const tol = 4 / Math.sqrt(n); // four standard errors
    expect(Math.abs(s / n)).toBeLessThan(tol);
    expect(Math.abs(s2 / n - 1)).toBeLessThan(tol * Math.SQRT2);
    expect(Math.abs(s4 / n - 3)).toBeLessThan(tol * Math.sqrt(96));
    expect(Math.abs(cross / n)).toBeLessThan(tol);
    expect(Math.abs(lag / n)).toBeLessThan(tol);
    expect(Math.abs(us / n - 0.5)).toBeLessThan(tol / Math.sqrt(12));
  });
});

describe("Karhunen–Loève by Nyström", () => {
  it("matches the exponential kernel's exact eigenvalues", () => {
    for (const ell of [0.1, 0.5, 2]) {
      const kl = karhunenLoeve("exponential", ell, 64);
      const ex = exponentialEigenvalues(ell, 2000), N = kl.x.length;
      // Unit variance on a unit interval: the trace is 1, less a tail past J
      // terms of ~ 2/(ℓπ²J), since λ_j → 2/(ℓ(jπ)²).
      const tail = 2 / (ell * Math.PI ** 2 * 2000);
      expect(Math.abs(ex.reduce((a, b) => a + b, 0) + tail - 1)).toBeLessThan(0.05 * tail);
      expect(kl.spectrum.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 12);
      // Second order in the node spacing, mode by mode — the kernel's kink.
      for (let j = 0; j < 64; j++)
        expect(Math.abs(kl.values[j] / ex[j] - 1)).toBeLessThan(1.2 * ((j + 1) / N) ** 2 + 2e-4);
    }
  });

  it("decays as the kernel's smoothness says", () => {
    const ell = 0.2, at = (k: Parameters<typeof karhunenLoeve>[0]) => karhunenLoeve(k, ell, 64).spectrum;
    const slope = (s: Float64Array, a: number, b: number) => Math.log(s[b] / s[a]) / Math.log((b + 1) / (a + 1));
    expect(slope(at("exponential"), 20, 40)).toBeCloseTo(-2, 0);
    expect(slope(at("matern32"), 20, 40)).toBeCloseTo(-4, 0);
    // Squared exponential: faster than any power, so the floor truncates it.
    const g = karhunenLoeve("gaussian", ell, 64);
    expect(termsOf(g)).toBeLessThan(40);
    expect(captured(g)).toBeCloseTo(1, 10);
  });

  it("gives orthonormal modes, and an interpolant that agrees with them at the nodes", () => {
    const kl = karhunenLoeve("matern32", 0.3, 12), N = kl.x.length, M = termsOf(kl);
    for (let a = 0; a < M; a++)
      for (let b = 0; b <= a; b++) {
        let s = 0;
        for (let i = 0; i < N; i++) s += kl.w[i] * kl.vectors[a * N + i] * kl.vectors[b * N + i];
        expect(s).toBeCloseTo(a === b ? 1 : 0, 10);
      }
    const at = modesAt(kl, kl.x.subarray(0, 50));
    for (let i = 0; i < 50; i++)
      for (let j = 0; j < M; j++)
        expect(at[i * M + j]).toBeCloseTo(Math.sqrt(kl.values[j]) * kl.vectors[j * N + i], 10);
  });

  it("reproduces the covariance from its modes", () => {
    const kl = karhunenLoeve("exponential", 0.3, 64), C = correlation("exponential", 0.3);
    const pts = [0.1, 0.35, 0.5, 0.9], f = new FieldAt(kl, pts);
    for (let a = 0; a < pts.length; a++)
      for (let b = 0; b < pts.length; b++) {
        let c = 0;
        for (let j = 0; j < f.M; j++) c += f.phi[a * f.M + j] * f.phi[b * f.M + j];
        // 64 terms of a j⁻² series: the tail is ~ 2ℓ/(π² · 64) ≈ 10⁻³ on the diagonal.
        expect(Math.abs(c - C(Math.abs(pts[a] - pts[b])))).toBeLessThan(a === b ? 2e-2 : 3e-3);
      }
  });
});

describe("random inputs", () => {
  const spec: FieldSpec = { kernel: "matern32", ell: 0.25, sigma: 0.5, terms: 20, massFollows: true, loadSigma: 0.3 };
  const kl = karhunenLoeve(spec.kernel, spec.ell, spec.terms);
  const M = termsOf(kl);

  it("is the same function of ω at every set of points — the coupling across levels", () => {
    const coarse = new FieldAt(kl, [0.125, 0.5, 0.875]), fine = new FieldAt(kl, [0.5, 0.0625, 0.125, 0.875]);
    const ones3 = new Float64Array(3).fill(1), ones4 = new Float64Array(4).fill(1);
    for (const index of [0, 17, 4242]) {
      const d = draw(spec, M, 3, index);
      const a = realise(coarse, spec, d, ones3, ones3), b = realise(fine, spec, d, ones4, ones4);
      expect(b.e[0]).toBeCloseTo(a.e[1], 13);
      expect(b.e[2]).toBeCloseTo(a.e[0], 13);
      expect(b.mu[3]).toBeCloseTo(a.mu[2], 13);
      expect(b.q![3]).toBeCloseTo(a.q![2], 13);
      expect(Math.abs(b.mu[0] ** 3 - b.e[0])).toBeLessThan(1e-13);
    }
  });

  it("has mean exactly e₀ at every point, and the lognormal quantiles it claims", () => {
    const pts = [0.05, 0.5, 0.95], f = new FieldAt(kl, pts), ones = new Float64Array(3).fill(1);
    const n = 20000, sum = [0, 0, 0], below = [0, 0, 0], q90 = lognormalQuantile(f, spec.sigma, 1.2815515655);
    for (let i = 0; i < n; i++) {
      const r = realise(f, spec, draw(spec, M, 11, i), ones, ones);
      r.e.forEach((v, q) => { sum[q] += v; if (v < q90[q]) below[q]++; });
    }
    for (let q = 0; q < 3; q++) {
      const cov = Math.sqrt(Math.expm1(spec.sigma ** 2 * f.s[q]));
      expect(Math.abs(sum[q] / n - 1)).toBeLessThan((4 * cov) / Math.sqrt(n));
      expect(Math.abs(below[q] / n - 0.9)).toBeLessThan(4 * Math.sqrt(0.09 / n));
      // Truncation keeps most, never more, of the variance.
      expect(f.s[q]).toBeLessThanOrEqual(1 + 1e-9);
      expect(f.s[q]).toBeGreaterThan(0.95);
    }
  });
});
