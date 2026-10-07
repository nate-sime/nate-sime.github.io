// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import { describe, expect, it } from "vitest";
import { SymBand, cholSolve, cholesky } from "../src/band";
import { bandSpectrum, denseGeneralizedEig, lowestModes, symEig } from "../src/eig";

/** A reproducible SPD band: diagonally dominant with random off-diagonals. */
function randomBand(n: number, bw: number, seed: number, shift = 0): SymBand {
  let s = seed;
  const rnd = () => ((s = (s * 1103515245 + 12345) % 2147483648) / 2147483648) - 0.5;
  const A = new SymBand(n, bw);
  for (let i = 0; i < n; i++) {
    for (let d = 1; d <= bw && d <= i; d++) A.add(i, i - d, rnd());
    A.add(i, i, 2 * bw + 1 + shift + rnd());
  }
  return A;
}

const norm = (v: ArrayLike<number>) => Math.sqrt(Array.from(v).reduce((a, b) => a + b * b, 0));

describe("banded Cholesky", () => {
  it("solves A x = b to round-off", () => {
    for (const [n, bw] of [[1, 0], [7, 2], [40, 3], [120, 6]]) {
      const A = randomBand(n, bw, 7 + n);
      const x = Float64Array.from({ length: n }, (_, i) => Math.sin(i + 1));
      const y = cholSolve(cholesky(A), A.matvec(x));
      expect(norm(y.map((v, i) => v - x[i]))).toBeLessThan(1e-12 * norm(x));
    }
  });

  it("refuses an indefinite matrix", () => {
    const A = new SymBand(2, 1);
    A.add(0, 0, 1); A.add(1, 1, 1); A.add(1, 0, 2);
    expect(() => cholesky(A)).toThrow(/positive definite/);
  });

  it("slices a principal submatrix", () => {
    const A = randomBand(10, 3, 3), S = A.slice(2, 8), D = A.toDense();
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) expect(S.get(i, j)).toBe(D[i + 2][j + 2]);
  });
});

describe("symmetric eigensolvers", () => {
  it("diagonalises a dense symmetric matrix", () => {
    for (const n of [1, 2, 5, 30]) {
      const A = randomBand(n, n - 1, 11 * n, -n).toDense();
      const { values, vectors } = symEig(A) as Required<ReturnType<typeof symEig>>;
      for (let i = 0; i < n; i++) {
        const Av = A.map((row) => row.reduce((s, a, j) => s + a * vectors[i][j], 0));
        expect(norm(Av.map((v, j) => v - values[i] * vectors[i][j]))).toBeLessThan(1e-12 * n);
        for (let j = 0; j <= i; j++) {
          const d = vectors[i].reduce((s, v, k) => s + v * vectors[j][k], 0);
          expect(d).toBeCloseTo(i === j ? 1 : 0, 12);
        }
        if (i > 0) expect(values[i]).toBeGreaterThanOrEqual(values[i - 1]);
      }
      // The values-only path skips accumulation but must agree.
      const only = symEig(A, false).values;
      only.forEach((v, i) => expect(v).toBeCloseTo(values[i], 11));
    }
  });

  it("solves the banded generalized problem three ways in agreement", () => {
    // A spread spectrum, as a beam's is (λ_i ~ i⁴): subspace iteration converges
    // at λ_nev / λ_m, which a bunched random spectrum would hold near 1.
    const n = 60, K = randomBand(n, 3, 5), M = randomBand(n, 3, 9, 4);
    for (let i = 0; i < n; i++) K.add(i, i, 0.01 * i ** 3);
    const all = bandSpectrum(K, M);
    const dense = denseGeneralizedEig(K.toDense(), M.toDense());
    const low = lowestModes(K, M, 5);
    for (let i = 0; i < n; i++) expect(all[i]).toBeCloseTo(dense.values[i], 11);
    for (let i = 0; i < 5; i++) {
      expect(low.values[i]).toBeCloseTo(all[i], 11);
      const x = low.vectors[i], Kx = K.matvec(x), Mx = M.matvec(x);
      expect(norm(Kx.map((v, j) => v - low.values[i] * Mx[j]))).toBeLessThan(1e-9 * norm(Kx));
      for (let j = 0; j <= i; j++) {
        const mij = low.vectors[j].reduce((s, v, k) => s + v * Mx[k], 0);
        expect(mij).toBeCloseTo(i === j ? 1 : 0, 10);
      }
    }
  });
});
