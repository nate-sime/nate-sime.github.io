// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Symmetric banded matrices, stored as their lower band in one flat array:
 *
 *   data[i·(bw + 1) + d] = A[i][i − d],   0 ≤ d ≤ bw.
 *
 * A 1D spline space of degree p couples functions at most p indices apart,
 * whatever the continuity, so its stiffness and mass matrices have half-
 * bandwidth p and Cholesky costs O(n p²) — the whole reason Monte Carlo over
 * thousands of beam solves is cheap. The layout is the one a batched GPU
 * factorisation would read: one sample per invocation, one flat buffer each.
 */

export class SymBand {
  readonly data: Float64Array;
  constructor(readonly n: number, readonly bw: number, data?: Float64Array) {
    this.data = data ?? new Float64Array(n * (bw + 1));
  }

  get(i: number, j: number): number {
    if (i < j) [i, j] = [j, i];
    return i - j > this.bw ? 0 : this.data[i * (this.bw + 1) + i - j];
  }

  /** A[i][j] += v (and by symmetry A[j][i]); |i − j| must lie inside the band. */
  add(i: number, j: number, v: number): void {
    if (i < j) [i, j] = [j, i];
    this.data[i * (this.bw + 1) + i - j] += v;
  }

  matvec(x: ArrayLike<number>, y = new Float64Array(this.n)): Float64Array {
    const { n, bw, data } = this, W = bw + 1;
    y.fill(0);
    for (let i = 0; i < n; i++) {
      y[i] += data[i * W] * x[i];
      for (let d = 1; d <= bw && d <= i; d++) {
        const a = data[i * W + d];
        y[i] += a * x[i - d];
        y[i - d] += a * x[i];
      }
    }
    return y;
  }

  /** |A| x — what bounds the round-off in A x, entry by entry. */
  absMatvec(x: ArrayLike<number>): Float64Array {
    const abs = new SymBand(this.n, this.bw, this.data.map(Math.abs));
    return abs.matvec(Array.from(x, Math.abs));
  }

  /** The principal submatrix on rows and columns [lo, hi). */
  slice(lo: number, hi: number): SymBand {
    const W = this.bw + 1;
    return new SymBand(hi - lo, this.bw, this.data.slice(lo * W, hi * W).map((v, idx) => {
      // Entries reaching left of `lo` belong to dropped columns.
      const i = Math.floor(idx / W), d = idx % W;
      return i - d < 0 ? 0 : v;
    }));
  }

  toDense(): Float64Array[] {
    const A = Array.from({ length: this.n }, () => new Float64Array(this.n));
    for (let i = 0; i < this.n; i++)
      for (let d = 0; d <= this.bw && d <= i; d++) A[i][i - d] = A[i - d][i] = this.data[i * (this.bw + 1) + d];
    return A;
  }
}

/** Banded Cholesky A = L Lᵀ; L is returned in the same lower-band layout. */
export function cholesky(A: SymBand): SymBand {
  const { n, bw } = A, W = bw + 1;
  const L = new SymBand(n, bw, A.data.slice());
  const l = L.data;
  for (let i = 0; i < n; i++) {
    const j0 = Math.max(0, i - bw);
    for (let j = j0; j <= i; j++) {
      let s = l[i * W + i - j];
      for (let k = Math.max(j0, j - bw); k < j; k++) s -= l[i * W + i - k] * l[j * W + j - k];
      if (i === j) {
        if (!(s > 0)) throw new Error(`matrix is not positive definite (pivot ${i}: ${s})`);
        l[i * W] = Math.sqrt(s);
      } else {
        l[i * W + i - j] = s / l[j * W];
      }
    }
  }
  return L;
}

/** y = L⁻¹ b. */
export function forwardSolve(L: SymBand, b: ArrayLike<number>): Float64Array {
  const { n, bw, data } = L, W = bw + 1, y = Float64Array.from(b);
  for (let i = 0; i < n; i++) {
    let s = y[i];
    for (let k = Math.max(0, i - bw); k < i; k++) s -= data[i * W + i - k] * y[k];
    y[i] = s / data[i * W];
  }
  return y;
}

/** x = L⁻ᵀ y. */
export function backSolve(L: SymBand, y: ArrayLike<number>): Float64Array {
  const { n, bw, data } = L, W = bw + 1, x = Float64Array.from(y);
  for (let i = n - 1; i >= 0; i--) {
    let s = x[i];
    for (let k = i + 1; k <= Math.min(n - 1, i + bw); k++) s -= data[k * W + k - i] * x[k];
    x[i] = s / data[i * W];
  }
  return x;
}

/** x = A⁻¹ b from A's Cholesky factor. */
export const cholSolve = (L: SymBand, b: ArrayLike<number>): Float64Array => backSolve(L, forwardSolve(L, b));
