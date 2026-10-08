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

/** A complex vector as two real ones. */
export interface Complex {
  readonly re: Float64Array;
  readonly im: Float64Array;
}

/**
 * x = A⁻¹ b for the complex symmetric (not Hermitian) A = a·K + b·M, a and b
 * complex scalars given as [re, im] — the dynamic stiffness of a harmonically
 * forced, Rayleigh-damped structure, K − Ω²M + iΩ(αM + βK), has exactly this
 * form. Banded LDLᵀ without pivoting, in complex arithmetic: Aᵀ = A, so the
 * factorisation keeps the band and the symmetry, and costs O(n p²) like the
 * real one.
 *
 * Without pivoting, because with damping the imaginary part ΩC is positive
 * definite: then −iA has a positive definite Hermitian part, which is the
 * classical condition for Gaussian elimination to need no pivots and to stay
 * stable (Higham, Accuracy and Stability, §10.4). Undamped at a resonance, A is
 * singular, and no factorisation would help.
 */
export function complexSolve(
  K: SymBand, M: SymBand, a: readonly [number, number], b: readonly [number, number], rhs: ArrayLike<number>,
): Complex {
  const { n, bw } = K, W = bw + 1;
  if (M.n !== n || M.bw !== bw) throw new Error("K and M must share their shape");
  const lr = new Float64Array(n * W), li = new Float64Array(n * W);
  for (let i = 0; i < K.data.length; i++) {
    lr[i] = a[0] * K.data[i] + b[0] * M.data[i];
    li[i] = a[1] * K.data[i] + b[1] * M.data[i];
  }
  // Row i: L[i][j] D[j] (j < i) into the band, then D[i] on the diagonal.
  for (let i = 0; i < n; i++) {
    const j0 = Math.max(0, i - bw);
    for (let j = j0; j <= i; j++) {
      let sr = lr[i * W + i - j], si = li[i * W + i - j];
      for (let k = Math.max(j0, j - bw); k < j; k++) {
        // (L[i][k] D[k]) · L[j][k], with L[i][k] D[k] already stored for k < j ≤ i.
        const ur = lr[i * W + i - k], ui = li[i * W + i - k];
        const vr = lr[j * W + j - k], vi = li[j * W + j - k];
        const dr = lr[k * W], di = li[k * W], d2 = dr * dr + di * di;
        // L[j][k] = (L[j][k] D[k]) / D[k].
        const xr = (vr * dr + vi * di) / d2, xi = (vi * dr - vr * di) / d2;
        sr -= ur * xr - ui * xi;
        si -= ur * xi + ui * xr;
      }
      lr[i * W + i - j] = sr;
      li[i * W + i - j] = si;
    }
    if (lr[i * W] === 0 && li[i * W] === 0) throw new Error(`singular dynamic stiffness (pivot ${i})`);
  }
  // L y = b (unit lower), z = D⁻¹ y, Lᵀ x = z; L[i][k] = stored(i, k) / D[k].
  const xr = Float64Array.from(rhs), xi = new Float64Array(n);
  const ldiv = (i: number, k: number): [number, number] => {
    const ur = lr[i * W + i - k], ui = li[i * W + i - k], dr = lr[k * W], di = li[k * W], d2 = dr * dr + di * di;
    return [(ur * dr + ui * di) / d2, (ui * dr - ur * di) / d2];
  };
  for (let i = 0; i < n; i++)
    for (let k = Math.max(0, i - bw); k < i; k++) {
      const [pr, pi] = ldiv(i, k);
      xr[i] -= pr * xr[k] - pi * xi[k];
      xi[i] -= pr * xi[k] + pi * xr[k];
    }
  for (let i = 0; i < n; i++) {
    const dr = lr[i * W], di = li[i * W], d2 = dr * dr + di * di, r = xr[i], m = xi[i];
    xr[i] = (r * dr + m * di) / d2;
    xi[i] = (m * dr - r * di) / d2;
  }
  for (let i = n - 1; i >= 0; i--)
    for (let k = i + 1; k <= Math.min(n - 1, i + bw); k++) {
      const [pr, pi] = ldiv(k, i);
      xr[i] -= pr * xr[k] - pi * xi[k];
      xi[i] -= pr * xi[k] + pi * xr[k];
    }
  return { re: xr, im: xi };
}
