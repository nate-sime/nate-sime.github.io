// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Symmetric eigensolvers for the vibration problem K φ = λ M φ.
 *
 * Two jobs, two methods:
 *
 * - the **whole spectrum** (`Beam.spectrum`, which the tests compare with the
 *   exact one, mode by mode) — Householder tridiagonalisation and implicit
 *   QL on the dense standard form L⁻¹ K L⁻ᵀ, M = L Lᵀ. O(n³), so it is for
 *   meshes a reader can count, not for sampling.
 * - the **lowest few modes** (every Monte Carlo sample) — subspace iteration on
 *   the banded factors, O(n p²) per sweep, converging at the ratio of the wanted
 *   eigenvalues to the first unwanted one. Iterating with m ≈ 2·nev vectors makes
 *   that ratio (nev / 2nev)⁴ for a beam, so a few sweeps reach round-off.
 *
 * The dense routine is the classical tred2/tql2 pair (Bowdler, Martin, Reinsch &
 * Wilkinson, Handbook for Automatic Computation II; as in EISPACK and JAMA).
 */

import { SymBand, cholSolve, cholesky, forwardSolve } from "./band";

export interface Eigen {
  /** Ascending. */
  readonly values: Float64Array;
  /** vectors[i] belongs to values[i]; absent when only values were asked for. */
  readonly vectors?: Float64Array[];
}

/** Eigen-decomposition of a dense symmetric matrix (A is not modified). */
export function symEig(A: ArrayLike<number>[], wantVectors = true): Eigen {
  const n = A.length;
  const V = Array.from({ length: n }, (_, i) => Float64Array.from(A[i]));
  const d = new Float64Array(n), e = new Float64Array(n);
  tred2(V, d, e, wantVectors);
  tql2(V, d, e, wantVectors);
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => d[a] - d[b]);
  const values = Float64Array.from(order, (i) => d[i]);
  if (!wantVectors) return { values };
  const vectors = order.map((j) => Float64Array.from({ length: n }, (_, i) => V[i][j]));
  return { values, vectors };
}

/** Householder reduction to tridiagonal form: diagonal into d, subdiagonal into e[1..]. */
function tred2(V: Float64Array[], d: Float64Array, e: Float64Array, accumulate: boolean): void {
  const n = V.length;
  for (let j = 0; j < n; j++) d[j] = V[n - 1][j];
  for (let i = n - 1; i > 0; i--) {
    let scale = 0, h = 0;
    for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
    if (scale === 0) {
      e[i] = d[i - 1];
      for (let j = 0; j < i; j++) { d[j] = V[i - 1][j]; V[i][j] = 0; V[j][i] = 0; }
    } else {
      for (let k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
      let f = d[i - 1], g = Math.sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (let j = 0; j < i; j++) e[j] = 0;
      for (let j = 0; j < i; j++) {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (let k = j + 1; k <= i - 1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
        e[j] = g;
      }
      f = 0;
      for (let j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
      const hh = f / (h + h);
      for (let j = 0; j < i; j++) e[j] -= hh * d[j];
      for (let j = 0; j < i; j++) {
        f = d[j]; g = e[j];
        for (let k = j; k <= i - 1; k++) V[k][j] -= f * e[k] + g * d[k];
        d[j] = V[i - 1][j];
        V[i][j] = 0;
      }
    }
    d[i] = h;
  }
  if (!accumulate) {
    // The tridiagonal's diagonal is left on V's, untouched by the reduction.
    for (let i = 0; i < n; i++) d[i] = V[i][i];
    e[0] = 0;
    return;
  }
  for (let i = 0; i < n - 1; i++) {
    V[n - 1][i] = V[i][i];
    V[i][i] = 1;
    const h = d[i + 1];
    if (h !== 0) {
      for (let k = 0; k <= i; k++) d[k] = V[k][i + 1] / h;
      for (let j = 0; j <= i; j++) {
        let g = 0;
        for (let k = 0; k <= i; k++) g += V[k][i + 1] * V[k][j];
        for (let k = 0; k <= i; k++) V[k][j] -= g * d[k];
      }
    }
    for (let k = 0; k <= i; k++) V[k][i + 1] = 0;
  }
  for (let j = 0; j < n; j++) { d[j] = V[n - 1][j]; V[n - 1][j] = 0; }
  V[n - 1][n - 1] = 1;
  e[0] = 0;
}

/** Implicit QL on the tridiagonal (d, e); rotations accumulated into V if asked. */
function tql2(V: Float64Array[], d: Float64Array, e: Float64Array, accumulate: boolean): void {
  const n = d.length;
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  let f = 0, tst1 = 0;
  const eps = 2 ** -52;
  for (let l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    let m = l;
    while (m < n - 1 && Math.abs(e[m]) > eps * tst1) m++;
    if (m > l) {
      let iter = 0;
      do {
        if (++iter > 60) throw new Error("tql2: no convergence");
        let g = d[l];
        let p = (d[l + 1] - g) / (2 * e[l]);
        let r = Math.hypot(p, 1);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const dl1 = d[l + 1];
        let h = g - d[l];
        for (let i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        p = d[m];
        let c = 1, c2 = c, c3 = c, s = 0, s2 = 0;
        const el1 = e[l + 1];
        for (let i = m - 1; i >= l; i--) {
          c3 = c2; c2 = c; s2 = s;
          g = c * e[i];
          h = c * p;
          r = Math.hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          if (accumulate)
            for (let k = 0; k < n; k++) {
              h = V[k][i + 1];
              V[k][i + 1] = s * V[k][i] + c * h;
              V[k][i] = c * V[k][i] - s * h;
            }
        }
        p = (-s * s2 * c3 * el1 * e[l]) / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (Math.abs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }
}

/** Dense Cholesky, lower factor; throws if A is not positive definite. */
function denseCholesky(A: ArrayLike<number>[]): Float64Array[] {
  const n = A.length, L = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++)
    for (let j = 0; j <= i; j++) {
      let s = A[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) {
        if (!(s > 0)) throw new Error(`matrix is not positive definite (pivot ${i})`);
        L[i][i] = Math.sqrt(s);
      } else L[i][j] = s / L[j][j];
    }
  return L;
}

/** K x = λ M x for small dense K, M (M SPD); vectors come out M-orthonormal. */
export function denseGeneralizedEig(K: ArrayLike<number>[], M: ArrayLike<number>[]): Required<Eigen> {
  const n = K.length, L = denseCholesky(M);
  // C = L⁻¹ K L⁻ᵀ, by two triangular solves.
  const W = Array.from({ length: n }, () => new Float64Array(n)); // W = L⁻¹ K, column by column
  for (let c = 0; c < n; c++)
    for (let i = 0; i < n; i++) {
      let s = K[i][c];
      for (let k = 0; k < i; k++) s -= L[i][k] * W[k][c];
      W[i][c] = s / L[i][i];
    }
  const C = Array.from({ length: n }, () => new Float64Array(n)); // C = L⁻¹ Wᵀ
  for (let c = 0; c < n; c++)
    for (let i = 0; i < n; i++) {
      let s = W[c][i];
      for (let k = 0; k < i; k++) s -= L[i][k] * C[k][c];
      C[i][c] = s / L[i][i];
    }
  symmetrise(C);
  const { values, vectors } = symEig(C, true) as Required<Eigen>;
  // x = L⁻ᵀ y.
  const xs = vectors.map((y) => {
    const x = Float64Array.from(y);
    for (let i = n - 1; i >= 0; i--) {
      let s = x[i];
      for (let k = i + 1; k < n; k++) s -= L[k][i] * x[k];
      x[i] = s / L[i][i];
    }
    return x;
  });
  return { values, vectors: xs };
}

function symmetrise(C: Float64Array[]): void {
  for (let i = 0; i < C.length; i++)
    for (let j = 0; j < i; j++) C[i][j] = C[j][i] = 0.5 * (C[i][j] + C[j][i]);
}

/** Every eigenvalue of K x = λ M x for banded K, M, via the dense standard form. */
export function bandSpectrum(K: SymBand, M: SymBand): Float64Array {
  const n = K.n, L = cholesky(M), Kd = K.toDense();
  // W = L⁻¹ K (rows of Kᵀ = K as columns), then C = L⁻¹ Wᵀ.
  const Wt = Kd.map((col) => forwardSolve(L, col)); // Wt[c] = column c of W
  const C = Array.from({ length: n }, () => new Float64Array(n));
  for (let i = 0; i < n; i++) {
    const row = new Float64Array(n);
    for (let c = 0; c < n; c++) row[c] = Wt[c][i]; // row i of W = column i of Wᵀ
    const col = forwardSolve(L, row);
    for (let r = 0; r < n; r++) C[r][i] = col[r];
  }
  symmetrise(C);
  return symEig(C, false).values;
}

/** Deterministic start vectors, so a solve is reproducible bit for bit. */
function startVectors(n: number, m: number): Float64Array[] {
  let s = 0x9e3779b9;
  const next = () => {
    s ^= s << 13; s ^= s >>> 17; s ^= s << 5; s >>>= 0;
    return s / 0x100000000 - 0.5;
  };
  return Array.from({ length: m }, () => Float64Array.from({ length: n }, next));
}

const dot = (a: ArrayLike<number>, b: ArrayLike<number>): number => {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
};

export interface Modes extends Required<Eigen> {
  readonly iterations: number;
}

/**
 * The `nev` lowest eigenpairs of K x = λ M x by subspace iteration, K and M
 * symmetric positive definite. Vectors are M-orthonormal.
 *
 * Two stopping tests, because the two converge at different rates: a Ritz
 * value's error is the square of its vector's, so a run that stops when the
 * values settle (`tol`) hands back vectors good to only √tol. The vectors are
 * therefore also required to settle (`vtol`) — or to stop improving, since the
 * floor they reach is round-off amplified by K's conditioning (~ h⁻⁴), which no
 * fixed tolerance can anticipate across meshes.
 *
 * The values have a floor of their own, and it is not ε. The Rayleigh quotient
 * xᵀKx sums entries of K far larger than λ₁ that cancel to it, so it carries an
 * error near ε·|x|ᵀ|K||x| — relative to λ₁, ε times the spread of the spectrum,
 * which on a fine cantilever is 10⁻¹⁰ and more. That bound is computed per mode
 * and `tol` is never asked to beat it.
 */
export function lowestModes(
  K: SymBand, M: SymBand, nev: number, tol = 1e-13, vtol = 1e-11, maxIt = 200,
): Modes {
  const n = K.n;
  nev = Math.min(nev, n);
  const m = Math.min(n, Math.max(2 * nev, nev + 8));
  if (m === n) {
    const full = denseGeneralizedEig(K.toDense(), M.toDense());
    return { values: full.values.slice(0, nev), vectors: full.vectors.slice(0, nev), iterations: 0 };
  }
  const LK = cholesky(K);
  // |K|, formed once: the noise bound reads it every sweep, and a plate's band is megabytes.
  const absK = new SymBand(n, K.bw, K.data.map(Math.abs));
  let X = startVectors(n, m);
  let prev = new Float64Array(nev).fill(Infinity);
  let prevDelta = Infinity;
  for (let it = 1; it <= maxIt; it++) {
    const MX = X.map((x) => M.matvec(x));
    const Y = MX.map((mx) => cholSolve(LK, mx));
    // Unit M-norm columns: K⁻¹ shrinks each direction by its eigenvalue, and the
    // spread would otherwise land in the conditioning of the projected Mr.
    const scale = new Float64Array(m);
    const MY = Y.map((y, j) => {
      const my = M.matvec(y), s = (scale[j] = 1 / Math.sqrt(dot(y, my)));
      y.forEach((v, i) => (y[i] = v * s));
      return my.map((v) => v * s);
    });
    // Yᵀ K Y without K: K Y = M X by construction, and M — positive entries,
    // spectrum bounded — forms its products without the cancellation K's would
    // suffer (Bathe's form of the iteration).
    const Kr = Y.map((yi) => Float64Array.from(MX, (mxj, j) => scale[j] * dot(yi, mxj)));
    const Mr = Y.map((yi) => Float64Array.from(MY, (myj) => dot(yi, myj)));
    symmetrise(Kr);
    symmetrise(Mr);
    const { values, vectors } = denseGeneralizedEig(Kr, Mr);
    const old = X;
    X = vectors.map((q) => {
      const x = new Float64Array(n);
      q.forEach((qj, j) => { for (let i = 0; i < n; i++) x[i] += qj * Y[j][i]; });
      return x;
    });
    let settled = true, delta = 0;
    for (let i = 0; i < nev; i++) {
      const x0 = X[i], noise = (16 * 2 ** -52 * dot(x0, absK.matvec(x0.map(Math.abs)))) / dot(x0, K.matvec(x0));
      if (Math.abs(values[i] - prev[i]) > Math.max(tol, noise) * Math.abs(values[i])) settled = false;
      // Change in direction, blind to the sign a Ritz vector comes back with —
      // as a difference of unit vectors, not √(1 − cos²), which cancels to √ε.
      const x = X[i], o = old[i], nx = Math.sqrt(dot(x, x));
      const so = Math.sign(dot(x, o)) / Math.sqrt(dot(o, o));
      let d2 = 0;
      for (let j = 0; j < n; j++) d2 += (x[j] / nx - so * o[j]) ** 2;
      delta = Math.max(delta, Math.sqrt(d2));
    }
    prev = values.slice(0, nev);
    if (settled && (delta < vtol || delta > 0.5 * prevDelta))
      return { values: prev, vectors: X.slice(0, nev), iterations: it };
    prevDelta = delta;
  }
  throw new Error(`subspace iteration did not converge in ${maxIt} sweeps`);
}
