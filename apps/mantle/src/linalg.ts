/** Small dense linear algebra. CPU/init only — f64, never in the GPU hot loop. */

export type LU = { A: Float64Array[]; piv: number[] };

export const mat = (m: number, n: number) =>
  Array.from({ length: m }, () => new Float64Array(n));

export function lu(A: Float64Array[]): LU {
  const n = A.length, piv = Array.from({ length: n }, (_, i) => i);
  for (let k = 0; k < n; k++) {
    let p = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i;
    if (p !== k) { [A[k], A[p]] = [A[p], A[k]]; [piv[k], piv[p]] = [piv[p], piv[k]]; }
    for (let i = k + 1; i < n; i++) {
      A[i][k] /= A[k][k];
      for (let j = k + 1; j < n; j++) A[i][j] -= A[i][k] * A[k][j];
    }
  }
  return { A, piv };
}

/** Thomas-algorithm factorisation of a tridiagonal system, reusable across solves. */
export type Tri = { a: Float64Array; cp: Float64Array; m: Float64Array };

export function triFactor(a: Float64Array, b: Float64Array, c: Float64Array): Tri {
  const n = b.length, cp = new Float64Array(n), m = new Float64Array(n);
  m[0] = b[0];
  cp[0] = c[0] / m[0];
  for (let i = 1; i < n; i++) {
    m[i] = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / m[i];
  }
  return { a, cp, m };
}

export function triSolve({ a, cp, m }: Tri, d: Float64Array): Float64Array {
  const n = m.length, x = new Float64Array(n);
  x[0] = d[0] / m[0];
  for (let i = 1; i < n; i++) x[i] = (d[i] - a[i] * x[i - 1]) / m[i];
  for (let i = n - 2; i >= 0; i--) x[i] -= cp[i] * x[i + 1];
  return x;
}

export function solve({ A, piv }: LU, b: Float64Array): Float64Array {
  const n = A.length, x = new Float64Array(n);
  for (let i = 0; i < n; i++) { let s = b[piv[i]]; for (let j = 0; j < i; j++) s -= A[i][j] * x[j]; x[i] = s; }
  for (let i = n - 1; i >= 0; i--) { let s = x[i]; for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j]; x[i] = s / A[i][i]; }
  return x;
}
