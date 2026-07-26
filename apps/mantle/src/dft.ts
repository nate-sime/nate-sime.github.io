/**
 * Azimuthal DFT along the second index. O(n²) — adequate for the CPU reference;
 * the GPU path replaces it with a Stockham autosort FFT.
 *
 * Convention: X̂[i][k] = (1/n) Σ_j X[i][j] e^{−2πijk/n}, so a field constant in φ
 * transforms to (value, 0, 0, …). This is the convention under which circulant
 * operators are diagonalised.
 */

import { mat } from "./linalg";

export function forward(X: Float64Array[]): { re: Float64Array[]; im: Float64Array[] } {
  const m = X.length, n = X[0].length, re = mat(m, n), im = mat(m, n);
  for (let i = 0; i < m; i++)
    for (let k = 0; k < n; k++) {
      let sr = 0, si = 0;
      for (let j = 0; j < n; j++) {
        const t = (-2 * Math.PI * j * k) / n;
        sr += X[i][j] * Math.cos(t);
        si += X[i][j] * Math.sin(t);
      }
      re[i][k] = sr / n;
      im[i][k] = si / n;
    }
  return { re, im };
}

export function inverse(re: Float64Array[], im: Float64Array[]): Float64Array[] {
  const m = re.length, n = re[0].length, X = mat(m, n);
  for (let i = 0; i < m; i++)
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) {
        const t = (2 * Math.PI * j * k) / n;
        s += re[i][k] * Math.cos(t) - im[i][k] * Math.sin(t);
      }
      X[i][j] = s;
    }
  return X;
}
