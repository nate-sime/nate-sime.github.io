// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * What is known about the uniform plate (D = μ = 1) without discretising it.
 *
 * Simply supported on all four edges, everything is a Navier series. The
 * modes are φ_mn = (2/√a) sin(mπx/a) sin(nπy), M-orthonormal on a × 1, with
 *
 *   λ_mn = ω̂_mn² = π⁴ ((m/a)² + n²)²,
 *
 * and a load with sine coefficients q_mn deflects the plate by
 * w_mn = q_mn / λ_mn in each — independent of ν, because the twisting term of
 * the energy integrates to a boundary term that vanishes on a simply supported
 * edge. A uniform load has q_mn = 16/(π² m n) for odd m, n; a unit point load
 * at (ξ, η) has q_mn = (4/a) sin(mπξ/a) sin(nπη).
 *
 * The uniform-load series falls as (mn)⁻¹(m² + n²)⁻² and converges fast. The
 * point-load series at the load point falls only as (m² + n²)⁻², and is
 * summed to m, n ≤ 1999: its tail is then ~10⁻⁷ of the answer, below any error
 * a mesh the app runs can reach — the point-load solution goes as r² log r
 * under the load, so its spline error falls only as h².
 *
 * Clamped on two opposite edges and simply supported on the others (CSCS), the
 * frequencies are Lévy's: exact, from a transcendental equation per n. Every
 * other support has no closed form; for those the app quotes reference
 * frequencies for the square plate (`REFERENCE`), which the tests check the
 * splines converge to.
 */

import type { Harmonic } from "../hierarchy";
import type { EdgeName } from "./plate";

export function navierEigenvalues(aspect: number, count: number): Float64Array {
  const mn: number[] = [];
  const top = Math.ceil(Math.sqrt(count)) + 2;
  for (let m = 1; m <= top * Math.ceil(aspect + 1); m++)
    for (let n = 1; n <= top * Math.ceil(1 / aspect + 1); n++) mn.push(Math.PI ** 4 * ((m / aspect) ** 2 + n * n) ** 2);
  return Float64Array.from(mn.sort((x, y) => x - y).slice(0, count));
}

/** (m, n) of the modes in ascending order, ties by m. */
export function navierModes(aspect: number, count: number): [number, number][] {
  const out: [number, number, number][] = [];
  const top = Math.ceil(Math.sqrt(count)) + 2;
  for (let m = 1; m <= top * Math.ceil(aspect + 1); m++)
    for (let n = 1; n <= top * Math.ceil(1 / aspect + 1); n++) out.push([m, n, (m / aspect) ** 2 + n * n]);
  return out.sort((x, y) => x[2] - y[2] || x[0] - y[0]).slice(0, count).map(([m, n]) => [m, n]);
}

export interface PlateStatic {
  /** w at one point. */
  readonly w: (x: number, y: number) => number;
  /** w on the grid xs × ys, flat [iy·xs.length + ix]; uniform load only (null for a point load). */
  readonly grid: ((xs: ArrayLike<number>, ys: ArrayLike<number>) => Float64Array) | null;
  /** ℓ(w), the work done by the load. */
  readonly compliance: number;
  /** ‖w‖ = (∫ w²)^½. */
  readonly rms: number;
}

/**
 * Terms per direction: m, n ≤ 1001 for the uniform load (a tail of ~10⁻¹² of
 * w anywhere), m, n ≤ 1999 for a point load.
 */
const TERMS = { uniform: 1001, point: 1999 } as const;

const sines = (N: number, t: number): Float64Array => Float64Array.from({ length: N + 1 }, (_, m) => Math.sin(m * Math.PI * t));

/** The simply supported plate under a unit uniform load, or a unit point load at (ξ, η). */
export function navierStatic(aspect: number, load: "uniform" | { x: number; y: number }): PlateStatic {
  const a = aspect, uniform = load === "uniform", N = uniform ? TERMS.uniform : TERMS.point;
  const lam = (m: number, n: number) => Math.PI ** 4 * ((m / a) ** 2 + n * n) ** 2;
  // Sine coefficients of the load, separable: q_mn = qx_m qy_n.
  const qx = new Float64Array(N + 1), qy = new Float64Array(N + 1);
  for (let m = 1; m <= N; m++) {
    qx[m] = uniform ? (m % 2 ? 4 / (Math.PI * m) : 0) : (4 / a) * Math.sin((m * Math.PI * load.x) / a);
    qy[m] = uniform ? (m % 2 ? 4 / (Math.PI * m) : 0) : Math.sin(m * Math.PI * load.y);
  }
  // ℓ(w) = Σ q_mn w_mn ∫sin∫sin for the uniform load (∫₀ᵃ sin(mπx/a) = 2a/(mπ) for odd m), = w(ξ, η) for a point load.
  let compliance = 0, sum2 = 0;
  for (let m = 1; m <= N; m++) {
    if (qx[m] === 0) continue;
    for (let n = 1; n <= N; n++) {
      if (qy[n] === 0) continue;
      const w = (qx[m] * qy[n]) / lam(m, n);
      sum2 += w * w;
      compliance += uniform ? w * ((2 * a) / (m * Math.PI)) * (2 / (n * Math.PI)) : w * (qx[m] * a / 4) * qy[n];
    }
  }
  const at = (sx: Float64Array, sy: Float64Array) => {
    let s = 0;
    for (let m = 1; m <= N; m++) {
      if (qx[m] === 0 || sx[m] === 0) continue;
      let t = 0;
      for (let n = 1; n <= N; n++) t += (qy[n] * sy[n]) / lam(m, n);
      s += qx[m] * sx[m] * t;
    }
    return s;
  };
  return {
    w: (x, y) => at(sines(N, x / a), sines(N, y)),
    grid: uniform ? (xs, ys) => {
      // w = Sx · (q_mn/λ_mn) · Syᵀ, odd terms only.
      const odd = Array.from({ length: (N + 1) / 2 }, (_, i) => 2 * i + 1), K = odd.length;
      const T = new Float64Array(ys.length * K); // T[iy·K + im] = Σ_n w_mn sin(nπy)
      for (let iy = 0; iy < ys.length; iy++) {
        const sy = sines(N, ys[iy]);
        for (let im = 0; im < K; im++) {
          const m = odd[im];
          let t = 0;
          for (const n of odd) t += (qy[n] * sy[n]) / lam(m, n);
          T[iy * K + im] = qx[m] * t;
        }
      }
      const out = new Float64Array(xs.length * ys.length);
      for (let ix = 0; ix < xs.length; ix++) {
        const sx = sines(N, xs[ix] / a);
        for (let iy = 0; iy < ys.length; iy++) {
          let s = 0;
          for (let im = 0; im < K; im++) s += T[iy * K + im] * sx[odd[im]];
          out[iy * xs.length + ix] = s;
        }
      }
      return out;
    } : null,
    compliance,
    rms: Math.sqrt((a / 4) * sum2),
  };
}

/**
 * Steady amplitude |u(x, y)| of the simply supported plate under the load
 * applied harmonically at Ω with Rayleigh damping C = aM + bK — modal
 * superposition over the exact modes, which Rayleigh damping leaves uncoupled:
 * u = Σ φ_mn (φ_mn, F) / (λ_mn − Ω² + iΩ(a + bλ_mn)).
 */
export function navierResponse(aspect: number, load: "uniform" | "point", at: { x: number; y: number }, h: Harmonic): number {
  const A = aspect, N = TERMS[load];
  let re = 0, im = 0;
  for (let m = 1; m <= N; m++) {
    const sx = Math.sin((m * Math.PI * at.x) / A);
    for (let n = 1; n <= N; n++) {
      // φ_mn at the point, and (φ_mn, F): ∫ φ for the uniform load, φ(at) for a point load there.
      const phi = (2 / Math.sqrt(A)) * sx * Math.sin(n * Math.PI * at.y);
      let f: number;
      if (load === "uniform") {
        if (m % 2 === 0 || n % 2 === 0) continue;
        f = (2 / Math.sqrt(A)) * ((2 * A) / (m * Math.PI)) * (2 / (n * Math.PI));
      } else f = phi;
      if (phi === 0 || f === 0) continue;
      const lam = Math.PI ** 4 * ((m / A) ** 2 + n * n) ** 2;
      const dr = lam - h.Omega ** 2, di = h.Omega * (h.a + h.b * lam), d2 = dr * dr + di * di;
      re += (phi * f * dr) / d2;
      im -= (phi * f * di) / d2;
    }
  }
  return Math.hypot(re, im);
}

/**
 * CSCS — clamped on x = 0 and x = a, simply supported on y = 0 and y = 1 — is
 * Lévy's case: w = X(x) sin(nπy) separates, and X solves
 * X⁗ − 2k²X″ + (k⁴ − λ)X = 0, k = nπ, with X = X′ = 0 at both ends. With
 * r₁ = (√λ + k²)^½, r₂ = (√λ − k²)^½ and ξ = x − a/2, h = a/2, the symmetric
 * modes A cosh r₁ξ + C cos r₂ξ and the antisymmetric B sinh r₁ξ + D sin r₂ξ
 * vanish with their slopes at ξ = h when
 *
 *   r₂ sin r₂h + r₁ tanh r₁h cos r₂h = 0,   r₂ tanh r₁h cos r₂h − r₁ sin r₂h = 0,
 *
 * both divided through by cosh r₁h so they stay O(1). Each is scanned in r₂
 * and its roots bisected; λ = (r₂² + k²)².
 */
export function levyEigenvalues(aspect: number, count: number): Float64Array {
  const h = aspect / 2, found: number[] = [];
  for (let n = 1; n <= count + 1; n++) {
    const k = n * Math.PI;
    const r1 = (r2: number) => Math.sqrt(r2 * r2 + 2 * k * k);
    const fs = (r2: number) => r2 * Math.sin(r2 * h) + r1(r2) * Math.tanh(r1(r2) * h) * Math.cos(r2 * h);
    const fa = (r2: number) => r2 * Math.tanh(r1(r2) * h) * Math.cos(r2 * h) - r1(r2) * Math.sin(r2 * h);
    for (const f of [fs, fa]) {
      // Roots are spaced about π/h apart in r₂; steps of a twentieth of that cannot skip a pair.
      const step = Math.PI / h / 20, top = (count + 2) * (Math.PI / h);
      let a = 1e-9, fa0 = f(a);
      for (let b = a + step; b <= top; b += step) {
        const fb = f(b);
        if (Math.sign(fb) !== Math.sign(fa0)) {
          let lo = a, hi = b, flo = fa0;
          for (let it = 0; it < 200 && hi - lo > 1e-15 * hi; it++) {
            const mid = 0.5 * (lo + hi), fm = f(mid);
            if (Math.sign(fm) === Math.sign(flo)) { lo = mid; flo = fm; } else hi = mid;
          }
          const r2 = 0.5 * (lo + hi);
          found.push((r2 * r2 + k * k) ** 2);
        }
        a = b;
        fa0 = fb;
      }
    }
  }
  return Float64Array.from(found.sort((x, y) => x - y).slice(0, count));
}

/** Exact eigenvalues where the edges allow them: Navier (SSSS) and Lévy (CSCS), any aspect ratio. */
export function exactPlateEigenvalues(edges: string, aspect: number, count: number): Float64Array | null {
  if (edges === "SSSS") return navierEigenvalues(aspect, count);
  if (edges === "CSCS") return levyEigenvalues(aspect, count);
  return null;
}

/**
 * Reference frequencies ω̂ = ω b² √(ρt/D) for the square plate without a closed
 * form, ν = 0.3 where it matters (not for CCCC: a clamped edge kills the
 * twisting term). CCCC and FFFF are Leissa's (Vibration of Plates, NASA SP-160,
 * 1969); the splines reproduce every digit he gives. CFFF is not: his table has
 * 3.4917, 8.5246, 21.429, and the splines converge from above to the values
 * below, 0.6% under his first — a conforming Rayleigh–Ritz method only
 * overestimates, so his are not the limit. These are p = 4 on 48 × 48 and
 * 64 × 64 elements, agreeing to the digits shown.
 */
export const REFERENCE: Partial<Record<EdgeName, { readonly omega: readonly number[]; readonly source: string }>> = {
  CCCC: { omega: [35.985, 73.394, 73.394, 108.22], source: "Leissa (1969)" },
  // The three rigid modes come first; these are the 4th, 5th and 6th.
  FFFF: { omega: [0, 0, 0, 13.468, 19.596, 24.270], source: "Leissa (1969)" },
  CFFF: { omega: [3.4710, 8.506, 21.284], source: "splines, converged" },
};
/** Timoshenko & Woinowsky-Krieger: the clamped square plate's centre deflection under a uniform load, w D / (q b⁴). */
export const CLAMPED_CENTRE = 0.00126;
