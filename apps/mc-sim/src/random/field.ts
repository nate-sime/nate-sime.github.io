// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The random inputs of one sample ω, at a fixed set of points.
 *
 * Stiffness is lognormal about the deterministic section,
 *
 *   e(x, ω) = e₀(x) · exp(σ g(x, ω) − ½ σ² s_M(x)),   s_M(x) = Σ_{j<M} λ_j φ_j(x)²,
 *
 * g the truncated KL field of `kl.ts`. Subtracting the truncated pointwise
 * variance s_M rather than the full σ² makes E[e(x)] = e₀(x) exactly, at every
 * x and for every M: refining the truncation changes the fluctuations, never
 * the mean. The coefficient of variation is √(e^{σ² s_M} − 1) — about σ for
 * small σ.
 *
 * Mass, optionally, follows the stiffness as a rectangular section of random
 * depth d would: I ∝ d³ and A ∝ d, so μ = μ₀ (e/e₀)^{1/3}. Frequencies, which go
 * as √(e/μ), then feel the field as (e/e₀)^{1/3} rather than (e/e₀)^{1/2}.
 *
 * Load, optionally, gets an independent Gaussian part with the same kernel:
 * q(x, ω) = 1 + σ_q g′(x, ω), g′ from the second normal channel. A point load
 * gets a random magnitude 1 + σ_q ξ′₀ instead. q is not lognormal — a
 * fluctuating pressure may change sign.
 *
 * `FieldAt` tabulates the modes at one set of points — a level's quadrature
 * points, or a plot's — once; a sample there is then a matrix–vector product.
 * Because every level tabulates the same expansion and every sample reads the
 * same ξ (`normals` is addressed by sample index), sample ω is the same function
 * on the fine and the coarse mesh: the coupling MLMC rests on.
 */

import { CHANNEL, normals } from "./philox";
import { modesAt, termsOf, type KLData, type Kernel } from "./kl";

export interface FieldSpec {
  readonly kernel: Kernel;
  /** Correlation length, in units of the span. */
  readonly ell: number;
  /** Standard deviation of log e (before truncation). */
  readonly sigma: number;
  /** KL terms asked for; the expansion may keep fewer (see `kl.ts`). */
  readonly terms: number;
  /** μ ∝ e^{1/3}: a random section depth rather than a random modulus. */
  readonly massFollows: boolean;
  /** Standard deviation of the load's random part; 0 for a deterministic load. */
  readonly loadSigma: number;
}

export class FieldAt {
  /** Number of points. */
  readonly n: number;
  readonly M: number;
  /** √λ_j φ_j(x_q), flat [q·M + j]. */
  readonly phi: Float64Array;
  /** s_M(x_q) = Σ_j λ_j φ_j(x_q)². */
  readonly s: Float64Array;

  constructor(kl: KLData, readonly points: ArrayLike<number>) {
    this.n = points.length;
    this.M = termsOf(kl);
    this.phi = modesAt(kl, points);
    this.s = new Float64Array(this.n);
    for (let q = 0; q < this.n; q++) {
      let v = 0;
      for (let j = 0; j < this.M; j++) v += this.phi[q * this.M + j] ** 2;
      this.s[q] = v;
    }
  }

  /** g(x_q) = Σ_j √λ_j φ_j(x_q) ξ_j; ξ may be longer than M (extra terms ignored). */
  gaussian(xi: ArrayLike<number>, out = new Float64Array(this.n)): Float64Array {
    const { M, phi } = this;
    for (let q = 0; q < this.n; q++) {
      let v = 0;
      for (let j = 0; j < M; j++) v += phi[q * M + j] * xi[j];
      out[q] = v;
    }
    return out;
  }

  /** The lognormal factor e/e₀ = exp(σ g − ½σ² s_M). */
  lognormal(sigma: number, xi: ArrayLike<number>, out = new Float64Array(this.n)): Float64Array {
    this.gaussian(xi, out);
    for (let q = 0; q < this.n; q++) out[q] = Math.exp(sigma * out[q] - 0.5 * sigma * sigma * this.s[q]);
    return out;
  }
}

/** The normals behind sample `index`: ξ for the stiffness and, if the load is random, ξ′ for it. */
export interface Draw {
  readonly xi: Float64Array;
  readonly eta: Float64Array | null;
}

export function draw(spec: FieldSpec, M: number, seed: number, index: number): Draw {
  return {
    xi: normals(seed, index, CHANNEL.stiffness, M),
    eta: spec.loadSigma > 0 ? normals(seed, index, CHANNEL.load, Math.max(M, 1)) : null,
  };
}

/** One sample's coefficients at a FieldAt's points, about the deterministic e₀, μ₀. */
export interface Realised {
  /** e(x_q). */
  readonly e: Float64Array;
  /** μ(x_q). */
  readonly mu: Float64Array;
  /** q(x_q), the distributed load per unit reference load; null if deterministic. */
  readonly q: Float64Array | null;
  /** Point-load magnitude per unit reference load. */
  readonly P: number;
}

export function realise(
  f: FieldAt, spec: FieldSpec, d: Draw, e0: Float64Array, mu0: Float64Array,
): Realised {
  const factor = f.lognormal(spec.sigma, d.xi);
  const e = new Float64Array(f.n), mu = new Float64Array(f.n);
  for (let q = 0; q < f.n; q++) {
    e[q] = e0[q] * factor[q];
    mu[q] = spec.massFollows ? mu0[q] * Math.cbrt(factor[q]) : mu0[q];
  }
  let q: Float64Array | null = null, P = 1;
  if (d.eta) {
    q = f.gaussian(d.eta);
    for (let i = 0; i < f.n; i++) q[i] = 1 + spec.loadSigma * q[i];
    P = 1 + spec.loadSigma * d.eta[0];
  }
  return { e, mu, q, P };
}

/** Pointwise quantiles of e/e₀: exp(z σ √s_M − ½σ² s_M) — lognormal, so exact. */
export function lognormalQuantile(f: FieldAt, sigma: number, z: number): Float64Array {
  return f.s.map((s) => Math.exp(z * sigma * Math.sqrt(s) - 0.5 * sigma * sigma * s));
}
