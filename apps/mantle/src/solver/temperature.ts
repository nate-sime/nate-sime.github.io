/**
 * Temperature transport on a uniform (r, φ) grid.
 *
 * `T` is carried as point values, not spline coefficients: semi-Lagrangian
 * advection, BFECC and monotone limiting are all point-value operations, and
 * coefficients would force a global interpolation solve every step.
 *
 *   advection — semi-Lagrangian, backward RK2 characteristic trace with the
 *               exactly divergence-free velocity, monotone-clamped bicubic
 *               interpolation, BFECC error correction. Unconditionally stable in
 *               the advective CFL, so dt is an accuracy knob.
 *   diffusion — implicit (backward Euler). Spectral in φ (the grid is uniform and
 *               periodic, so the DFT diagonalises ∂_φφ exactly as −k²) and
 *               2nd-order FD in r, giving one tridiagonal solve per mode. This
 *               removes the dt ~ h² limit and reuses the mode-decoupling of §6.1.
 *
 * Isothermal boundaries are Dirichlet and constant in time, so in mode space they
 * contribute only to k = 0; all other modes see homogeneous data.
 */

import { mat, triFactor, triSolve, type Tri } from "../linalg";
import * as dft from "../dft";

export type Velocity = (r: number, phi: number) => { ur: number; up: number };

/**
 * Thomas factors of (I − dt ∇²) for azimuthal modes k = 0 … na/2.
 *
 * ∂_φφ contributes −k² and the radial part is dt-independent, so the factors
 * depend on k only through k² — identical for k and na−k. Storing the half range
 * and indexing by `min(k, na−k)` mirrors the Stokes blocks (`radialOperator`).
 * Built in f64 once per dt; the GPU path fixes dt so this happens at init only.
 */
export function diffusionFactors(
  nr: number, na: number, ri: number, dr: number, dt: number,
): Tri[] {
  const ni = nr - 2, out: Tri[] = [];
  for (let k = 0; k <= na / 2; k++) {
    const a = new Float64Array(ni), b = new Float64Array(ni), c = new Float64Array(ni);
    for (let i = 0; i < ni; i++) {
      const r = ri + (i + 1) * dr;
      a[i] = -dt * (1 / (dr * dr) - 1 / (2 * r * dr));
      b[i] = 1 - dt * (-2 / (dr * dr) - (k * k) / (r * r));
      c[i] = -dt * (1 / (dr * dr) + 1 / (2 * r * dr));
    }
    out.push(triFactor(a, b, c));
  }
  return out;
}

/** Catmull–Rom, clamped to the bracketing values — monotone, no over/undershoot. */
function cubic(p0: number, p1: number, p2: number, p3: number, t: number): number {
  const v = p1 + 0.5 * t * (p2 - p0 + t * (2 * p0 - 5 * p1 + 4 * p2 - p3
    + t * (3 * (p1 - p2) + p3 - p0)));
  return Math.min(Math.max(v, Math.min(p1, p2)), Math.max(p1, p2));
}

export class Temperature {
  T: Float64Array[];
  readonly dr: number;
  readonly dphi: number;
  private tri: { dt: number; f: Tri[] } | null = null;

  constructor(
    readonly nr: number, readonly na: number,
    readonly ri: number, readonly ro: number,
    readonly tIn = 1, readonly tOut = 0,
  ) {
    this.dr = (ro - ri) / (nr - 1);
    this.dphi = (2 * Math.PI) / na;
    this.T = mat(nr, na);
    this.reset();
  }

  r(i: number): number { return this.ri + i * this.dr; }

  /** Steady conduction profile ln(r_o/r)/ln(r_o/r_i) — the Nu = 1 reference. */
  conduction(i: number): number {
    return this.tOut + (this.tIn - this.tOut)
      * (Math.log(this.ro / this.r(i)) / Math.log(this.ro / this.ri));
  }

  /** Conductive profile plus an azimuthal perturbation to seed convection. */
  reset(amp = 0.05, mode = 4): void {
    for (let i = 0; i < this.nr; i++) {
      const s = i / (this.nr - 1);
      for (let j = 0; j < this.na; j++)
        this.T[i][j] = this.conduction(i)
          + amp * Math.sin(Math.PI * s) * Math.cos(mode * j * this.dphi);
    }
    this.applyBC();
  }

  applyBC(): void {
    this.T[0].fill(this.tIn);
    this.T[this.nr - 1].fill(this.tOut);
  }

  /** Monotone bicubic sample; r clamped to the domain, φ wrapped. */
  sample(src: Float64Array[], r: number, phi: number): number {
    const { nr, na, dr, dphi } = this;
    const x = Math.min(nr - 1, Math.max(0, (r - this.ri) / dr));
    const i = Math.min(nr - 2, Math.floor(x)), tr = x - i;
    let y = phi / dphi;
    y -= Math.floor(y / na) * na;
    const j = Math.floor(y), tp = y - j;

    const col = new Float64Array(4);
    for (let m = -1; m <= 2; m++) {
      const row = src[Math.min(nr - 1, Math.max(0, i + m))];
      col[m + 1] = cubic(row[(j - 1 + na) % na], row[j % na],
        row[(j + 1) % na], row[(j + 2) % na], tp);
    }
    return cubic(col[0], col[1], col[2], col[3], tr);
  }

  /** Local extrema of the cell containing (r, φ) — the BFECC limiter bracket. */
  private bracket(src: Float64Array[], r: number, phi: number): [number, number] {
    const { nr, na, dr, dphi } = this;
    const x = Math.min(nr - 1, Math.max(0, (r - this.ri) / dr));
    const i = Math.min(nr - 2, Math.floor(x));
    let y = phi / dphi;
    y -= Math.floor(y / na) * na;
    const j = Math.floor(y);
    let lo = Infinity, hi = -Infinity;
    for (let m = 0; m <= 1; m++)
      for (let l = 0; l <= 1; l++) {
        const v = src[i + m][(j + l) % na];
        lo = Math.min(lo, v); hi = Math.max(hi, v);
      }
    return [lo, hi];
  }

  /** Backward RK2 trace from (r, φ); dt < 0 traces forward. */
  private departure(u: Velocity, r: number, phi: number, dt: number): [number, number] {
    const a = u(r, phi);
    const rm = r - 0.5 * dt * a.ur, pm = phi - 0.5 * dt * (a.up / r);
    const b = u(Math.min(this.ro, Math.max(this.ri, rm)), pm);
    return [r - dt * b.ur, phi - dt * (b.up / r)];
  }

  private advectOnce(src: Float64Array[], u: Velocity, dt: number): Float64Array[] {
    const out = mat(this.nr, this.na);
    for (let i = 1; i < this.nr - 1; i++) {
      const r = this.r(i);
      for (let j = 0; j < this.na; j++) {
        const [rd, pd] = this.departure(u, r, j * this.dphi, dt);
        out[i][j] = this.sample(src, rd, pd);
      }
    }
    return out;
  }

  /**
   * BFECC: a forward then reverse SL pass estimates the (dissipative) error,
   * which is subtracted before the real advection — recovering high order for
   * roughly triple the cost. `bfecc = false` gives plain semi-Lagrangian.
   */
  advect(u: Velocity, dt: number, bfecc = true): void {
    let src = this.T;
    if (bfecc) {
      const T1 = this.advectOnce(this.T, u, dt);
      T1[0].fill(this.tIn); T1[this.nr - 1].fill(this.tOut);
      const T2 = this.advectOnce(T1, u, -dt);
      const corrected = mat(this.nr, this.na);
      for (let i = 0; i < this.nr; i++)
        for (let j = 0; j < this.na; j++)
          corrected[i][j] = this.T[i][j] + 0.5 * (this.T[i][j] - T2[i][j]);
      corrected[0].fill(this.tIn); corrected[this.nr - 1].fill(this.tOut);
      src = corrected;
    }

    const out = this.advectOnce(src, u, dt);
    if (bfecc)
      for (let i = 1; i < this.nr - 1; i++) {
        const r = this.r(i);
        for (let j = 0; j < this.na; j++) {
          const [rd, pd] = this.departure(u, r, j * this.dphi, dt);
          const [lo, hi] = this.bracket(this.T, rd, pd);
          out[i][j] = Math.min(Math.max(out[i][j], lo), hi); // reject new extrema
        }
      }
    this.T = out;
    this.applyBC();
  }

  /** Implicit diffusion: (I − dt ∇²) Tⁿ⁺¹ = Tⁿ, one tridiagonal solve per mode. */
  diffuse(dt: number): void {
    const { nr, na, dr } = this, ni = nr - 2;
    if (!this.tri || this.tri.dt !== dt)
      this.tri = { dt, f: diffusionFactors(nr, na, this.ri, dr, dt) };

    const { re, im } = dft.forward(this.T);
    const xr = mat(nr, na), xi = mat(nr, na);
    for (let k = 0; k < na; k++) {
      const br = new Float64Array(ni), bi = new Float64Array(ni);
      for (let i = 0; i < ni; i++) { br[i] = re[i + 1][k]; bi[i] = im[i + 1][k]; }
      // Dirichlet data is constant in φ, so it enters the k = 0 mode only.
      if (k === 0) {
        const r1 = this.r(1), rn = this.r(nr - 2);
        br[0] += dt * (1 / (dr * dr) - 1 / (2 * r1 * dr)) * this.tIn;
        br[ni - 1] += dt * (1 / (dr * dr) + 1 / (2 * rn * dr)) * this.tOut;
      }
      const f = this.tri.f[Math.min(k, na - k)];
      const sr = triSolve(f, br), si = triSolve(f, bi);
      for (let i = 0; i < ni; i++) { xr[i + 1][k] = sr[i]; xi[i + 1][k] = si[i]; }
    }
    this.T = dft.inverse(xr, xi);
    this.applyBC();
  }

  /**
   * Nusselt number at each boundary: total heat flux over the conductive flux
   * 2π/ln(r_o/r_i). Agreement between the two is a global heat-balance check.
   */
  nusselt(): { inner: number; outer: number } {
    const { nr, na, dr, dphi } = this;
    const norm = Math.log(this.ro / this.ri) / (2 * Math.PI);
    // 4th-order one-sided ∂_r. Nu is the benchmark quantity, so the boundary
    // flux must not be the accuracy bottleneck: a 2nd-order stencil leaves an
    // O(dr²) bias of ~1e-4 that would swamp the comparison.
    const d4 = (f: (m: number) => number) =>
      (-25 * f(0) + 48 * f(1) - 36 * f(2) + 16 * f(3) - 3 * f(4)) / (12 * dr);
    let qi = 0, qo = 0;
    for (let j = 0; j < na; j++) {
      qi += -d4((m) => this.T[m][j]);
      qo += d4((m) => this.T[nr - 1 - m][j]);
    }
    return { inner: norm * this.ri * qi * dphi, outer: norm * this.ro * qo * dphi };
  }

  /** Max |u| dt / cell — used to size the step for accuracy, not stability. */
  maxSpeed(u: Velocity): number {
    let m = 0;
    for (let i = 1; i < this.nr - 1; i++) {
      const r = this.r(i);
      for (let j = 0; j < this.na; j++) {
        const { ur, up } = u(r, j * this.dphi);
        m = Math.max(m, Math.abs(ur) / this.dr, Math.abs(up) / (r * this.dphi));
      }
    }
    return m;
  }
}
