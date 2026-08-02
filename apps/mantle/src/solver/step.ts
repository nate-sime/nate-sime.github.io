/**
 * Coupled Boussinesq convection loop. Per step:
 *
 *   1. buoyancy load from T          4. semi-Lagrangian + BFECC advection of T
 *   2. Stokes solve → ψ              5. implicit diffusion of T
 *   3. recover u = ∇×ψ
 *
 * Stokes is quasi-static — a constraint re-solved each step, not a time
 * evolution. With SL advection and implicit diffusion both unconditionally
 * stable, dt is limited by accuracy alone.
 */

import { ANNULUS, type Geometry } from "../geometry";
import { Axis, P, clampedAxis, periodicAxis, Field } from "../spline";
import { mat } from "../linalg";
import { gauss } from "../quad";
import { StokesSolver, VariableStokes } from "./stokes";
import { viscosityAt, strainRate } from "./assembly";
import {
  meanViscosity, viscosity, strainScale, depthAt,
  meanTackleyViscosity, tackleyViscosity,
} from "./rheology";
import { Temperature, type Velocity } from "./temperature";

/**
 * Buoyancy load ℓ(v) = ∫ Ra T (ĝ·u[v]) dx = Ra ∫∫ T B_m N_l' dr dφ.
 *
 * The 1/h of the along-gravity velocity `u_r[v] = v_φ/h` cancels the h of the
 * area element, so there is no metric weight at all — **this assembly is
 * identical in both geometries**, which is the one place the box needed nothing
 * doing to it.
 */
export function buoyancyLoad(
  rAx: Axis, aAx: Axis, T: (r: number, phi: number) => number, Ra: number,
): Float64Array[] {
  const b = mat(rAx.n, aAx.n);
  const aq = aAx.elements().flatMap(([a, c]) =>
    gauss(a, c).map(({ x, w }) => ({ x, w, ...aAx.ders(x, 1) })));

  for (const [ea, eb] of rAx.elements())
    for (const { x: r, w: wr } of gauss(ea, eb)) {
      const R = rAx.ders(r, 0);
      for (const A of aq) {
        const f = Ra * T(r, A.x) * wr * A.w;
        for (let p = 0; p <= P; p++) {
          const i = rAx.dof(R.span, p), c = f * R.N[0][p];
          for (let q = 0; q <= P; q++) b[i][aAx.dof(A.span, q)] += c * A.N[1][q];
        }
      }
    }
  return b;
}

export interface SimOptions {
  nr?: number; na?: number;      // spline space for ψ
  gnr?: number; gna?: number;    // grid for T
  /** Domain and metric; defaults to the annulus. See `src/geometry.ts`. */
  geom?: Geometry;
  /**
   * Initial perturbation of the conduction profile: amplitude, and a count of
   * cells around the domain. Exposed because the mode is what selects *which*
   * linear-stability wavenumber a run is testing — see the critical-Ra check in
   * `tests/temperature.test.ts` — and `Temperature`'s default is a mantle-like 4.
   */
  seed?: { amp?: number; mode?: number };
  Ra?: number; cfl?: number; dtMax?: number;
  /** Tier 2: μ(T, d) by matrix-free PCG instead of the direct DFT solve. */
  variable?: boolean;
  /** ln of the viscosity contrast across the temperature range. */
  gamma?: number;
  /**
   * ln of the viscosity contrast across the depth of the layer — c in
   * `μ ∝ exp(c(d − ½))`. Positive stiffens the deep interior, which is the sign
   * the mantle has and the sign the Blankenbach 2b case states.
   */
  cz?: number;
  /** Krylov budget per solve, in the variable-μ tier. */
  iters?: number;
  /** Power-law index. 1 is μ(T, d) exactly; ≈3 is dislocation creep. */
  n?: number;
  /** Rheology updates per step. 1 is pure time-lagging. */
  picard?: number;
  /** Tackley pseudo-plastic law instead of the power law, in the variable-μ tier. */
  tackley?: boolean;
  /** Constant ductile yield stress, Tackley law. */
  sigmaY?: number;
  /** Gradient of brittle yield stress with depth, Tackley law. */
  sigmaB?: number;
  /** Minimum plastic viscosity, Tackley law. */
  etaStar?: number;
}

export class Simulation {
  readonly geom: Geometry;
  readonly rAx: Axis;
  readonly aAx: Axis;
  /** Tier 1: the direct DFT solve. Null when μ varies in φ. */
  readonly stokes: StokesSolver | null;
  /** Tier 2: matrix-free PCG. Null for constant μ. */
  readonly variable: VariableStokes | null;
  readonly psi: Field;
  readonly temp: Temperature;
  readonly Ra: number;
  readonly gamma: number;
  readonly cz: number;
  readonly iters: number;
  readonly n: number;
  readonly picard: number;
  readonly tackley: boolean;
  readonly sigmaY: number;
  readonly sigmaB: number;
  readonly etaStar: number;
  private readonly cfl: number;
  private readonly dtMax: number;
  time = 0;

  constructor(o: SimOptions = {}) {
    const { nr = 32, na = 64, gnr = 65, gna = 128, geom = ANNULUS } = o;
    this.geom = geom;
    this.Ra = o.Ra ?? 1e4;
    this.gamma = o.gamma ?? 0;
    this.cz = o.cz ?? 0;
    this.iters = o.iters ?? 12;
    this.n = o.n ?? 1;
    this.picard = o.picard ?? 1;
    this.tackley = o.tackley ?? false;
    this.sigmaY = o.sigmaY ?? 1;
    this.sigmaB = o.sigmaB ?? 1;
    this.etaStar = o.etaStar ?? 1e-3;
    this.cfl = o.cfl ?? 1.0;
    this.dtMax = o.dtMax ?? 1e-3;
    this.rAx = clampedAxis(nr, geom.lo, geom.hi);
    this.aAx = periodicAxis(na, geom.span);
    // The tier is a property of the *mode*, not of γ: γ = 0 in the variable tier
    // is the isoviscous limiting case, and must still take the Krylov path or it
    // would verify nothing.
    const variable = o.variable ?? false;
    this.stokes = variable
      ? null : new StokesSolver(this.rAx, this.aAx, () => 1, false, geom);
    this.variable = variable
      ? new VariableStokes(this.rAx, this.aAx,
        this.tackley ? meanTackleyViscosity(geom) : meanViscosity(geom, this.gamma, this.cz),
        geom)
      : null;
    this.psi = new Field(this.rAx, this.aAx, geom);
    this.temp = new Temperature(geom, gnr, gna);
    if (o.seed) this.temp.reset(o.seed.amp ?? 0.05, o.seed.mode ?? 4);
    this.solveFlow();
  }

  /**
   * Steps 1–3: buoyancy → ψ.
   *
   * **The rheology is time-lagged**: μ is evaluated from the ψ
   * the previous solve left, which is what keeps the system CG sees linear and
   * symmetric even when the law is not. At n = 1 that lag is exact — μ does not
   * depend on ψ at all — and at n > 1 it is one Picard step. Asking for more
   * `picard` sweeps re-lags against the ψ just computed, at a full Krylov budget
   * each; one is the default because Stokes is quasi-static, so the previous
   * frame's strain rate is already an O(dt) guess at this one's.
   */
  private solveFlow(): void {
    const T = (r: number, phi: number) => this.temp.sample(this.temp.T, r, phi);
    const load = buoyancyLoad(this.rAx, this.aAx, T, this.Ra);
    if (!this.variable) {
      const c = this.stokes!.solve(load);
      for (let i = 0; i < this.rAx.n; i++) this.psi.c[i].set(c[i]);
      return;
    }
    for (let sweep = 0; sweep < this.picard; sweep++) {
      const e = strainRate(this.variable.tables, this.psi.c);
      // Tackley reads ε̇ raw: the yield stress is an absolute threshold, so
      // normalising it away (as the power law does, for scale-invariance)
      // would defeat the point.
      let mu: Float64Array;
      if (this.tackley) {
        mu = viscosityAt(this.variable.tables, T,
          (t, s, r) => tackleyViscosity(
            t, depthAt(this.geom, r), s, this.sigmaY, this.sigmaB, this.etaStar), e);
      } else {
        const { d, g } = strainScale(e);
        mu = viscosityAt(this.variable.tables, T,
          (t, s, r) => viscosity(t, this.gamma, (s + d) / g, this.n,
            depthAt(this.geom, r), this.cz), e);
      }
      // ψ is passed in as the initial guess and updated in place — see
      // `VariableStokes.solve` on why the previous frame is the right start.
      this.variable.solve(load, mu, this.psi.c, this.iters);
    }
  }

  get velocity(): Velocity {
    return (r, phi) => this.psi.velocity(r, phi);
  }

  /** Advance one step; returns the dt taken. */
  step(dt?: number): number {
    const u = this.velocity;
    const s = this.temp.maxSpeed(u);
    const h = dt ?? Math.min(this.dtMax, s > 0 ? this.cfl / s : this.dtMax);
    this.temp.advect(u, h);
    this.temp.diffuse(h);
    this.solveFlow();
    this.time += h;
    return h;
  }

  /** Run until the Nusselt number settles, or `maxSteps` is reached. */
  run(maxSteps = 2000, tol = 1e-6): { steps: number; nu: number; converged: boolean } {
    let prev = this.temp.nusselt().outer;
    for (let n = 1; n <= maxSteps; n++) {
      this.step();
      if (n % 20 === 0) {
        const nu = this.temp.nusselt().outer;
        if (Math.abs(nu - prev) < tol * Math.max(1, Math.abs(nu)))
          return { steps: n, nu, converged: true };
        prev = nu;
      }
    }
    return { steps: maxSteps, nu: this.temp.nusselt().outer, converged: false };
  }
}
