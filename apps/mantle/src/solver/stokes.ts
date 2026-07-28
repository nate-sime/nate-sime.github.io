/**
 * Stokes solvers, both tiers of the cost ladder. CPU reference, f64.
 *
 *   `StokesSolver`  — DFT in φ + per-mode radial direct solve. The whole solve
 *                     when μ depends on r alone (tier 1, constant μ), and the
 *                     **preconditioner** when it does not (tier 2). One kernel,
 *                     two jobs — the unification the whole scheme is built around.
 *   `VariableStokes` — matrix-free preconditioned conjugate gradients on the true
 *                     variable-μ operator, using the above as M⁻¹ (tier 2).
 *
 * Per azimuthal mode k the radial operator is
 *
 *   A_k = 4 S1(k) R1 + S0(k) R2 + S1(k) R3 + S2(k) R4,
 *
 * assembled and factorised **once in f64**: it is time-invariant, so every
 * subsequent solve is an application only. The GPU path uploads the dense
 * inverses of the same blocks — see `modeInverses` for why inverses rather
 * than factors, and `A_k = A_{n−k}` for why only half the modes are stored.
 *
 * Boundary conditions: `ψ = const` on each radius, imposed by
 * dropping the first and last radial DOF (open knots make them interpolatory).
 * `σ_rφ = 0` is natural to the dissipation form. That elimination also removes
 * the k = 0 kernel outright — see `modeInverses` — so the mean mode is an
 * ordinary nonsingular block, solved when `k0` is set and skipped when it is
 * not (constant μ leaves it unforced).
 */

import { ANNULUS, type Geometry } from "../geometry";
import { Axis, P } from "../spline";
import { mat, lu, solve, type LU } from "../linalg";
import { gauss } from "../quad";
import { radialBlocks, azimuthalSymbols, radialOperator } from "./operators";
import { applyOperator, operatorTables, type OperatorTables } from "./assembly";

export class StokesSolver {
  readonly nr: number;
  readonly na: number;
  private readonly fac: (LU | null)[]; // k = 0 … na/2; see radialOperator

  /**
   * `mu` is a radial profile: the constant μ for tier 1, the azimuthal mean
   * μ̄(r) when this object is a tier-2 preconditioner. `k0` solves the azimuthal
   * mean mode instead of zeroing it.
   */
  constructor(
    readonly rAx: Axis, readonly aAx: Axis,
    mu: (r: number) => number = () => 1, k0 = false,
    readonly geom: Geometry = ANNULUS,
  ) {
    this.nr = rAx.n;
    this.na = aAx.n;
    const R = radialBlocks(rAx, mu, geom), S = azimuthalSymbols(aAx);
    this.fac = Array.from({ length: this.na / 2 + 1 }, (_, k) =>
      k === 0 && !k0 ? null : lu(radialOperator(R, S, k)));
  }

  /** Solve a(ψ, v) = ⟨load, v⟩ for the spline coefficients of ψ. */
  solve(load: Float64Array[]): Float64Array[] {
    const { nr, na } = this, ni = nr - 2;
    const re = mat(nr, na), im = mat(nr, na);
    for (let i = 0; i < nr; i++)
      for (let k = 0; k < na; k++) {
        let sr = 0, si = 0;
        for (let j = 0; j < na; j++) {
          const t = (-2 * Math.PI * j * k) / na;
          sr += load[i][j] * Math.cos(t);
          si += load[i][j] * Math.sin(t);
        }
        re[i][k] = sr / na;
        im[i][k] = si / na;
      }

    const xr = mat(nr, na), xi = mat(nr, na);
    for (let k = 0; k < na; k++) {
      const f = this.fac[Math.min(k, na - k)];
      if (!f) continue;
      const br = new Float64Array(ni), bi = new Float64Array(ni);
      for (let i = 0; i < ni; i++) { br[i] = re[i + 1][k]; bi[i] = im[i + 1][k]; }
      const sr = solve(f, br), si = solve(f, bi);
      for (let i = 0; i < ni; i++) { xr[i + 1][k] = sr[i]; xi[i + 1][k] = si[i]; }
    }

    const psi = mat(nr, na);
    for (let i = 0; i < nr; i++)
      for (let j = 0; j < na; j++) {
        let s = 0;
        for (let k = 0; k < na; k++) {
          const t = (2 * Math.PI * j * k) / na;
          s += xr[i][k] * Math.cos(t) - xi[i][k] * Math.sin(t);
        }
        psi[i][j] = s;
      }
    return psi;
  }
}

/** Load vector ⟨S, v⟩ = ∫∫ S B_m N_l h dr dφ. */
export function loadVector(
  rAx: Axis, aAx: Axis, S: (r: number, phi: number) => number,
  geom: Geometry = ANNULUS,
): Float64Array[] {
  const b = mat(rAx.n, aAx.n);
  const aq = aAx.elements().flatMap(([a, c]) =>
    gauss(a, c).map(({ x, w }) => ({ x, w, ...aAx.ders(x, 0) })));

  for (const [ea, eb] of rAx.elements())
    for (const { x: r, w: wr } of gauss(ea, eb)) {
      const R = rAx.ders(r, 0);
      for (const A of aq) {
        const f = S(r, A.x) * wr * A.w * geom.h(r);
        for (let p = 0; p <= P; p++) {
          const i = rAx.dof(R.span, p), c = f * R.N[0][p];
          for (let q = 0; q <= P; q++) b[i][aAx.dof(A.span, q)] += c * A.N[0][q];
        }
      }
    }
  return b;
}

// ---- tier 2: variable μ ------------------------------------------------------

const dot = (a: Float64Array[], b: Float64Array[]): number => {
  let s = 0;
  for (let i = 0; i < a.length; i++)
    for (let j = 0; j < a[i].length; j++) s += a[i][j] * b[i][j];
  return s;
};

/**
 * Variable-μ Stokes by preconditioned conjugate gradients (tier 2).
 *
 * CG rather than a flexible GMRES: the dissipation form
 * is symmetric positive definite on the constrained space, and the
 * preconditioner — the μ̄(r) Galerkin operator inverted exactly per mode — is
 * SPD and *fixed*, so nothing here needs the extra machinery. Flexible variants
 * exist to tolerate a preconditioner that changes between iterations; this one
 * does not.
 *
 * **Warm start.** `solve` iterates from the ψ it is handed rather than from
 * zero. Stokes is quasi-static, so consecutive frames differ only by the O(dt)
 * change in T, and the previous ψ is an excellent guess — that is what makes a
 * small fixed budget (rather than an iteration-to-tolerance loop, which would
 * make frame time depend on the flow) enough to stay converged. The price is one
 * extra operator apply per frame to form the initial residual.
 *
 * **Fixed budget, no convergence test.** The GPU twin cannot test a residual
 * without a readback in the hot loop, so the CPU reference does not
 * either: both run exactly `iters` iterations, which is also what makes them
 * comparable step for step. `residual` exists for the verification suite and is
 * never called from the stepping path.
 */
export class VariableStokes {
  readonly tables: OperatorTables;
  readonly pre: StokesSolver;

  constructor(
    readonly rAx: Axis, readonly aAx: Axis, muBar: (r: number) => number,
    readonly geom: Geometry = ANNULUS,
  ) {
    this.tables = operatorTables(rAx, aAx, geom);
    this.pre = new StokesSolver(rAx, aAx, muBar, true, geom);
  }

  /** A ψ, with μ given at the tensor grid of quadrature points. */
  apply(c: Float64Array[], mu: Float64Array): Float64Array[] {
    return applyOperator(this.tables, this.rAx.n, this.aAx.n, c, mu, this.geom);
  }

  /**
   * ‖b − Aψ‖₂, **recomputed**. Deliberately a separate call rather than
   * something `solve` returns, for two reasons that point the same way.
   *
   * It must never be read off CG's recursion. The iteration carries
   * `r ← r − αAp`, and in finite precision that drifts below the true residual
   * without bound once the true one stalls at the arithmetic's floor — in f32 it
   * reaches 1e-17 against a real 1e-3, so a convergence claim made from it is not
   * a claim about the answer. Recomputing means an extra operator apply.
   *
   * And that apply is pure verification: the coupled loop wants ψ, not a
   * diagnostic. Charging it only to the caller who asks keeps it off the
   * reference's stepping path, where it is otherwise paid every step and read by
   * nobody.
   */
  residual(load: Float64Array[], mu: Float64Array, x: Float64Array[]): number {
    const nr = this.rAx.n, na = this.aAx.n, Ax = this.apply(x, mu);
    let s = 0;
    for (let i = 1; i < nr - 1; i++)
      for (let j = 0; j < na; j++) s += (load[i][j] - Ax[i][j]) ** 2;
    return Math.sqrt(s);
  }

  /** PCG. `x` is both the initial guess and the output — updated in place. */
  solve(
    load: Float64Array[], mu: Float64Array, x: Float64Array[], iters: number,
  ): void {
    const nr = this.rAx.n, na = this.aAx.n;

    // The load's boundary rows are assembled but not in the trial space; leaving
    // them in would give the residual a component the operator can never remove.
    const r = mat(nr, na);
    const Ax = this.apply(x, mu);
    for (let i = 1; i < nr - 1; i++)
      for (let j = 0; j < na; j++) r[i][j] = load[i][j] - Ax[i][j];

    let z = this.pre.solve(r);
    const p = z.map((row) => Float64Array.from(row));
    let rz = dot(r, z);

    for (let n = 0; n < iters; n++) {
      const Ap = this.apply(p, mu);
      const pAp = dot(p, Ap);
      // The Krylov space is spent. Reachable rather than defensive: with γ = 0
      // the preconditioner *is* the operator, so the residual collapses after
      // one iteration. In f64 it lands at round-off rather than exactly zero and
      // the loop usually runs on harmlessly; on the GPU, in f32, the same
      // situation does produce a zero denominator, which is why the update
      // kernels carry the same guard.
      if (!(pAp > 0) || !(rz > 0)) break;
      const alpha = rz / pAp;
      for (let i = 1; i < nr - 1; i++)
        for (let j = 0; j < na; j++) {
          x[i][j] += alpha * p[i][j];
          r[i][j] -= alpha * Ap[i][j];
        }
      z = this.pre.solve(r);
      const rzNew = dot(r, z), beta = rzNew / rz;
      for (let i = 1; i < nr - 1; i++)
        for (let j = 0; j < na; j++) p[i][j] = z[i][j] + beta * p[i][j];
      rz = rzNew;
    }
  }
}
