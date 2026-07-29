/**
 * Separable Galerkin operators for constant-viscosity Stokes.
 *
 * The bilinear form is the viscous dissipation
 *
 *   a(ψ, v) = ∫ 2μ ε(u[ψ]) : ε(u[v]) dx,     u[·] = ∇×(· ẑ),
 *
 * which is the correct variational statement for free-slip: `ψ = const` is
 * essential (imposed by dropping the boundary DOFs) and `σ_rφ = 0` — the
 * condition a curved boundary actually demands — is *natural*, so no boundary
 * term need be derived. Pressure is absent because `∇·u[v] ≡ 0` pointwise.
 *
 * With `u_r = ψ_φ/h`, `u_φ = −ψ_r` and incompressibility (`ε_φφ = −ε_rr`),
 *
 *   2ε:ε = 4(ε_rr² + ε_rφ²),
 *   ε_rr = ∂_φ A[ψ],   A[ψ] = ψ_r/h − h′ψ/h²,
 *   ε_rφ = −½ C[ψ],    C[ψ] = ψ_rr − (h′/h) ψ_r − ψ_φφ/h²,
 *
 * `h` being the transverse metric of `geometry.ts` — `r` on the annulus, where
 * these are the familiar curvilinear expressions, and `1` in a box, where the
 * `h′` terms drop and they collapse to `ψ_rz` and `ψ_zz − ψ_xx`. Both are the
 * *same* form with a different Jacobian, which is why one assembly serves both.
 *
 * Every radial coefficient is φ-independent, so in the tensor basis
 * `ψ = Σ c_ij B_i(r) N_j(φ)` the form separates into four radial blocks paired
 * with the azimuthal Gram matrices ∫N_jN_l, ∫N_j'N_l', ∫N_j''N_l''. Uniform
 * periodic knots make those circulant, hence exactly diagonalised by the DFT
 * — the per-mode multipliers are the *discrete symbols* below, not
 * the analytic k², k⁴. Being read off the knot vector, they carry the period of
 * φ with them and need no correction for a box of arbitrary length.
 */

import { ANNULUS, type Geometry } from "../geometry";
import { Axis, P } from "../spline";
import { mat, lu, solve } from "../linalg";
import { gauss } from "../quad";

/**
 * R1..R4: radial blocks pairing with symbols S1, S0, S1, S2 respectively.
 *
 * `mu` is a **radial** viscosity profile, integrated inside the quadrature
 * rather than multiplying the assembled blocks. For tier 1 it is the constant
 * μ; for tiers 2–3 it is the azimuthal mean μ̄(r) of `rheology.ts`, which is the
 * only extra ingredient the preconditioner needs — every coefficient stays
 * φ-independent, so the circulant structure (and with it the exact DFT
 * decoupling) survives untouched. That is what lets one kernel serve as both
 * the whole solve and the preconditioner.
 *
 * The area element is `h dr dφ`, and `b = B/h²` is the coefficient pairing with
 * −N″ — so in a box the four integrands are `B′B′`, `B″B″`, `B″B + BB″` and
 * `BB`, with no negative power of a coordinate anywhere. That is what lets the
 * box put its lower boundary at z = 0, which the annulus form could not.
 */
export function radialBlocks(
  ax: Axis, mu: (r: number) => number = () => 1, geom: Geometry = ANNULUS,
): Float64Array[][] {
  const n = ax.n, R = [0, 1, 2, 3].map(() => mat(n, n));
  for (const [ea, eb] of ax.elements())
    for (const { x: r, w: w0 } of gauss(ea, eb)) {
      const { span, N } = ax.ders(r, 2);
      const J = geom.h(r), ih = 1 / J, ih2 = ih * ih, dh = geom.dh;
      const w = w0 * mu(r);
      for (let p = 0; p <= P; p++) {
        const i = ax.dof(span, p);
        const b = N[0][p] * ih2;
        const a = N[1][p] * ih - dh * N[0][p] * ih2, g = N[2][p] - dh * N[1][p] * ih;
        for (let q = 0; q <= P; q++) {
          const m = ax.dof(span, q);
          const bq = N[0][q] * ih2;
          const aq = N[1][q] * ih - dh * N[0][q] * ih2;
          const gq = N[2][q] - dh * N[1][q] * ih;
          R[0][i][m] += w * a * aq * J;                  // ∫ μ a_i a_m h dr
          R[1][i][m] += w * g * gq * J;                  // ∫ μ g_i g_m h dr
          R[2][i][m] += w * (g * bq + b * gq) * J;       // ∫ μ (g_ib_m + b_ig_m) h dr
          R[3][i][m] += w * b * bq * J;                  // ∫ μ b_i b_m h dr
        }
      }
    }
  return R;
}

/** DFT symbols of the circulant azimuthal Gram matrices [S0, S1, S2]. */
export function azimuthalSymbols(ax: Axis): Float64Array[] {
  const n = ax.n, row = [0, 1, 2].map(() => new Float64Array(n));
  for (const [ea, eb] of ax.elements())
    for (const { x, w } of gauss(ea, eb)) {
      const { span, N } = ax.ders(x, 2);
      for (let p = 0; p <= P; p++) {
        if (ax.dof(span, p) !== 0) continue; // circulant: first row suffices
        for (let q = 0; q <= P; q++) {
          const l = ax.dof(span, q);
          for (let d = 0; d < 3; d++) row[d][l] += w * N[d][p] * N[d][q];
        }
      }
    }
  // Real symmetric circulant ⇒ real symbols λ_k = Σ_m c_m cos(2πkm/n).
  return row.map((c) => {
    const s = new Float64Array(n);
    for (let k = 0; k < n; k++) {
      let t = 0;
      for (let m = 0; m < n; m++) t += c[m] * Math.cos((2 * Math.PI * k * m) / n);
      s[k] = t;
    }
    return s;
  });
}

/**
 * Interior radial block for azimuthal mode k:
 *
 *   A_k = 4 S1(k) R1 + S0(k) R2 + S1(k) R3 + S2(k) R4,
 *
 * with the first and last radial DOF dropped — open knots make them
 * interpolatory, so that is exactly the essential condition ψ = 0.
 * The viscosity is already inside R (see `radialBlocks`).
 *
 * The symbols are even (`S[k] = S[n−k]`, being cosine sums), so `A_k = A_{n−k}`
 * and only `k = 0 … n/2` are distinct. Both the CPU factorisations and the GPU
 * inverse table are built over that half range and indexed by `min(k, n−k)`.
 *
 * **k = 0 is not singular in this space** — see `modeInverses`.
 */
export function radialOperator(
  R: Float64Array[][], S: Float64Array[], k: number,
): Float64Array[] {
  const [R1, R2, R3, R4] = R, [S0, S1, S2] = S;
  const ni = R1.length - 2, A = mat(ni, ni);
  for (let i = 0; i < ni; i++)
    for (let j = 0; j < ni; j++)
      A[i][j] = 4 * S1[k] * R1[i + 1][j + 1] + S0[k] * R2[i + 1][j + 1]
        + S1[k] * R3[i + 1][j + 1] + S2[k] * R4[i + 1][j + 1];
  return A;
}

/**
 * Dense per-mode inverses `A_k⁻¹`, k = 0 … n/2, flattened for GPU upload.
 * The GPU applies these as a matvec in f32.
 *
 * *Why inverses and not the LU factors.* Both are computed here in f64, so both
 * are accurate; the difference is what the f32 *application* costs. Triangular
 * substitution against f32 factors reintroduces the operator's κ ~ h⁻⁴ into the
 * hot loop, whereas a matvec is backward stable with respect to the stored
 * matrix — the conditioning is spent once, at init, in f64. A matvec is also
 * perfectly parallel, while substitution is sequential in r.
 *
 * **`k0` — the azimuthal mean mode.** In the unconstrained space the k = 0 block
 * is rank-deficient by two, its kernel spanned by the gauge constant and by rigid
 * rotation `ψ = −Ωr²/2` (stress-free, hence genuinely free-slip). Both are
 * already removed here, by the *essential* condition rather than by any explicit
 * constraint: dropping both boundary DOFs pins ψ(r_i) = ψ(r_o) = 0, and a
 * parabola vanishing at two distinct radii is zero. (In a box the same two modes
 * are the constant and `ψ = −Uz`, uniform horizontal translation, and a linear
 * function vanishing at z = 0 and z = 1 is zero — the argument is the same one.)
 * The block is measurably as well conditioned as k = 1, so no angular-momentum
 * side condition is needed — for the operator or the preconditioner.
 *
 * That choice is not arbitrary: ψ(r_o) − ψ(r_i) is the net transverse volume
 * flux, so fixing it at zero is the gauge "no net circulation". It is consistent
 * because the buoyancy load is *exactly* orthogonal to that mode — its velocity
 * is purely transverse and the load pairs only with the component along gravity
 * — so the physics never asks for what is being suppressed.
 *
 * Tier 1 still passes `k0 = false`: with constant μ the modes do not couple and
 * the k = 0 forcing vanishes identically (∮∂_φT dφ = 0), so ψ̂₀ ≡ 0 and solving
 * it would only amplify f32 transform noise by κ(A_0). Tier 2 passes `k0 = true`
 * — variable μ couples the modes, so the residual the preconditioner is handed
 * has genuine mean-mode content.
 */
export function modeInverses(
  rAx: Axis, aAx: Axis, mu: (r: number) => number = () => 1, k0 = false,
  geom: Geometry = ANNULUS,
): Float32Array {
  const R = radialBlocks(rAx, mu, geom), S = azimuthalSymbols(aAx);
  const ni = rAx.n - 2, nk = aAx.n / 2 + 1;
  const out = new Float32Array(nk * ni * ni);
  const e = new Float64Array(ni);
  for (let k = k0 ? 0 : 1; k < nk; k++) {
    const f = lu(radialOperator(R, S, k));
    for (let c = 0; c < ni; c++) {
      e.fill(0); e[c] = 1;
      const col = solve(f, e); // = A⁻¹ e_c
      for (let a = 0; a < ni; a++) out[k * ni * ni + a * ni + c] = col[a];
    }
  }
  return out;
}
