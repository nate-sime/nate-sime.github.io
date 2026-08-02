/**
 * Temperature-, depth- and strain-rate-dependent viscosity.
 *
 *   μ(T, d, ε̇) = clamp( exp(−γ(T − ½) + c(d − ½)) · ŝ^((1−n)/n),  μ_lo, μ_hi ),
 *   ŝ = (ε̇ + δ) / G                                       (`strainScale`)
 *
 * A Frank–Kamenetskii thermal term, an exponential depth term, and a regularised
 * power law, with a viscosity floor and ceiling. Five choices here are
 * load-bearing.
 *
 * **The exponents are centred on T = ½ and d = ½**, so the *geometric* mean of
 * each is exactly 1 across [0, 1]: raising a contrast stiffens the cold lid (or
 * the deep interior) and weakens the hot interior (or the shallow layer)
 * symmetrically instead of rescaling the whole problem, and the effective
 * Rayleigh number stays comparable to the isoviscous run at the same Ra. That is
 * what makes the sliders' effect legible.
 *
 * It is also what lets the Blankenbach variable-viscosity cases be run *exactly*
 * rather than approximately. Those state the law as μ = exp(−b T + c d) with the
 * reference at the cold surface, which is this one times the constant
 * `exp((γ − c)/2)`; a constant factor on μ is a rescaling of Ra and nothing else,
 * since the Stokes problem is linear in μ and the load is linear in Ra. So the
 * benchmark is reproduced by solving at `Ra · exp((γ − c)/2)` — see
 * `tests/blankenbach.test.ts`, where 2a and 2b are checked that way.
 *
 * **n = 1 removes the power law identically**, not approximately: the exponent
 * is then 0 and the factor is exactly 1 for any strain rate, while the thermal
 * and depth terms already lie inside [μ_lo, μ_hi] by construction so the clamp is
 * inactive too. μ(T, d, ε̇) at n = 1 *is* μ(T, d), bit for bit, so "n → 1
 * recovers the linear tier" is an exact identity rather than a limit — and the
 * two laws share one code path and one set of pipelines.
 * Likewise γ → 0 and c → 0 give μ ≡ 1 exactly, the isoviscous check, and c → 0
 * alone gives back the μ(T) law bit for bit.
 *
 * **Depth is the one dependence the preconditioner represents exactly.** μ̄(r)
 * is a radial profile (below), and `d` is a function of r alone — so unlike the
 * thermal term, which is only radial in the conductive state, and unlike the
 * power law, which has no radial profile at all, the depth term is carried by the
 * preconditioner with no approximation whatever the flow does. At γ = 0 that
 * makes μ̄ *the* operator rather than a spectral match to it, so the FFT solve of
 * tier 1 is exact and PCG converges in a single step at any depth contrast —
 * which is what `tests/rheology.test.ts` pins, and why raising this contrast does
 * not cost Krylov iterations the way raising the thermal one does.
 *
 * **The floor and ceiling are the interval the thermal and depth terms already
 * span**, `[e^{−(γ+c)/2}, e^{+(γ+c)/2}]`. The two exponents reach their extremes
 * at opposite corners — (T = 0, d = 1) and (T = 1, d = 0) — so that interval is
 * attained and never exceeded, and "contrast" means the *product* of the two the
 * sliders ask for. The power law then redistributes μ within that range rather
 * than widening it, so the solver never sees a contrast the preconditioner was
 * not sized for. Without a ceiling the `(1−n)/n < 0` exponent makes μ unbounded
 * as ε̇ → 0; `ε̇_min` regularises that, and the clamp bounds what the
 * regularisation leaves. At n = 3 the clamp is active on ~30% of the quadrature
 * points, which is the floor and ceiling doing their job on the tails of a field
 * that spans decades, not a sign they are set too tight.
 *
 * **The strain rate is normalised by the flow's own geometric mean**
 * (`strainScale`), not by a fixed material constant. Two things force that, and
 * the second is the one that was measured rather than assumed.
 *
 * *Scale.* With a fixed reference the prefactor scales like the flow amplitude
 * to the −(n−1)/n, so sweeping Ra across the four decades the slider offers
 * would drive μ into one clamp or the other and the rheology would read as
 * constant almost everywhere. Normalising makes the law scale-invariant in ψ, so
 * `n` controls how sharply μ localises into weak shear zones and strong stagnant
 * regions — the phenomenon the mode exists to show — while Ra keeps meaning Ra.
 * What is given up is the absolute coupling between flow speed and viscosity
 * level: a stronger flow does not become globally weaker, only more *unevenly*
 * weak. That is the same trade the T = ½ centring makes, for the same reason.
 *
 * *Centre.* Which average is used is not a detail. The exponent is negative, so
 * μ is a convex function of ε̇ and dividing by the **RMS** — the obvious first
 * choice — biases the factor upward: ε̇ is peaked, most of the domain sits well
 * below its RMS, and almost everywhere gets stiffened. Measured at Ra = 2×10⁴,
 * n = 3 and a 10² contrast, that moved ⟨log μ⟩ from the isoviscous run's 0.21 to
 * 0.79 and collapsed max |u| from 27 to 6.6 — the power law switching the
 * convection off, which is precisely the silent rescaling of the effective Ra
 * that the T = ½ centring exists to prevent. Dividing by the **geometric** mean
 * instead makes ⟨log factor⟩ zero by construction: measured ⟨log μ⟩ 0.19 against
 * the same 0.21, with max |u| rising to 51 — shear-thinning speeding the flow
 * up, which is the right sign and the right magnitude.
 */

import type { Geometry } from "../geometry";

/**
 * Regularisation of the strain rate, relative to its own RMS. Bounds μ where the
 * flow is momentarily stagnant, and keeps `log ε̇` finite at a stagnation point
 * so the geometric mean is well defined.
 */
export const EPS_MIN = 1e-3;

/**
 * μ(T, ŝ, d). `strain` is the **normalised** strain rate `ŝ` from `strainScale`;
 * the default 1 is the unstrained case, and `n = 1` the μ(T, d) law unchanged.
 * `depth` is the normalised depth from `depthAt`, and its default ½ is the
 * neutral point of the depth term — so at `cz = 0` the argument is irrelevant and
 * every μ(T) call site reads unchanged.
 *
 * T is clamped to [0, 1] first. The monotone limiter and the isothermal
 * boundaries already bound T, but only to within the advection scheme's own
 * round-off; without the clamp a stray 1 + 1e-7 would be exponentiated, and at
 * γ = 9 that is a silent perturbation of the coefficient rather than an error.
 * `depth` needs no such guard: it is a function of the quadrature abscissa, not
 * of the solution, so it is in [0, 1] by construction.
 */
export const viscosity = (
  T: number, gamma: number, strain = 1, n = 1, depth = 0.5, cz = 0,
): number => {
  const hi = Math.exp(0.5 * (gamma + cz));
  const mu = Math.exp(-gamma * (Math.min(1, Math.max(0, T)) - 0.5) + cz * (depth - 0.5))
    * strain ** ((1 - n) / n);
  return Math.min(hi, Math.max(1 / hi, mu));
};

/**
 * γ from the viscosity contrast μ(0)/μ(1) — the quantity the UI exposes. The
 * depth contrast μ|_{d=1}/μ|_{d=0} maps to c through the same logarithm, which is
 * why there is one function and not two.
 */
export const gammaFor = (contrast: number): number => Math.log(contrast);

/**
 * Normalised depth: 0 at the cold boundary, 1 at the hot one.
 *
 * The solver's radial coordinate runs the other way — `lo` is the hot inner
 * radius, or the floor of the box — so this is `(hi − r)/(hi − lo)` and not the
 * axis coordinate. Depth *below the cold surface* is the variable a viscosity
 * profile is stated in, in the literature and in the benchmarks, and it is the
 * one that means the same thing in both geometries: in a box it is exactly
 * `1 − z`, and on the annulus it is the fractional depth through the shell.
 */
export const depthAt = (g: Geometry, r: number): number =>
  (g.hi - r) / (g.hi - g.lo);

/**
 * The normalisation of the strain-rate field: the regularisation offset
 * `δ = ε̇_min · rms(ε̇)` and the geometric mean `G` of `ε̇ + δ`. Normalised
 * strain is then `ŝ = (ε̇ + δ)/G`, whose geometric mean is exactly 1 — so the
 * power-law factor has geometric mean exactly 1 too, whatever `n` is, and the
 * law redistributes μ without shifting it.
 *
 * Both scales are relative, so ψ can be scaled freely; a *zero* field — the
 * state the very first solve is warm-started from — has no scale at all and
 * returns `(1, 1)`, making ŝ ≡ 1 and that solve a plain μ(T) one. It is not a
 * fallback that needs to be good: the next step has a flow.
 */
export const strainScale = (e: Float64Array): { d: number; g: number } => {
  let s = 0;
  for (let i = 0; i < e.length; i++) s += e[i] * e[i];
  const rms = Math.sqrt(s / e.length);
  if (!(rms > 0)) return { d: 1, g: 1 };
  const d = EPS_MIN * rms;
  let l = 0;
  for (let i = 0; i < e.length; i++) l += Math.log(e[i] + d);
  return { d, g: Math.exp(l / e.length) };
};

/**
 * μ̄(r) for the FFT preconditioner (tier 2): the azimuthal mean
 * viscosity, evaluated on the geometry's steady conduction profile —
 * `ln(r_o/r)/ln(r_o/r_i)` on the annulus, `1 − z` in a box — and on the depth
 * term, which is radial exactly and needs no profile assumed for it at all.
 *
 * **This is deliberately a fixed profile, not the running mean.** Rebuilding the
 * preconditioner means re-inverting one dense radial block per azimuthal mode in
 * f64 — the same O(n³) job that dominates start-up — so it cannot
 * happen per frame, and a preconditioner does not need to: it only has to be
 * spectrally close, and it is refreshed whenever γ or c changes, which is a UI
 * event.
 *
 * The *thermal* half of the profile is exact at t = 0 and through the subcritical
 * regime. Once convection develops, the interior mixes towards T ≈ ½ (where the
 * thermal term ≈ 1) and the variation concentrates into thin boundary layers, so
 * μ̄ overstates the radial spread and the iteration count rises — which is the
 * honest failure mode for a preconditioner, and is what the fixed iteration
 * budget is sized against. The *depth* half is under no such assumption: it
 * depends on r and nothing else, so it is exact for every state of the flow, and
 * a run at γ = 0 with any depth contrast is preconditioned by its own operator.
 *
 * The power law is deliberately **not** represented here. It has no radial
 * profile to average — its variation is where the flow deforms, which is
 * azimuthal as much as radial — and the clamp keeps μ inside the same interval
 * the thermal term spans, so μ̄ remains as good a spectral match at n = 3 as it
 * is at n = 1 (measured).
 */
export const meanViscosity = (
  geom: Geometry, gamma: number, cz = 0,
): ((r: number) => number) => (r) =>
  viscosity(geom.conduction(r), gamma, 1, 1, depthAt(geom, r), cz);

/**
 * Tackley (2000) pseudo-plastic rheology: an Arrhenius (diffusion-creep) branch
 * in parallel with a Bingham (yielding) branch, so the weaker of the two sets
 * the local viscosity — the mechanism that localises the cold lid into
 * plate-like shear zones rather than deforming it uniformly.
 *
 *   η_lin(T, d)   = A₀(d) exp(27.631/(T + 0.88)) · 5.86052e-13
 *   η_plast(d, ε̇) = η* + (σ_Y + σ_b d)/ε̇
 *   η              = (η_lin⁻¹ + η_plast⁻¹)⁻¹
 *
 * A₀ steps from 1 to 30 at the 670 km discontinuity — upper/lower mantle — and
 * is a property of depth alone, not a slider. σ_Y, σ_b and η* are.
 */

/** Nondimensional depth of the 660/670 km discontinuity (mantle thickness 2885 km). */
export const TACKLEY_TRANSITION_DEPTH = 670 / 2885;

/** Floor under ε̇ in the plastic branch — a stress divided by zero strain rate is not a physical limit, just an unset one. */
export const TACKLEY_STRAIN_FLOOR = 1e-8;

/** η_lin(T, d): activation energy 27.631, offset 0.88, prefactor 5.86052e-13 — nondimensional Arrhenius diffusion creep, stiffened 30× below 670 km. */
export const tackleyLinear = (T: number, depth: number): number => {
  const A0 = depth < TACKLEY_TRANSITION_DEPTH ? 1 : 30;
  const Tc = Math.min(1, Math.max(0, T));
  return A0 * Math.exp(27.631 / (Tc + 0.88)) * 5.86052e-13;
};

/** η_plast(d, ε̇): ductile floor η* in series with a depth-dependent yield stress over the strain rate. */
export const tackleyPlastic = (
  depth: number, strain: number, sigmaY: number, sigmaB: number, etaStar: number,
): number =>
  etaStar + (sigmaY + sigmaB * depth) / Math.max(strain, TACKLEY_STRAIN_FLOOR);

/** η(T, d, ε̇): the two branches as resistors in parallel. */
export const tackleyViscosity = (
  T: number, depth: number, strain: number,
  sigmaY: number, sigmaB: number, etaStar: number,
): number => {
  const lin = tackleyLinear(T, depth);
  const plast = tackleyPlastic(depth, strain, sigmaY, sigmaB, etaStar);
  return 1 / (1 / lin + 1 / plast);
};

/**
 * μ̄(r) for the FFT preconditioner: the ε̇ → 0 limit, where the plastic
 * branch's resistance diverges and only η_lin remains — independent of
 * σ_Y, σ_b, η*, so changing them is a pure uniform write, no re-inversion.
 */
export const meanTackleyViscosity = (geom: Geometry): ((r: number) => number) =>
  (r) => tackleyLinear(geom.conduction(r), depthAt(geom, r));

/**
 * Brandenburg pseudo-plastic rheology: the same parallel combination as
 * Tackley's — a linear branch and a Bingham (yielding) branch, the weaker
 * setting η — but with the linear branch replaced by the μ(T, d) exponential
 * law of this file's own `viscosity`, times a depth-stepped prefactor A₀
 * instead of Tackley's fixed Arrhenius form:
 *
 *   η_lin(T, d)   = A₀(d) exp(−b(T − ½) + c(d − ½))
 *   η_plast(d, ε̇) = η* + (σ_Y + σ_b d)/ε̇
 *   η              = (η_lin⁻¹ + η_plast⁻¹)⁻¹
 *
 * η_plast is `tackleyPlastic` unchanged — the two laws state the identical
 * yielding branch, so there is one implementation of it, not two. A₀ steps
 * between `aUpper` and `aLower` at a threshold depth `d0`, all three supplied
 * by the caller rather than fixed the way Tackley's 1/30 at 670 km are.
 */

/** A₀(d): `aUpper` above the threshold depth `d0`, `aLower` below it. */
export const brandenburgA0 = (depth: number, aUpper: number, aLower: number, d0: number): number =>
  depth < d0 ? aUpper : aLower;

/** η_lin(T, d): the μ(T, d) exponential, unclamped, times the depth-stepped prefactor A₀. */
export const brandenburgLinear = (
  T: number, depth: number, b: number, c: number, aUpper: number, aLower: number, d0: number,
): number => {
  const Tc = Math.min(1, Math.max(0, T));
  return brandenburgA0(depth, aUpper, aLower, d0) * Math.exp(-b * (Tc - 0.5) + c * (depth - 0.5));
};

/** η(T, d, ε̇): the two branches as resistors in parallel. */
export const brandenburgViscosity = (
  T: number, depth: number, strain: number, b: number, c: number,
  aUpper: number, aLower: number, d0: number,
  sigmaY: number, sigmaB: number, etaStar: number,
): number => {
  const lin = brandenburgLinear(T, depth, b, c, aUpper, aLower, d0);
  const plast = tackleyPlastic(depth, strain, sigmaY, sigmaB, etaStar);
  return 1 / (1 / lin + 1 / plast);
};

/**
 * μ̄(r) for the FFT preconditioner: the ε̇ → 0 limit, where only η_lin
 * remains — independent of σ_Y, σ_b, η* exactly as Tackley's is, but *not*
 * independent of b, c, aUpper, aLower or d0, so changing any of those five
 * re-inverts the radial blocks (see `GpuSimulation.setContrast` and
 * `.setBrandenburgProfile`).
 */
export const meanBrandenburgViscosity = (
  geom: Geometry, b: number, c: number, aUpper: number, aLower: number, d0: number,
): ((r: number) => number) =>
  (r) => brandenburgLinear(geom.conduction(r), depthAt(geom, r), b, c, aUpper, aLower, d0);
