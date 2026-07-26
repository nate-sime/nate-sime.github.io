/**
 * Temperature- and strain-rate-dependent viscosity.
 *
 *   μ(T, ε̇) = clamp( exp(−γ(T − ½)) · ŝ^((1−n)/n),  μ_lo, μ_hi ),
 *   ŝ = (ε̇ + δ) / G                                       (`strainScale`)
 *
 * A Frank–Kamenetskii thermal term times a regularised power law, with a
 * viscosity floor and ceiling. Four choices here are load-bearing.
 *
 * **The exponent is centred on T = ½**, so the *geometric* mean of the thermal
 * term is exactly 1 across T ∈ [0, 1]: raising the contrast stiffens the cold lid
 * and weakens the hot interior symmetrically instead of rescaling the whole
 * problem, and the effective Rayleigh number stays comparable to the isoviscous
 * run at the same Ra. That is what makes the slider's effect legible.
 *
 * **n = 1 removes the power law identically**, not approximately: the exponent
 * is then 0 and the factor is exactly 1 for any strain rate, while the thermal
 * term already lies inside [μ_lo, μ_hi] by construction so the clamp is inactive
 * too. μ(T, ε̇) at n = 1 *is* μ(T), bit for bit, so "n → 1 recovers the linear
 * tier" is an exact identity rather than a limit — and the two laws share one
 * code path and one set of pipelines.
 * Likewise γ → 0 gives μ ≡ 1 exactly, the isoviscous check.
 *
 * **The floor and ceiling are the interval the thermal term already spans**,
 * `[e^{−γ/2}, e^{+γ/2}]`. So "contrast" means the total viscosity contrast under
 * *both* laws — the power law redistributes μ within the range the user asked
 * for rather than widening it — and the solver never sees a contrast the
 * preconditioner was not sized for. Without a ceiling the `(1−n)/n < 0` exponent
 * makes μ unbounded as ε̇ → 0; `ε̇_min` regularises that, and the clamp bounds
 * what the regularisation leaves. At n = 3 the clamp is active on ~30% of the
 * quadrature points, which is the floor and ceiling doing their job on the tails
 * of a field that spans decades, not a sign they are set too tight.
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

/**
 * Regularisation of the strain rate, relative to its own RMS. Bounds μ where the
 * flow is momentarily stagnant, and keeps `log ε̇` finite at a stagnation point
 * so the geometric mean is well defined.
 */
export const EPS_MIN = 1e-3;

/**
 * μ(T, ŝ). `strain` is the **normalised** strain rate `ŝ` from `strainScale`;
 * the default 1 is the unstrained case, and `n = 1` the μ(T) law unchanged.
 *
 * T is clamped to [0, 1] first. The monotone limiter and the isothermal
 * boundaries already bound T, but only to within the advection scheme's own
 * round-off; without the clamp a stray 1 + 1e-7 would be exponentiated, and at
 * γ = 9 that is a silent perturbation of the coefficient rather than an error.
 */
export const viscosity = (
  T: number, gamma: number, strain = 1, n = 1,
): number => {
  const hi = Math.exp(0.5 * gamma);
  const mu = Math.exp(-gamma * (Math.min(1, Math.max(0, T)) - 0.5))
    * strain ** ((1 - n) / n);
  return Math.min(hi, Math.max(1 / hi, mu));
};

/** γ from the viscosity contrast μ(0)/μ(1) — the quantity the UI exposes. */
export const gammaFor = (contrast: number): number => Math.log(contrast);

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
 * viscosity, evaluated on the steady conduction profile
 * `T = ln(r_o/r)/ln(r_o/r_i)`.
 *
 * **This is deliberately a fixed profile, not the running mean.** Rebuilding the
 * preconditioner means re-inverting one dense radial block per azimuthal mode in
 * f64 — the same O(n³) job that dominates start-up — so it cannot
 * happen per frame, and a preconditioner does not need to: it only has to be
 * spectrally close, and it is refreshed whenever γ changes, which is a UI event.
 *
 * The profile is exact at t = 0 and through the subcritical regime. Once
 * convection develops, the interior mixes towards T ≈ ½ (where μ ≈ 1) and the
 * variation concentrates into thin boundary layers, so μ̄ overstates the radial
 * spread and the iteration count rises — which is the honest failure mode for a
 * preconditioner, and is what the fixed iteration budget is sized against.
 *
 * The power law is deliberately **not** represented here. It has no radial
 * profile to average — its variation is where the flow deforms, which is
 * azimuthal as much as radial — and the clamp keeps μ inside the same interval
 * the thermal term spans, so μ̄ remains as good a spectral match at n = 3 as it
 * is at n = 1 (measured).
 */
export const meanViscosity = (
  ri: number, ro: number, gamma: number,
): ((r: number) => number) => {
  const d = Math.log(ro / ri);
  return (r) => viscosity(Math.log(ro / r) / d, gamma);
};
