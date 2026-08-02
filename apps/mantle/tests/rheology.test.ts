/**
 * μ(T, d) mode. CPU reference.
 *
 * The ladder here is the same shape as the constant-μ suite's: check the
 * *operator* before the solve, and the limiting case before the general one. Two
 * rungs carry most of the weight.
 *
 *   **A(A⁻¹b) = b.** The matrix-free apply and the tier-1 direct solve are
 *   independent code paths — one gathers over quadrature points, the other
 *   inverts assembled radial blocks per azimuthal mode — so composing them is a
 *   genuine cross-check that they are the same operator. Nothing else in this
 *   file would catch a mis-paired symbol or a dropped curvature term.
 *
 *   **γ → 0 reproduces tier 1.** With μ ≡ 1 the preconditioner *is* the
 *   operator, so PCG must converge in one step to the direct solution. A failure
 *   localises to the Krylov wrapper rather than the physics.
 */

import { describe, it, expect } from "vitest";
import { ANNULUS, annulus, box } from "../src/geometry";
import { clampedAxis, periodicAxis, Field } from "../src/spline";
import { StokesSolver, VariableStokes, loadVector } from "../src/solver/stokes";
import { viscosityAt, strainRate, operatorTables } from "../src/solver/assembly";
import { radialBlocks, azimuthalSymbols, radialOperator } from "../src/solver/operators";
import {
  viscosity, gammaFor, meanViscosity, strainScale, depthAt, EPS_MIN,
  tackleyViscosity, tackleyLinear, tackleyPlastic, meanTackleyViscosity,
  TACKLEY_TRANSITION_DEPTH, TACKLEY_STRAIN_FLOOR,
} from "../src/solver/rheology";
import { mat, lu, solve } from "../src/linalg";
import { source, Ri, Ro } from "../src/mms";

const NA = 64;
const axes = (nr: number) => [clampedAxis(nr, Ri, Ro), periodicAxis(NA)] as const;
/**
 * The MMS domain [Ri, Ro] is its own choice, independent of the app's default
 * `ANNULUS` — so anything that needs a `Geometry` matching these axes (the
 * radial viscosity profile, `depthAt`) must be built from Ri/Ro directly
 * rather than borrowing `ANNULUS`, or the profile is evaluated outside the
 * domain it is meant to describe the moment the two diverge.
 */
const MMS_ANNULUS = annulus(Ri, Ro);

/** A φ-dependent temperature field, so μ genuinely couples the azimuthal modes. */
const Tfield = (r: number, phi: number) =>
  Math.log(Ro / r) / Math.log(Ro / Ri) + 0.2 * Math.sin(3 * phi);

/** A deterministic interior field — boundary rows zero, as the trial space wants. */
const field = (nr: number, f: (i: number, j: number) => number) => {
  const m = mat(nr, NA);
  for (let i = 1; i < nr - 1; i++)
    for (let j = 0; j < NA; j++) m[i][j] = f(i, j);
  return m;
};

describe("rheology", () => {
  it("reduces to μ ≡ 1 exactly as the activation vanishes", () => {
    for (const T of [0, 0.25, 0.5, 0.9, 1]) expect(viscosity(T, 0)).toBe(1);
  });

  it("realises the requested viscosity contrast", () => {
    for (const c of [10, 1e3, 1e5]) {
      const g = gammaFor(c);
      expect(viscosity(0, g) / viscosity(1, g)).toBeCloseTo(c, 6);
      // Centred on T = ½, so the geometric mean is 1 and raising the contrast
      // does not quietly rescale the effective Rayleigh number.
      expect(Math.sqrt(viscosity(0, g) * viscosity(1, g))).toBeCloseTo(1, 12);
    }
    // T outside [0,1] is clamped, not extrapolated.
    expect(viscosity(1 + 1e-3, gammaFor(1e4))).toBe(viscosity(1, gammaFor(1e4)));
  });
});

/**
 * The depth term. Two things have to be exact rather than close.
 *
 * **c = 0 is the μ(T) law bit for bit**, so every existing call site, test and
 * measured number above stays a statement about the same function — the depth
 * dependence is an addition to the law, not a re-parameterisation of it.
 *
 * **The depth coordinate is not the axis coordinate.** `r` runs from the hot
 * boundary to the cold one in both geometries, so a depth term written in `r`
 * would stiffen the wrong half of the layer and still look plausible on screen:
 * a viscous lid and a viscous floor both suppress convection. Pinning the ends
 * is what distinguishes them.
 */
describe("depth dependence", () => {
  const C = gammaFor(1e2);

  it("is the μ(T) law exactly at zero depth contrast", () => {
    for (const T of [0, 0.3, 0.5, 1])
      for (const d of [0, 0.25, 1])
        expect(viscosity(T, gammaFor(1e3), 1, 1, d, 0))
          .toBe(viscosity(T, gammaFor(1e3)));
    // …and μ ≡ 1 when neither term is switched on, whatever the depth.
    for (const d of [0, 0.5, 1]) expect(viscosity(0.7, 0, 1, 1, d, 0)).toBe(1);
  });

  it("realises the requested depth contrast, centred on d = ½", () => {
    for (const c of [10, 1e2, 1e3]) {
      const g = gammaFor(c);
      // Deep over shallow, at fixed T: the ratio the slider names.
      expect(viscosity(0.5, 0, 1, 1, 1, g) / viscosity(0.5, 0, 1, 1, 0, g))
        .toBeCloseTo(c, 6);
      // Centred, so the geometric mean over the layer is 1 and raising the depth
      // contrast does not quietly rescale the effective Rayleigh number — the
      // same property the T = ½ centring buys, for the same reason.
      expect(Math.sqrt(viscosity(0.5, 0, 1, 1, 1, g) * viscosity(0.5, 0, 1, 1, 0, g)))
        .toBeCloseTo(1, 12);
    }
    // Sign: positive c stiffens the *deep* interior.
    expect(viscosity(0.5, 0, 1, 1, 1, C)).toBeGreaterThan(viscosity(0.5, 0, 1, 1, 0, C));
  });

  /**
   * The two terms multiply, and the clamp is widened to exactly the interval
   * they jointly span — attained at the opposite corners (T = 0, d = 1) and
   * (T = 1, d = 0). So "total contrast" is the product of the two the sliders
   * ask for, which is what the power law is then confined to redistribute
   * within, and what the fixed Krylov budget is sized against.
   */
  it("spans the product of the two contrasts, and no more", () => {
    const g = gammaFor(1e3), hi = Math.sqrt(1e3 * 1e2);
    expect(viscosity(0, g, 1, 1, 1, C)).toBeCloseTo(hi, 6);
    expect(viscosity(1, g, 1, 1, 0, C)).toBeCloseTo(1 / hi, 9);
    // The power law cannot widen it, at any strain rate.
    for (const s of [1e-12, 1e-3, 1, 1e6])
      for (const [T, d] of [[0, 1], [0.5, 0.5], [1, 0]] as const) {
        const mu = viscosity(T, g, s, 3, d, C);
        expect(mu).toBeLessThanOrEqual(hi * (1 + 1e-12));
        expect(mu).toBeGreaterThanOrEqual((1 / hi) * (1 - 1e-12));
      }
  });

  // Depth is measured from the *cold* boundary in both geometries, which is the
  // opposite end from the axis origin in the box and the outer radius on the
  // annulus. Both directions are pinned; a sign error would swap them.
  it("measures depth from the cold boundary in both geometries", () => {
    expect(depthAt(ANNULUS, ANNULUS.hi)).toBeCloseTo(0, 12);
    expect(depthAt(ANNULUS, ANNULUS.lo)).toBeCloseTo(1, 12);
    const b = box(2);
    expect(depthAt(b, 1)).toBeCloseTo(0, 12);      // z = 1 is the cold ceiling
    expect(depthAt(b, 0)).toBeCloseTo(1, 12);      // z = 0 is the hot floor
    expect(depthAt(b, 0.25)).toBeCloseTo(0.75, 12);
  });

  /**
   * The preconditioner carries the depth term *exactly*, because μ̄ is a radial
   * profile and depth is a function of r alone. That is not true of the thermal
   * term (exact only in the conductive state) or of the power law (not radial at
   * all), and it is the reason the depth slider does not cost Krylov iterations.
   */
  it("is represented exactly by the μ̄(r) profile", () => {
    for (const g of [ANNULUS, box(2)]) {
      const mean = meanViscosity(g, 0, C);
      for (const f of [0, 0.3, 0.5, 1]) {
        const r = g.lo + (g.hi - g.lo) * f;
        expect(mean(r)).toBeCloseTo(viscosity(0.5, 0, 1, 1, depthAt(g, r), C), 12);
      }
      // Stiff at depth, weak near the surface — the profile, not just its ends.
      expect(mean(g.lo)).toBeGreaterThan(mean(g.hi));
    }
  });
});

/**
 * The regularised power law, and its reduction to the linear tier as n → 1. That
 * reduction is asserted as an *identity* rather than a limit, which is what lets
 * both laws run through one set of kernels with `n` as a uniform: if it were
 * only approximate, the μ(T) tier would silently become a second numerical path.
 */
describe("power law", () => {
  const G = gammaFor(1e3);

  it("is the μ(T) law exactly at n = 1", () => {
    for (const T of [0, 0.3, 0.5, 1])
      for (const s of [1e-3, 0.5, 1, 7, 1e3])
        expect(viscosity(T, G, s, 1)).toBe(viscosity(T, G));
    // …and μ ≡ 1 at zero activation stays exact for every n, the isoviscous
    // check carried forward rather than weakened.
    for (const n of [1, 2, 3, 5]) expect(viscosity(0.2, 0, 4, n)).toBe(1);
  });

  it("realises the power-law exponent between the clamps", () => {
    // A mild contrast and strains near 1, so neither clamp is active and what is
    // measured is the exponent itself.
    const g = gammaFor(4);
    for (const n of [2, 3, 5]) {
      const ratio = viscosity(0.5, g, 1.3, n) / viscosity(0.5, g, 0.7, n);
      expect(ratio).toBeCloseTo((1.3 / 0.7) ** ((1 - n) / n), 12);
    }
    // Shear-thinning: n > 1 weakens what deforms fastest.
    expect(viscosity(0.5, g, 4, 3)).toBeLessThan(viscosity(0.5, g, 0.25, 3));
  });

  it("never exceeds the requested contrast", () => {
    for (const c of [10, 1e3, 1e5]) {
      const g = gammaFor(c), hi = Math.sqrt(c);
      for (const T of [0, 0.5, 1])
        for (const s of [1e-12, 1e-3, 1, 1e6]) {
          const mu = viscosity(T, g, s, 3);
          expect(mu).toBeLessThanOrEqual(hi * (1 + 1e-12));
          expect(mu).toBeGreaterThanOrEqual((1 / hi) * (1 - 1e-12));
        }
    }
  });

  /**
   * The centring, and the reason the reference is a *geometric* mean.
   *
   * The exponent is negative, so μ is convex in ε̇ and normalising by the RMS —
   * the obvious choice — stiffens almost everywhere, since ε̇ is peaked and most
   * of the domain sits below its RMS. Measured in the coupled loop that moved
   * ⟨log μ⟩ from 0.21 to 0.79 and cut max |u| from 27 to 6.6: the power law
   * switching convection off, which is exactly the silent rescaling of the
   * effective Ra that the T = ½ centring exists to prevent. Dividing by the
   * geometric mean makes ⟨log ŝ⟩ zero *by construction*, which is what this pins.
   */
  it("centres the normalised strain rate on a geometric mean of 1", () => {
    const e = Float64Array.from({ length: 4096 }, (_, i) =>
      Math.abs(Math.sin(0.7 * i) * Math.exp(2 * Math.cos(0.11 * i))));
    const { d, g } = strainScale(e);
    let l = 0;
    for (const v of e) l += Math.log((v + d) / g);
    expect(l / e.length).toBeCloseTo(0, 12);
    // The offset is relative, so the whole thing is scale-invariant in ψ: this
    // is what lets Ra move over four decades without the rheology drifting.
    const scaled = strainScale(e.map((v) => 1e4 * v));
    expect(scaled.g / g).toBeCloseTo(1e4, 6);
    expect(scaled.d / d).toBeCloseTo(1e4, 6);
    // A field with no scale at all — the state before the first solve.
    expect(strainScale(new Float64Array(64))).toEqual({ d: 1, g: 1 });
    expect(EPS_MIN).toBeGreaterThan(0);
  });

  // The tabulated strain rate rides on the operator's own eval tables, so a
  // mis-paired weight would show up as a wrong μ rather than as a wrong
  // operator — invisible to every μ(T) test. This checks it against the
  // independent definition in `Field`.
  it("computes the same strain rate as the field's second invariant", () => {
    const nr = 20, [rAx, aAx] = axes(nr);
    const f = new Field(rAx, aAx);
    f.interpolate((r, phi) =>
      Math.sin(4 * (r - Ri)) * Math.cos(3 * phi) + 0.2 * r * r);
    const t = operatorTables(rAx, aAx);
    const e = strainRate(t, f.c);

    let err = 0, sc = 0;
    for (let qr = 0; qr < t.rx.length; qr += 7)
      for (let qa = 0; qa < t.ax.length; qa += 13) {
        const ref = f.strainRate(t.rx[qr], t.ax[qa]);
        err = Math.max(err, Math.abs(e[qr * t.ax.length + qa] - ref));
        sc = Math.max(sc, Math.abs(ref));
      }
    expect(sc).toBeGreaterThan(0);
    expect(err / sc).toBeLessThan(1e-12);
  });
});

/**
 * The nonlinear solve. μ depends on ψ, so the linear PCG of the μ(T) tier is
 * wrapped in a time-lagged (Picard) update: solve with μ frozen, recompute μ
 * from the ψ that produced, repeat. What has to be true is that the lagging
 * *contracts* — otherwise the frame loop's single sweep approximates nothing.
 */
describe("nonlinear solve", () => {
  const nr = 24, G = gammaFor(1e3), N = 3;

  /** μ from ψ under the lagged law — the shared half of a sweep. */
  const lag = (vs: VariableStokes, x: Float64Array[], gamma: number) => {
    const e = strainRate(vs.tables, x);
    const { d, g } = strainScale(e);
    return viscosityAt(vs.tables, Tfield,
      (t, s) => viscosity(t, gamma, (s + d) / g, N), e);
  };

  /** One time-lagged rheology update and solve, exactly as `step.ts` does it. */
  const sweep = (
    vs: VariableStokes, load: Float64Array[], x: Float64Array[],
    iters: number, gamma = G,
  ) => vs.solve(load, lag(vs, x, gamma), x, iters);

  it("contracts under Picard sweeps", () => {
    const [rAx, aAx] = axes(nr);
    const load = loadVector(rAx, aAx, source);
    const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, G));

    const x = mat(nr, NA);
    const res: number[] = [];
    for (let k = 0; k < 6; k++) {
      sweep(vs, load, x, 60);
      // ‖b − A(μ[ψ])ψ‖ with μ *rebuilt from* ψ: the nonlinear residual, not the
      // linear one the solve just drove down against a frozen μ.
      res.push(vs.residual(load, lag(vs, x, G), x));
    }
    // The first sweep starts from ψ = 0, where the law has no scale to normalise
    // against, so the transient is not the contraction rate. Measure it on the
    // tail, where it is ~0.5 per sweep.
    expect(res[5]).toBeLessThan(0.5 * res[3]);
    expect(res[3]).toBeLessThan(0.5 * res[1]);
  });

  // Structural, and it must stay structural under a rheology that depends on the
  // solution itself: u = ∇×ψ is divergence-free for *any* coefficients.
  it("keeps the nonlinear velocity divergence-free", () => {
    const [rAx, aAx] = axes(nr);
    const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, G));
    const x = mat(nr, NA);
    for (let k = 0; k < 3; k++) sweep(vs, loadVector(rAx, aAx, source), x, 40);

    const f = new Field(rAx, aAx);
    for (let i = 0; i < nr; i++) f.c[i].set(x[i]);
    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++) {
      const r = Ri + ((Ro - Ri) * i) / 30;
      for (let j = 0; j < 40; j++) {
        const phi = (2 * Math.PI * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(r, phi)));
        const v = f.velocity(r, phi);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    }
    expect(speed).toBeGreaterThan(0);
    expect(div / speed).toBeLessThan(1e-12);
  });

  /**
   * The dissipation form carries the curved-boundary free-slip condition
   * automatically, and that extends to variable μ by moving μ inside the
   * integral — no re-derivation. The μ(T) tier confirms it at the
   * `h^{p−1} = h²` rate ψ'' dictates. This is the same claim where μ depends on
   * ψ's own second derivatives.
   *
   * **The rate does not survive; the condition does.** Measured at a 10²
   * contrast with the Picard iteration converged and every solve driven to a
   * 1e-12 residual, the wall shear rate falls 8.2e-1 → 3.5e-1 → 2.1e-1 over
   * 16×32 → 24×48 → 32×64 and then *rises* to 2.6e-1 at 48×96 — roughly first
   * order and then non-monotone, against a clean second order at n = 1. The
   * cause is the coefficient's regularity, not the boundary condition: the
   * viscosity clamp puts a kink in μ along a level set of ε̇ whose position is
   * itself part of the discrete solution, so μ is Lipschitz rather than smooth
   * and its kink locus moves with the mesh. Ruled out as explanations, each by
   * measurement: an unconverged Krylov solve (identical answers at 200, 800 and
   * 2500 iterations), an unconverged Picard iteration (identical at 8 and 20
   * sweeps), and the width of the stiff wall layer the regularisation creates
   * (unchanged across three decades of ε̇_min).
   *
   * So what is asserted is the statement that distinguishes free-slip from the
   * wrong `ω = 0` condition — that the wall traction is *small against the
   * interior* and shrinks under refinement — rather than a rate the model does
   * not support. Under `ω = 0` the wall value would be `2u_φ/r`, comparable to
   * the interior shear rate, and this would fail at every resolution.
   */
  it("still satisfies free-slip against a nonlinear μ", () => {
    const g2 = gammaFor(1e2);
    const ratio = ([[16, 32], [32, 64]] as const).map(([n, na]) => {
      const rAx = clampedAxis(n, Ri, Ro), aAx = periodicAxis(na);
      const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, g2));
      const load = loadVector(rAx, aAx, source);
      const x = mat(n, na);
      for (let k = 0; k < 6; k++) sweep(vs, load, x, 80, g2);

      const f = new Field(rAx, aAx);
      for (let i = 0; i < n; i++) f.c[i].set(x[i]);
      let wall = 0, interior = 0;
      for (let j = 0; j < 64; j++) {
        const phi = (2 * Math.PI * j) / 64;
        wall = Math.max(wall,
          Math.abs(f.shearRate(Ri, phi)), Math.abs(f.shearRate(Ro, phi)));
        for (let i = 1; i < 40; i++)
          interior = Math.max(interior,
            Math.abs(f.shearRate(Ri + ((Ro - Ri) * i) / 40, phi)));
      }
      expect(interior).toBeGreaterThan(0);
      return wall / interior;
    });
    expect(ratio[1]).toBeLessThan(0.05);
    expect(ratio[1]).toBeLessThan(0.5 * ratio[0]);
  });
});

describe("variable-μ operator", () => {
  // In the unconstrained space the k = 0 block is rank-deficient by two. It is
  // not, once ψ = 0 is imposed at *both* radii: rigid rotation ψ = −Ωr²/2 cannot vanish at
  // two distinct radii, so the essential condition removes the whole kernel and
  // no angular-momentum side constraint is needed. This pins that.
  it("has a nonsingular k = 0 radial block", () => {
    for (const nr of [16, 32]) {
      const [rAx, aAx] = axes(nr);
      const R = radialBlocks(rAx), S = azimuthalSymbols(aAx);
      const A = radialOperator(R, S, 0), n = A.length;
      const x = Float64Array.from({ length: n }, (_, i) => Math.sin(3 * i));
      const b = new Float64Array(n);
      for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++) b[i] += A[i][j] * x[j];
      const y = solve(lu(radialOperator(R, S, 0)), b);
      let e = 0, sc = 0;
      for (let i = 0; i < n; i++) {
        e = Math.max(e, Math.abs(y[i] - x[i])); sc = Math.max(sc, Math.abs(x[i]));
      }
      expect(e / sc).toBeLessThan(1e-9);
    }
  });

  // The decisive one: the matrix-free gather must be the operator the assembled
  // per-mode blocks invert, or every downstream number is measuring two
  // different problems agreeing with themselves.
  it("inverts to the tier-1 direct solve for constant μ", () => {
    const nr = 24, [rAx, aAx] = axes(nr);
    const load = loadVector(rAx, aAx, source);
    const vs = new VariableStokes(rAx, aAx, () => 1);
    const mu = viscosityAt(vs.tables, () => 0.5, () => 1);

    const psi = new StokesSolver(rAx, aAx, () => 1, true).solve(load);
    psi[0].fill(0); psi[nr - 1].fill(0);          // the trial space
    const back = vs.apply(psi, mu);

    let e = 0, sc = 0;
    for (let i = 1; i < nr - 1; i++)
      for (let j = 0; j < NA; j++) {
        e = Math.max(e, Math.abs(back[i][j] - load[i][j]));
        sc = Math.max(sc, Math.abs(load[i][j]));
      }
    expect(e / sc).toBeLessThan(1e-11);
  });

  // Symmetry is what makes CG legal at all. It is not automatic: the boundary
  // rows have to be eliminated on both sides, and the C[ψ]C[v] cross terms have
  // to pair N with N'' the same way in each argument.
  it("is symmetric under a strongly varying viscosity", () => {
    const nr = 20, [rAx, aAx] = axes(nr);
    const vs = new VariableStokes(rAx, aAx, () => 1);
    const mu = viscosityAt(vs.tables, Tfield, (t) => viscosity(t, gammaFor(1e4)));
    const u = field(nr, (i, j) => Math.sin(11 * i + 7 * j) * 0.5);
    const v = field(nr, (i, j) => Math.cos(3 * i - 5 * j) + 0.2 * i);

    const Au = vs.apply(u, mu), Av = vs.apply(v, mu);
    let a = 0, b = 0;
    for (let i = 0; i < nr; i++)
      for (let j = 0; j < NA; j++) { a += Au[i][j] * v[i][j]; b += u[i][j] * Av[i][j]; }
    expect(Math.abs(a - b) / Math.abs(a)).toBeLessThan(1e-12);
  });
});

describe("variable-μ solve", () => {
  it("reproduces the tier-1 solution as the activation vanishes", () => {
    const nr = 24, [rAx, aAx] = axes(nr);
    const load = loadVector(rAx, aAx, source);
    const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, 0));
    const mu = viscosityAt(vs.tables, Tfield, (t) => viscosity(t, 0));

    // A budget of *one*: μ ≡ 1 makes the preconditioner the exact inverse of the
    // operator, so PCG must land on the answer in a single step. Asking for one
    // iteration is what turns this from "the answer eventually agrees" into a
    // statement that the two tiers are the same linear problem.
    const x = mat(nr, NA);
    const b = vs.residual(load, mu, x);   // x = 0, so this is ‖b‖
    vs.solve(load, mu, x, 1);
    expect(vs.residual(load, mu, x) / b).toBeLessThan(1e-12);

    const ref = new StokesSolver(rAx, aAx).solve(load);
    let e = 0, sc = 0;
    for (let i = 0; i < nr; i++)
      for (let j = 0; j < NA; j++) {
        e = Math.max(e, Math.abs(x[i][j] - ref[i][j]));
        sc = Math.max(sc, Math.abs(ref[i][j]));
      }
    expect(e / sc).toBeLessThan(1e-10);
  });

  /**
   * A **purely depth-dependent** μ is still a radial profile, so the tier-1 DFT
   * solve is not merely a good preconditioner for it — it is the operator. The
   * modes do not couple, the per-mode blocks carry μ(r) inside their own
   * quadrature, and PCG must therefore land on the direct solution in a single
   * iteration at any depth contrast.
   *
   * That is the sharpest statement available about where the depth term is
   * carried, and it is one the thermal term cannot make: μ(T) is only radial in
   * the conductive state, so its preconditioner is an approximation from the
   * first convecting step onwards. It is also the reason the depth slider does
   * not cost iterations — the budget sized for γ is untouched by c.
   */
  it("is preconditioned exactly by μ̄(r) when only depth varies", () => {
    const nr = 24, [rAx, aAx] = axes(nr);
    const c = gammaFor(1e3);
    const load = loadVector(rAx, aAx, source);
    const mean = meanViscosity(MMS_ANNULUS, 0, c);
    const vs = new VariableStokes(rAx, aAx, mean);
    // γ = 0, so T drops out and μ is the depth profile alone — evaluated at the
    // quadrature points by the operator, and integrated into the radial blocks
    // by the preconditioner. The two must be the same operator.
    const mu = viscosityAt(vs.tables, Tfield,
      (t, _s, r) => viscosity(t, 0, 1, 1, depthAt(MMS_ANNULUS, r), c));

    const x = mat(nr, NA);
    const b = vs.residual(load, mu, x);   // x = 0, so this is ‖b‖
    vs.solve(load, mu, x, 1);
    expect(vs.residual(load, mu, x) / b).toBeLessThan(1e-12);

    const ref = new StokesSolver(rAx, aAx, mean, true).solve(load);
    let e = 0, sc = 0;
    for (let i = 0; i < nr; i++)
      for (let j = 0; j < NA; j++) {
        e = Math.max(e, Math.abs(x[i][j] - ref[i][j]));
        sc = Math.max(sc, Math.abs(ref[i][j]));
      }
    expect(sc).toBeGreaterThan(0);
    expect(e / sc).toBeLessThan(1e-10);
    // …and the profile really was varying: a stiff floor slows the flow it is
    // solved with, so this is not one step to the *isoviscous* answer.
    const iso = new StokesSolver(rAx, aAx, () => 1, true).solve(load);
    let d = 0;
    for (let i = 0; i < nr; i++)
      for (let j = 0; j < NA; j++) d = Math.max(d, Math.abs(iso[i][j] - ref[i][j]));
    expect(d / sc).toBeGreaterThan(0.1);
  });

  // The preconditioner is the point of the tier: without it the biharmonic's
  // κ ~ h⁻⁴ would make CG useless. Asserting a *rate* rather than a final
  // residual keeps this meaningful if the mesh or the contrast changes.
  it("converges under the μ̄(r) preconditioner at a 10³ contrast", () => {
    const nr = 24, [rAx, aAx] = axes(nr);
    const g = gammaFor(1e3);
    const load = loadVector(rAx, aAx, source);
    const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, g));
    const mu = viscosityAt(vs.tables, Tfield, (t) => viscosity(t, g));

    const b0 = vs.residual(load, mu, mat(nr, NA));   // ‖b‖, the cold residual
    const drop = (n: number) => {
      const x = mat(nr, NA);
      vs.solve(load, mu, x, n);
      return vs.residual(load, mu, x) / b0;
    };
    expect(drop(8)).toBeLessThan(0.2);
    expect(drop(16)).toBeLessThan(1e-2);
    expect(drop(32)).toBeLessThan(1e-5);

    // Warm-started — which is what the frame loop does — a budget that would not
    // converge from cold still *holds* the residual, which is the property the
    // fixed per-frame budget actually relies on.
    const x = mat(nr, NA);
    vs.solve(load, mu, x, 32);
    vs.solve(load, mu, x, 12);
    expect(vs.residual(load, mu, x) / b0).toBeLessThan(1e-8);
  });

  // Structural, and it must stay structural: u = ∇×ψ is divergence-free for
  // *any* coefficients, so the rheology cannot degrade it. This is the property
  // the whole formulation exists to buy.
  it("keeps the variable-μ velocity divergence-free", () => {
    const nr = 24, [rAx, aAx] = axes(nr);
    const g = gammaFor(1e3);
    const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, g));
    const mu = viscosityAt(vs.tables, Tfield, (t) => viscosity(t, g));
    const x = mat(nr, NA);
    vs.solve(loadVector(rAx, aAx, source), mu, x, 24);

    const f = new Field(rAx, aAx);
    for (let i = 0; i < nr; i++) f.c[i].set(x[i]);
    let div = 0, speed = 0;
    for (let i = 1; i < 30; i++) {
      const r = Ri + ((Ro - Ri) * i) / 30;
      for (let j = 0; j < 40; j++) {
        const phi = (2 * Math.PI * j) / 40;
        div = Math.max(div, Math.abs(f.divergence(r, phi)));
        const v = f.velocity(r, phi);
        speed = Math.max(speed, Math.abs(v.ur), Math.abs(v.up));
      }
    }
    expect(speed).toBeGreaterThan(0);
    expect(div / speed).toBeLessThan(1e-12);
  });

  /**
   * The dissipation form was chosen because its natural boundary condition
   * is σ_rφ = 0 *with the curvature term carried automatically* — which carries
   * over to variable μ by moving μ inside the integral, with no
   * re-derivation. This is that claim under test: ε_rφ must still go to zero at
   * both radii, at the h^(p−1) = h² rate ψ'' dictates. An ω = 0 implementation,
   * or a μ that had been pulled outside the derivatives, would plateau instead.
   */
  it("still satisfies free-slip, ε_rφ → 0 at O(h²), with variable μ", () => {
    const g = gammaFor(1e2);
    const err = [16, 32, 64].map((nr) => {
      const [rAx, aAx] = axes(nr);
      const vs = new VariableStokes(rAx, aAx, meanViscosity(MMS_ANNULUS, g));
      const mu = viscosityAt(vs.tables, Tfield, (t) => viscosity(t, g));
      const x = mat(nr, NA);
      vs.solve(loadVector(rAx, aAx, source), mu, x, 60);

      const f = new Field(rAx, aAx);
      for (let i = 0; i < nr; i++) f.c[i].set(x[i]);
      let e = 0;
      for (let j = 0; j < NA; j++) {
        const phi = (2 * Math.PI * j) / NA;
        e = Math.max(e, Math.abs(f.shearRate(Ri, phi)), Math.abs(f.shearRate(Ro, phi)));
      }
      return e;
    });
    for (let i = 1; i < err.length; i++)
      expect(Math.log2(err[i - 1] / err[i])).toBeGreaterThan(1.8);
  });
});

/**
 * Tackley (2000): an Arrhenius branch and a Bingham (yield-stress) branch
 * combined as resistors in parallel, so the weaker sets η.
 */
describe("Tackley viscosity", () => {
  const T = 0.5, sigmaY = 1, sigmaB = 1, etaStar = 1e-3;

  it("steps A0 30x at the 670 km discontinuity", () => {
    const above = tackleyLinear(T, TACKLEY_TRANSITION_DEPTH - 1e-6);
    const below = tackleyLinear(T, TACKLEY_TRANSITION_DEPTH + 1e-6);
    expect(below / above).toBeCloseTo(30, 6);
  });

  it("clamps T to [0, 1], not extrapolated", () => {
    expect(tackleyLinear(1 + 1e-3, 0)).toBe(tackleyLinear(1, 0));
    expect(tackleyLinear(-1e-3, 0)).toBe(tackleyLinear(0, 0));
  });

  it("is the parallel combination of its two branches, so it never exceeds either", () => {
    for (const strain of [1e-6, 1e-2, 1, 100])
      for (const d of [0, 0.5, 1]) {
        const lin = tackleyLinear(T, d);
        const plast = tackleyPlastic(d, strain, sigmaY, sigmaB, etaStar);
        const eta = tackleyViscosity(T, d, strain, sigmaY, sigmaB, etaStar);
        expect(eta).toBeLessThanOrEqual(Math.min(lin, plast) * (1 + 1e-12));
        expect(eta).toBeCloseTo(1 / (1 / lin + 1 / plast), 12);
      }
  });

  // The floor rescues the division from ε̇ = 0; it is not part of the physics,
  // which already has a finite limit there (η_plast → ∞, η → η_lin).
  it("reduces to η_lin as ε̇ → 0", () => {
    const lin = tackleyLinear(T, 0.5);
    const eta = tackleyViscosity(T, 0.5, TACKLEY_STRAIN_FLOOR, sigmaY, sigmaB, etaStar);
    expect(eta / lin).toBeCloseTo(1, 3);
  });

  it("reduces η_plast to η* as ε̇ → ∞", () => {
    expect(tackleyPlastic(1, 1e12, sigmaY, sigmaB, etaStar) / etaStar).toBeCloseTo(1, 6);
  });

  it("shear-thins: η is non-increasing in ε̇", () => {
    let prev = Infinity;
    for (const s of [1e-3, 1e-1, 1, 10, 1e3]) {
      const eta = tackleyViscosity(T, 0.5, s, sigmaY, sigmaB, etaStar);
      expect(eta).toBeLessThanOrEqual(prev);
      prev = eta;
    }
  });

  it("is a pure floor when the yield stress is switched off", () => {
    for (const strain of [1e-3, 1, 1e3])
      expect(tackleyPlastic(0.7, strain, 0, 0, etaStar)).toBe(etaStar);
  });

  // μ̄(r) does not take σ_Y, σ_b or η* as arguments at all — the ε̇ → 0 limit
  // that makes it exact does not depend on them, so there is nothing to check
  // them against here.
  it("is η_lin evaluated on the conduction profile, exactly", () => {
    for (const g of [ANNULUS, box(2)]) {
      const mean = meanTackleyViscosity(g);
      for (const f of [0, 0.3, 0.5, 1]) {
        const r = g.lo + (g.hi - g.lo) * f;
        expect(mean(r)).toBe(tackleyLinear(g.conduction(r), depthAt(g, r)));
      }
    }
  });
});
