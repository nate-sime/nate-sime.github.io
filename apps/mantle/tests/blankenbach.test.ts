/**
 * The Blankenbach et al. (1989) **variable-viscosity** cases, 2a and 2b.
 *
 * Case 1a — constant μ — is checked in `temperature.test.ts`, where it pins the
 * advection scheme and the boundary-flux stencil against a published number.
 * These two are that check made again with the rheology in the loop: 2a is
 * temperature dependence alone at a 10³ contrast, and 2b adds depth dependence
 * on top of a 1.6×10⁴ one — which makes it the one published number this code's
 * depth term can be held to, and 2a the one its thermal term can. That is why
 * they are worth minutes of the suite's runtime.
 *
 * **The translation to the centred law is exact, and it is one factor.** The
 * benchmark states μ = exp(−b T + c d) with the reference viscosity at the cold
 * surface (T = 0, d = 0); `rheology.ts` centres both exponents,
 * μ = exp(−γ(T−½) + c(d−½)). The two differ by the constant exp((γ−c)/2), and a
 * constant on μ is a rescaling of Ra and nothing else — the Stokes problem is
 * linear in μ, the load is linear in Ra, so μ → Kμ with Ra → K·Ra leaves ψ, u and
 * T *identical*. Running at `Ra·exp((γ−c)/2)` therefore does not approximate the
 * benchmark, it is the benchmark: ×31.6 for 2a, ×16 for 2b.
 *
 * **What the depth term is worth here is not subtle.** Take it out of 2b and the
 * reference viscosity moves by e^c = 64, so the effective Rayleigh number moves
 * by 64 — and reversing its sign moves it by 64 the other way. Neither is a few
 * percent from the benchmark; both are a different problem. Agreement to the
 * tolerances below is therefore a statement about the term itself, not about
 * being in the right ballpark.
 *
 * **The tolerances are resolution, and the resolution is a runtime decision.**
 * Case 1a reproduces to 0.2% on a 33×64 grid; these do not on theirs, and the
 * reason is physical rather than a defect. At Nu ≈ 10 the thermal boundary
 * layers are half as thick as case 1a's, and a viscosity contrast of 10³ (2a) or
 * 1.6×10⁴ (2b) makes one of them thinner and faster still — so the same grid
 * resolves it far less well, and the *surface flux* is the quantity that suffers,
 * being a one-sided gradient taken across it. Measured on the grids committed
 * here, marched well past the step counts the tests use, and on ψ/T grids 1.5×
 * finer in each direction:
 *
 *            committed grids                   1.5× finer
 *   2a   Nu  10.67  (+6.0%)  V 517  (+7.7%)    Nu  9.99  (−0.8%)  V 483 (+0.5%)
 *   2b   Nu   7.93 (+14.5%)  V 179  (+4.0%)    Nu ≤7.20  (+3.9%)  V 177 (+3.1%)
 *
 * — 2b's finer figure was still falling when it was stopped, so it bounds that
 * grid's error rather than stating it. Both converge on the published values, at
 * about the first order Nu converges at here. The refinement also closes the gap
 * between the two *boundary fluxes*, which is the internal evidence that this is
 * resolution and not a wrong answer: at a true steady state they are one number,
 * and 2b's disagree by 9% on the committed grid against 1.5% on the finer one.
 * The finer runs take four to ten minutes each, which is what the suite is not
 * paying; the committed pair costs about 4½ minutes, and the assertions are set
 * just outside what they measure.
 *
 * Both are run in the **walled** box the benchmark states — the mirrored domain
 * of `geometry.ts` — seeded with the single cell its reference setups state as
 * their initial condition, half a cosine across the width.
 */

import { describe, it, expect } from "vitest";
import { box } from "../src/geometry";
import { Simulation } from "../src/solver/step";
import { gammaFor } from "../src/solver/rheology";

/**
 * V_rms over the domain — the benchmark's second published quantity, and worth
 * having precisely because it is not the first: Nu is a temperature gradient at
 * one boundary, V_rms a volume integral of the velocity. A rheology that was
 * wrong in a way that happened to leave the surface flux alone would have to be
 * wrong compatibly in both.
 *
 * A midpoint average of |u|² over a tensor grid: spectrally accurate across the
 * period, second order in depth, and at 129 points either error is orders below
 * the ones being measured. In a walled box the sampled span is the mirrored
 * period, whose mean is the mean over the visible half.
 */
const vrms = (sim: Simulation, m = 128): number => {
  const g = sim.geom;
  let s = 0, n = 0;
  for (let i = 0; i <= m; i++) {
    const r = g.lo + ((g.hi - g.lo) * i) / m;
    for (let j = 0; j < m; j++) {
      const { ur, up } = sim.psi.velocity(r, (g.span * j) / m);
      s += ur * ur + up * up;
      n++;
    }
  }
  return Math.sqrt(s / n);
};

/**
 * Ra to solve at, so that the centred law is the benchmark's uncentred one. See
 * the header: this is the whole of the translation.
 */
const benchRa = (Ra: number, gamma: number, cz: number): number =>
  Ra * Math.exp((gamma - cz) / 2);

/**
 * A fixed step count rather than a convergence test on Nu. `Simulation.run`
 * stops when Nu stops moving *between polls*, and these cases approach their
 * steady state monotonically and slowly enough that it would stop during the
 * approach and report a number still on its way — a passing test measuring the
 * transient. The counts below are measured instead: another 1000 steps moves Nu
 * by 0.2% (2a) and 0.4% (2b), against tolerances of 10% and 18%.
 *
 * `cfl = 4` sizes the step. dt is an accuracy parameter here, not a stability
 * limit — semi-Lagrangian advection and implicit diffusion are both
 * unconditionally stable — so what it buys is runtime and what it costs is
 * accuracy, and both were measured rather than assumed on case 2a: halving it to
 * `cfl = 2` doubles the step count and moves the settled Nu by 0.7%, an order
 * below the resolution error being reported, while `cfl = 8` is *not*
 * interchangeable — it takes a visibly different trajectory and sits 7% lower at
 * the same simulated time, still drifting.
 */
const settle = (sim: Simulation, steps: number): Simulation => {
  for (let n = 0; n < steps; n++) sim.step();
  return sim;
};

describe("Blankenbach variable-viscosity benchmarks", () => {
  /**
   * **Case 2a**: unit square, free-slip, Ra = 10⁴, μ = exp(−b T) with
   * b = ln(10³). Accepted Nu = 10.0660, V_rms = 480.4334.
   *
   * This is the μ(T) tier's first check against a published number — everything
   * else about it is verified against the manufactured solution, against the
   * constant-μ tier, or against its own preconditioner, all of which are checks
   * that the code solves the equations it says it does. This one says those
   * equations are the ones the field agreed on.
   */
  it("reproduces case 2a: temperature-dependent viscosity", () => {
    const gamma = gammaFor(1e3);
    const sim = settle(new Simulation({
      geom: box(1, "free-slip"), nr: 16, na: 32, gnr: 33, gna: 64,
      Ra: benchRa(1e4, gamma, 0), gamma, variable: true, iters: 12,
      cfl: 4, dtMax: 1e-3, seed: { amp: 0.1, mode: 1 },
    }), 2000);

    // Measured 10.64 and 515 here — this grid's own answer, 6% and 7% above the
    // published pair, and converging on them under refinement (see the header).
    const nu = sim.temp.nusselt();
    expect(Math.abs(nu.outer / 10.0660 - 1)).toBeLessThan(0.10);
    expect(Math.abs(vrms(sim) / 480.4334 - 1)).toBeLessThan(0.10);
  }, 300_000);

  /**
   * **Case 2b**: a 2.5 × 1 box, free-slip, Ra = 10⁴, and the case this rheology
   * exists to be able to run — μ = exp(−b T + c d) with b = ln(16384) and
   * c = ln(64), so the viscosity spans four decades in temperature *and*
   * stiffens 64-fold with depth. Accepted Nu = 6.9299, V_rms = 171.755.
   *
   * Depth is measured from the cold surface, which is the opposite sense to the
   * solver's radial coordinate — see `depthAt`. That is the one place a sign
   * error would survive every unit test and still produce a plausible picture,
   * and it is what this number is sensitive to: with the sign reversed the
   * reference viscosity moves by 64 and so does the effective Ra.
   */
  it("reproduces case 2b: temperature- and depth-dependent viscosity", () => {
    const gamma = gammaFor(16384), cz = gammaFor(64);
    const sim = settle(new Simulation({
      geom: box(2.5, "free-slip"), nr: 16, na: 64, gnr: 33, gna: 128,
      Ra: benchRa(1e4, gamma, cz), gamma, cz, variable: true, iters: 12,
      cfl: 4, dtMax: 1e-3, seed: { amp: 0.1, mode: 1 },
    }), 1200);

    // Measured 7.97 and 180. The surface flux is this configuration's weakest
    // quantity — a one-sided gradient across a boundary layer under a lid 256
    // times stiffer than the floor — so it carries the loosest bound here, while
    // V_rms, an integral over the whole domain, is held to 8%. Both bounds are
    // far tighter than the *factor* of 64 that dropping the depth term, or
    // reversing its sign, would move this case by.
    const nu = sim.temp.nusselt();
    expect(Math.abs(nu.outer / 6.9299 - 1)).toBeLessThan(0.18);
    expect(Math.abs(vrms(sim) / 171.755 - 1)).toBeLessThan(0.08);
  }, 300_000);
});
