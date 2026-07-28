/** Temperature transport and coupled convection. */

import { describe, it, expect } from "vitest";
import { ANNULUS, box } from "../src/geometry";
import { Temperature, type Velocity } from "../src/solver/temperature";
import { Simulation } from "../src/solver/step";

const RI = ANNULUS.lo, RO = ANNULUS.hi;

describe("temperature transport", () => {
  it("measures boundary heat flux to 4th order", () => {
    // Nu is the benchmark quantity, so its stencil must not dominate the error.
    const err = [33, 65, 129].map((nr) => {
      const T = new Temperature(ANNULUS, nr, 16);
      for (let i = 0; i < nr; i++) T.T[i].fill(T.conduction(i));
      const nu = T.nusselt();
      return Math.max(Math.abs(nu.inner - 1), Math.abs(nu.outer - 1));
    });
    for (let i = 1; i < err.length; i++)
      expect(Math.log2(err[i - 1] / err[i])).toBeGreaterThan(3.5);
  });

  it("relaxes by diffusion to the analytic conduction profile", () => {
    const T = new Temperature(ANNULUS, 65, 32);
    T.reset(0.3, 3);
    for (let n = 0; n < 400; n++) T.diffuse(2e-3);
    let e = 0;
    for (let i = 0; i < 65; i++)
      for (let j = 0; j < 32; j++) e = Math.max(e, Math.abs(T.T[i][j] - T.conduction(i)));
    expect(e).toBeLessThan(1e-5);
    expect(T.nusselt().outer).toBeCloseTo(1, 4);
  });

  it("advects a field around one rigid rotation, BFECC beating plain SL", () => {
    // u_φ = r ⇒ dφ/dt = 1, so the field must return to itself after t = 2π.
    const u: Velocity = (r) => ({ ur: 0, up: r });
    const run = (bfecc: boolean) => {
      const T = new Temperature(ANNULUS, 49, 96, 0, 0);
      for (let i = 0; i < 49; i++)
        for (let j = 0; j < 96; j++)
          T.T[i][j] = Math.sin((Math.PI * i) / 48) * Math.cos(2 * j * T.dphi);
      T.applyBC();
      const T0 = T.T.map((r) => Float64Array.from(r));
      for (let n = 0; n < 40; n++) T.advect(u, (2 * Math.PI) / 40, bfecc);
      let e = 0;
      for (let i = 0; i < 49; i++)
        for (let j = 0; j < 96; j++) e = Math.max(e, Math.abs(T.T[i][j] - T0[i][j]));
      return e;
    };
    const plain = run(false), bfecc = run(true);
    expect(bfecc).toBeLessThan(plain);
    expect(bfecc).toBeLessThan(1e-2);
  });
});

describe("coupled convection", () => {
  it("stays conductive when unforced (Ra = 0)", () => {
    const sim = new Simulation({ nr: 20, na: 32, gnr: 33, gna: 64, Ra: 0, dtMax: 2e-3 });
    sim.run(200, 1e-8);
    expect(sim.temp.nusselt().outer).toBeCloseTo(1, 3);
  });

  it("convects above critical Ra and conserves heat globally", () => {
    const sim = new Simulation({ nr: 20, na: 32, gnr: 33, gna: 64, Ra: 1e4, dtMax: 2e-3 });
    const res = sim.run(600, 1e-7);
    const nu = sim.temp.nusselt();
    expect(res.converged).toBe(true);
    expect(nu.outer).toBeGreaterThan(1.2);            // genuine convective transport
    expect(nu.inner).toBeCloseTo(nu.outer, 2);        // heat in = heat out at steady state
  });

  it("keeps the convecting velocity exactly divergence-free", () => {
    const sim = new Simulation({ nr: 20, na: 32, gnr: 33, gna: 64, Ra: 1e4, dtMax: 2e-3 });
    for (let n = 0; n < 40; n++) sim.step();
    let d = 0;
    for (let i = 1; i < 30; i++)
      for (let j = 0; j < 40; j++)
        d = Math.max(d, Math.abs(sim.psi.divergence(
          RI + ((RO - RI) * i) / 30, (2 * Math.PI * j) / 40)));
    expect(d).toBeLessThan(1e-12);
  });
});

/**
 * The Cartesian box.
 *
 * These are the same properties as above, restated on the second geometry — and
 * that is the point of them. The box is not a second solver: it is the same
 * discretisation with `h = 1` instead of `h = r` (see `src/geometry.ts`), so
 * anything that holds in one and not the other is a metric term in the wrong
 * place rather than a difference in physics. The critical-Rayleigh check at the
 * end is the one that would notice.
 */
describe("Cartesian box", () => {
  const grid = { nr: 20, na: 32, gnr: 33, gna: 64 } as const;

  it("measures boundary heat flux to 4th order", () => {
    const g = box(4);
    const err = [33, 65, 129].map((nr) => {
      const T = new Temperature(g, nr, 16);
      for (let i = 0; i < nr; i++) T.T[i].fill(T.conduction(i));
      const nu = T.nusselt();
      return Math.max(Math.abs(nu.inner - 1), Math.abs(nu.outer - 1));
    });
    // The conduction profile is *linear* here, so the one-sided stencil is exact
    // on it up to round-off — there is no O(dr⁴) trend to measure, and demanding
    // one would be measuring noise. What matters is that both ends read 1.
    for (const e of err) expect(e).toBeLessThan(1e-12);
  });

  it("relaxes by diffusion to the linear conduction profile", () => {
    const T = new Temperature(box(4), 65, 32);
    T.reset(0.3, 3);
    for (let n = 0; n < 400; n++) T.diffuse(2e-3);
    let e = 0;
    for (let i = 0; i < 65; i++)
      for (let j = 0; j < 32; j++) e = Math.max(e, Math.abs(T.T[i][j] - T.conduction(i)));
    expect(e).toBeLessThan(1e-5);
    expect(T.nusselt().outer).toBeCloseTo(1, 4);
  });

  it("stays conductive when unforced (Ra = 0)", () => {
    const sim = new Simulation({ geom: box(4), ...grid, Ra: 0, dtMax: 2e-3 });
    sim.run(200, 1e-8);
    expect(sim.temp.nusselt().outer).toBeCloseTo(1, 3);
  });

  it("convects above critical Ra and conserves heat globally", () => {
    const sim = new Simulation({ geom: box(4), ...grid, Ra: 1e4, dtMax: 2e-3 });
    const res = sim.run(600, 1e-7);
    const nu = sim.temp.nusselt();
    expect(res.converged).toBe(true);
    expect(nu.outer).toBeGreaterThan(1.2);
    expect(nu.inner).toBeCloseTo(nu.outer, 2);
  });

  it("keeps the convecting velocity exactly divergence-free", () => {
    const g = box(4);
    const sim = new Simulation({ geom: g, ...grid, Ra: 1e4, dtMax: 2e-3 });
    for (let n = 0; n < 40; n++) sim.step();
    let d = 0;
    for (let i = 1; i < 30; i++)
      for (let j = 0; j < 40; j++)
        d = Math.max(d, Math.abs(sim.psi.divergence(i / 30, (g.span * j) / 40)));
    expect(d).toBeLessThan(1e-12);
  });

  /**
   * **The critical Rayleigh number**, which is the sharpest single statement
   * available about whether the Cartesian metric is right.
   *
   * Linear theory for a free-slip layer of unit depth gives a growth rate
   * `σ = Ra k²/q² − q` with `q = π² + k²`, so a mode is marginal at
   * `Ra_c(k) = q³/k²`. That is minimised at `k = π/√2`, where
   * `Ra_c = 27π⁴/4 = 657.5`. A box of length `2√2` is exactly one wavelength of
   * that mode, so seeding a single roll into it puts the layer at the minimum
   * and its stability is decided by Ra alone — 400 must decay, 1200 must grow.
   *
   * Every metric term is implicated in that number. A stray `1/h` in the Stokes
   * operator, a wrong Jacobian in the assembly, a transverse wavenumber that
   * still assumed a 2π period in the diffusion: each moves `Ra_c`, and none of
   * them would be visible in a run that merely looked like convection. It is
   * also the check that the box is being *driven* correctly, since the buoyancy
   * load is the one assembly that needed no change at all.
   */
  it("puts the onset of convection at the analytic critical Rayleigh number", () => {
    // One wavelength of the k = π/√2 mode, so mode 1 is the critical roll.
    const excess = (Ra: number) => {
      const sim = new Simulation({
        geom: box(2 * Math.SQRT2), ...grid,
        Ra, dtMax: 2e-3, seed: { amp: 0.01, mode: 1 },
      });
      for (let n = 0; n < 400; n++) sim.step();
      return sim.temp.nusselt().outer - 1;
    };
    // Ra_c = 27π⁴/4 = 657.5 sits between these two.
    expect((27 * Math.PI ** 4) / 4).toBeCloseTo(657.5, 1);
    // Subcritical: σ < 0, the roll decays, and Nu − 1 is quadratic in what is
    // left of it — so this is a very long way below the supercritical value.
    expect(excess(400)).toBeLessThan(1e-3);
    // Supercritical: σ > 0 and the roll saturates into a finite-amplitude cell.
    expect(excess(1200)).toBeGreaterThan(0.02);
  });

  /**
   * **Blankenbach et al. (1989), case 1a**: a unit square with free-slip,
   * isothermal top and bottom at Ra = 10⁴, whose accepted Nusselt number is
   * 4.884409. It is the reference every Cartesian mantle convection code is
   * checked against, and the onset test above cannot substitute for it — that
   * one pins the *linear* problem, this one pins finite-amplitude heat transport
   * including the advection scheme and the boundary-flux stencil.
   *
   * The domain here is periodic and the benchmark's is walled, and that is not a
   * discrepancy to wave away: a periodic box of length 2 carrying **one full
   * wavelength** is a symmetric pair of counter-rotating unit cells, and the
   * symmetry planes at x = 0, 1, 2 satisfy u_x = 0 and ∂u_z/∂x = 0 — which is
   * exactly free-slip. So the half-cell *is* the benchmark domain, reached by the
   * only boundary condition a radix-2 transform can impose.
   *
   * The tolerance is 1%, against a measured 0.2% at this resolution: this is a
   * physics check, not a convergence study, and it is deliberately run on the
   * coarsest grid that resolves the boundary layers so the suite stays usable.
   */
  it("reproduces the Blankenbach 1a benchmark Nusselt number", () => {
    const sim = new Simulation({
      geom: box(2), nr: 16, na: 32, gnr: 33, gna: 64,
      Ra: 1e4, dtMax: 1e-3, seed: { amp: 0.1, mode: 1 },
    });
    const res = sim.run(2500, 1e-8);
    const nu = sim.temp.nusselt();
    expect(res.converged).toBe(true);
    expect(nu.outer).toBeCloseTo(4.884409, 1);
    expect(Math.abs(nu.outer / 4.884409 - 1)).toBeLessThan(0.01);
    // At a genuine steady state the two boundary fluxes are the same number
    // twice over, from independent reductions at opposite ends of the layer.
    expect(nu.inner).toBeCloseTo(nu.outer, 8);
  });
});
