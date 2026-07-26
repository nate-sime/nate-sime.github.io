/** Temperature transport and coupled convection. */

import { describe, it, expect } from "vitest";
import { Temperature, type Velocity } from "../src/solver/temperature";
import { Simulation } from "../src/solver/step";

const RI = 0.55, RO = 1.0;

describe("temperature transport", () => {
  it("measures boundary heat flux to 4th order", () => {
    // Nu is the benchmark quantity, so its stencil must not dominate the error.
    const err = [33, 65, 129].map((nr) => {
      const T = new Temperature(nr, 16, RI, RO);
      for (let i = 0; i < nr; i++) T.T[i].fill(T.conduction(i));
      const nu = T.nusselt();
      return Math.max(Math.abs(nu.inner - 1), Math.abs(nu.outer - 1));
    });
    for (let i = 1; i < err.length; i++)
      expect(Math.log2(err[i - 1] / err[i])).toBeGreaterThan(3.5);
  });

  it("relaxes by diffusion to the analytic conduction profile", () => {
    const T = new Temperature(65, 32, RI, RO);
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
      const T = new Temperature(49, 96, RI, RO, 0, 0);
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
