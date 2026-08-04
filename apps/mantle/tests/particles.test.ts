/**
 * The f64 tracer cloud: the marker-in-cell projection on its own, and the
 * thermochemical coupling it feeds — `Ra·T − Rb·C` in the buoyancy load. The RK2
 * push itself (the pathline machinery) is exercised implicitly by every test
 * below that steps a coupled `Simulation`; what these tests are really after
 * is the two things §5 of the plan adds on top of it: does the projected
 * composition converge to the field it was drawn from, and does turning the
 * coupling on (or leaving it off) do exactly what it says.
 */

import { describe, it, expect } from "vitest";
import { Simulation } from "../src/solver/step";
import { Particles } from "../src/solver/particles";
import { ANNULUS, box } from "../src/geometry";
import { depthOf, SPECIES_CONDITIONS } from "../src/particles";

describe("composition projection", () => {
  // A smooth field, not the two-species step function the default initial
  // condition paints on: the step function's own discontinuity would make
  // "how close does the marker-in-cell estimate come to the field it was
  // drawn from" a question about Gibbs ringing rather than about the
  // projection, so this overwrites `c` with something the bicubic grid
  // representation can actually resolve exactly at infinite N.
  const analytic = (geom: typeof ANNULUS, r: number, phi: number) =>
    0.5 + 0.3 * Math.cos(phi) * Math.sin(Math.PI * depthOf(geom, r));

  it("recovers a smooth composition field within its bias-plus-noise budget", () => {
    const p = new Particles(ANNULUS, 65, 128, { seed: 11 });
    for (let i = 0; i < p.r.length; i++) p.c[i] = analytic(ANNULUS, p.r[i], p.phi[i]);
    p.project();

    let sum = 0, n = 0, worst = 0;
    for (let i = 0; i < p.cnr; i++) {
      const r = ANNULUS.lo + i * p.cdr;
      for (let j = 0; j < p.cna; j++) {
        const phi = j * p.cdphi;
        const err = Math.abs(p.C[i][j] - analytic(ANNULUS, r, phi));
        sum += err; n++; worst = Math.max(worst, err);
      }
    }
    // TARGET_PARTICLES_PER_CELL = 16 sets the noise floor; the mean absolute
    // error averages that down across ~8000 nodes, so it is the tighter and
    // more reliable of the two bounds. `worst` only guards against something
    // structurally broken (wrong grid, wrong weights), not against the one
    // unlucky cell a tracer method is allowed to have.
    expect(sum / n).toBeLessThan(0.05);
    expect(worst).toBeLessThan(0.6);
  });

  it("a sparser cloud measurably raises the noise floor", () => {
    const rms = (count: number) => {
      const p = new Particles(ANNULUS, 65, 128, { seed: 21, count });
      for (let i = 0; i < p.r.length; i++) p.c[i] = analytic(ANNULUS, p.r[i], p.phi[i]);
      p.project();
      let sq = 0, n = 0;
      for (let i = 0; i < p.cnr; i++) {
        const r = ANNULUS.lo + i * p.cdr;
        for (let j = 0; j < p.cna; j++) {
          const phi = j * p.cdphi;
          const e = p.C[i][j] - analytic(ANNULUS, r, phi);
          sq += e * e; n++;
        }
      }
      return Math.sqrt(sq / n);
    };
    const cnr = Math.floor((65 + 1) / 2), cna = 128 / 2, cells = cnr * cna;
    // A cloud two orders below the design target (16/cell): the bias this
    // smooth, low-wavenumber field leaves on a 33×64 grid is small enough
    // that it, not the tracer noise, dominates near the design density (the
    // regression above measures that regime) — so this checks the noise
    // term the other way, by starving it until it is forced to dominate.
    const starved = rms(Math.floor(0.1 * cells));
    const design = rms(16 * cells);
    expect(starved).toBeGreaterThan(3 * design);
  });
});

describe("thermochemical coupling", () => {
  const opts = { nr: 12, na: 24, gnr: 17, gna: 32, geom: ANNULUS, Ra: 1e4 } as const;

  it("defaults Rb to 0 — particles present but not felt by the flow", () => {
    const sim = new Simulation({ ...opts, particles: { seed: 3 } });
    expect(sim.Rb).toBe(0);
  });

  it("Rb = 0 reproduces a run with no particles at all, exactly", () => {
    const plain = new Simulation(opts);
    const withParticles = new Simulation({ ...opts, particles: { seed: 3, layerDepth: 0.3 } });
    for (let n = 0; n < 8; n++) { plain.step(1e-4); withParticles.step(1e-4); }
    for (let i = 0; i < plain.psi.c.length; i++)
      for (let j = 0; j < plain.psi.c[i].length; j++)
        expect(withParticles.psi.c[i][j]).toBe(plain.psi.c[i][j]);
  });

  // The standard thermochemical setup: a dense layer sitting on the hot
  // (basal) boundary, exactly where the thermal instability itself sits, so
  // the two buoyancy sources are co-located and can be made to cancel or
  // reverse rather than merely superpose. A large enough buoyancy ratio must
  // visibly stabilise the layer against the thermal drive that would
  // otherwise lift it, which is the entrainment-vs-stable-blanket question
  // §5 of the plan exists to let particles answer.
  it("a dense basal layer damps convective vigour as the buoyancy ratio increases", () => {
    const seeded = { ...opts, seed: { amp: 0.1, mode: 4 } };
    const uncoupled = new Simulation({
      ...seeded, particles: { seed: 5, layerDepth: 0.3 },
    });
    const coupled = new Simulation({
      ...seeded, particles: { seed: 5, layerDepth: 0.3 }, Rb: 4,
    });
    for (let n = 0; n < 300; n++) { uncoupled.step(); coupled.step(); }
    const vUncoupled = uncoupled.temp.rmsVelocity(uncoupled.velocity);
    const vCoupled = coupled.temp.rmsVelocity(coupled.velocity);
    expect(vUncoupled).toBeGreaterThan(0);
    expect(vCoupled).toBeLessThan(vUncoupled);
  });

  it("is mutable at runtime, like the GPU twin's uniform", () => {
    const sim = new Simulation({ ...opts, particles: { seed: 3, layerDepth: 0.3 } });
    for (let n = 0; n < 5; n++) sim.step();
    const before = sim.psi.c.map((r) => Float64Array.from(r));
    sim.Rb = 2;
    sim.step();
    let moved = 0;
    for (let i = 0; i < sim.psi.c.length; i++)
      for (let j = 0; j < sim.psi.c[i].length; j++)
        moved = Math.max(moved, Math.abs(sim.psi.c[i][j] - before[i][j]));
    expect(moved).toBeGreaterThan(0);
  });
});

/**
 * The purely compositional limit van Keken et al. (1997)'s Rayleigh–Taylor
 * case needs: Ra = 0, so the thermal term drops out of the load entirely and
 * only Rb is left to drive the flow. This is the whole reason the load reads
 * `Ra·T − Rb·C` rather than `Ra·(T − B·C)` — a single ratio B = Rb/Ra is
 * undefined exactly where this case lives.
 */
describe("isothermal (Ra = 0) buoyancy", () => {
  const opts = { nr: 12, na: 24, gnr: 17, gna: 32, geom: ANNULUS, Ra: 0 } as const;

  it("Ra = 0 alone leaves ψ at zero — nothing left in the load to drive the flow", () => {
    const sim = new Simulation(opts);
    for (const row of sim.psi.c) for (const v of row) expect(v).toBeCloseTo(0, 12);
  });

  it("Rb alone still drives the flow when Ra = 0", () => {
    const sim = new Simulation({
      ...opts, particles: { seed: 9, layerDepth: 0.3 }, Rb: 1,
    });
    let max = 0;
    for (const row of sim.psi.c) for (const v of row) max = Math.max(max, Math.abs(v));
    expect(max).toBeGreaterThan(0);
  });

  // Tier 1 (isoviscous) solves a linear system for a fixed operator, so
  // doubling the load must double ψ exactly — the sharpest available check
  // that the compositional term is scaled by Rb alone, with nothing left of
  // the old B = Rb/Ra ratio quietly capping it.
  it("is linear in Rb, exactly, in the isoviscous tier", () => {
    const make = (Rb: number) =>
      new Simulation({ ...opts, particles: { seed: 9, layerDepth: 0.3 }, Rb });
    const one = make(1), two = make(2);
    for (let i = 0; i < one.psi.c.length; i++)
      for (let j = 0; j < one.psi.c[i].length; j++)
        expect(two.psi.c[i][j]).toBeCloseTo(2 * one.psi.c[i][j], 9);
  });
});

/**
 * van Keken et al. (1997)'s own Rayleigh–Taylor interface: a lighter fluid
 * (φ = 0) underlying a heavier one (φ = 1) across a boundary perturbed by
 * one cosine half-wavelength across the domain's width.
 */
describe("van Keken interface composition", () => {
  const g = box(2);
  const cond = SPECIES_CONDITIONS["van Keken interface"];
  const layerDepth = 0.2;

  it("switches from 0 to 1 exactly at the perturbed interface, at several φ", () => {
    for (const phi of [0, g.width / 4, g.width / 2, (3 * g.width) / 4]) {
      const boundary = layerDepth + 0.02 * Math.cos((Math.PI * phi) / g.width);
      expect(cond.composition(g, boundary - 1e-6, phi, layerDepth)).toBe(0);
      expect(cond.composition(g, boundary + 1e-6, phi, layerDepth)).toBe(1);
    }
  });

  // The same depth reads on opposite sides of the interface depending on φ
  // alone — the sharpest available check that the cosine term is actually
  // wired in, rather than a perturbation that silently evaluates to zero
  // everywhere and leaves a flat interface no φ argument could distinguish.
  it("the interface height genuinely depends on φ", () => {
    const atZero = layerDepth + 0.02 * Math.cos(0);              // φ = 0: +amplitude
    const atHalf = layerDepth + 0.02 * Math.cos(Math.PI / 2);    // φ = width/2: cos = 0
    expect(atZero).not.toBeCloseTo(atHalf, 6);
    const r = (atZero + atHalf) / 2;
    expect(cond.composition(g, r, 0, layerDepth)).toBe(0);              // below the φ = 0 interface
    expect(cond.composition(g, r, g.width / 2, layerDepth)).toBe(1);    // above the φ = width/2 interface
  });
});
