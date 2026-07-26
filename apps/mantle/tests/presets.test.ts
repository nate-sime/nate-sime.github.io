/**
 * The control tables. No GPU and no DOM: these are the invariants the
 * solver assumes and the pane silently depends on, and every one of them fails
 * at *selection* time rather than at start-up, so a bad entry could ship
 * unnoticed behind a list item nobody clicked.
 */

import { describe, it, expect } from "vitest";
import {
  MESH, PRESETS, SPEEDS, VISCOSITY, DEFAULT_PRESET, defaultState,
} from "../src/ui/presets";
import { P, clampedAxis, periodicAxis } from "../src/spline";

const pow2 = (n: number) => Number.isInteger(Math.log2(n));
const entries = Object.entries(PRESETS);

describe("resolution presets", () => {
  it.each(entries)("%s is solvable by the radix-2 azimuthal FFT", (_name, p) => {
    // Uniform periodic knots make the azimuthal operator circulant,
    // and the transform that diagonalises it is radix-2.
    expect(pow2(p.na)).toBe(true);
    expect(pow2(p.gna)).toBe(true);
    // 4·N floats of workgroup memory per transform against the 16 KB limit.
    expect(Math.max(p.na, p.gna)).toBeLessThanOrEqual(1024);
  });

  it.each(entries)("%s leaves interior degrees of freedom to solve for", (_name, p) => {
    expect(p.nr).toBeGreaterThan(P + 1);   // ψ = const drops the two boundary DOFs
    expect(p.gnr).toBeGreaterThan(5);      // the 5-point one-sided Nusselt stencil
  });

  // dt is an accuracy parameter, so it has to fall with the grid — otherwise a
  // refinement silently buys resolution and spends it on a worse Courant number.
  it("reduces dt monotonically along the ladder", () => {
    const dr = entries.map(([, p]) => 1 / (p.gnr - 1));
    const dt = entries.map(([, p]) => p.dt);
    for (let i = 1; i < entries.length; i++) {
      expect(dr[i]).toBeLessThan(dr[i - 1]);
      expect(dt[i]).toBeLessThan(dt[i - 1]);
      // and roughly in proportion, so the Courant number stays comparable
      expect(dt[i] / dt[i - 1] / (dr[i] / dr[i - 1])).toBeCloseTo(1, 0);
    }
  });

  it("starts from a preset that exists", () => {
    expect(PRESETS[DEFAULT_PRESET]).toBeDefined();
    expect(defaultState().dt).toBe(PRESETS[DEFAULT_PRESET].dt);
    expect(PRESETS[defaultState().resolution]).toBeDefined();
  });
});

/**
 * The viscosity table decides which changes cost a rebuild. `main.ts`
 * rebuilds on a change of `variable` and writes a uniform otherwise, so an entry
 * claiming strain-rate dependence without the Krylov tier would send `n` to a
 * solver with no μ buffer to put it in.
 */
describe("viscosity table", () => {
  it("only offers strain-rate dependence inside the Krylov tier", () => {
    const laws = Object.values(VISCOSITY);
    for (const { variable, strainRate } of laws)
      expect(strainRate && !variable).toBe(false);
    // One law per tier-and-law combination that exists, and a constant one: a
    // duplicate would make the list's two rebuild classes ambiguous.
    expect(laws.filter((l) => !l.variable)).toHaveLength(1);
    expect(new Set(laws.map((l) => `${l.variable}${l.strainRate}`)).size)
      .toBe(laws.length);
  });

  it("starts from a law that exists, with a usable exponent", () => {
    const s = defaultState();
    expect(VISCOSITY[s.viscosity]).toBeDefined();
    // n = 1 is the Newtonian identity; below it the exponent (1−n)/n flips sign
    // and the law would shear-*thicken*, which no slider should reach silently.
    expect(s.n).toBeGreaterThanOrEqual(1);
    expect(s.picard).toBeGreaterThanOrEqual(1);
  });
});

/**
 * The mesh overlay. Its values are not labels — the render shader discriminates
 * on them numerically (`> 0` for on, `< 1.5` for the spline mesh), so a
 * renumbered table would silently draw the wrong discretisation.
 */
describe("mesh table", () => {
  it("numbers the modes the way the shader reads them", () => {
    expect(MESH.off).toBe(0);
    const on = Object.entries(MESH).filter(([, v]) => v > 0);
    expect(on).toHaveLength(2);
    expect(new Set(Object.values(MESH)).size).toBe(Object.keys(MESH).length);
    // The shader's split point: 1 is the spline mesh, everything above it the
    // grid. Both must land on the side their label claims.
    expect(MESH["ψ elements"]).toBeLessThan(1.5);
    expect(MESH["T grid"]).toBeGreaterThan(1.5);
  });

  // The overlay divides the annulus into these counts, which the shader derives
  // from nothing — they are uploaded. This is the relation `gpu/sim.ts` reads off
  // the axes, stated where a change to `spline.ts` would be noticed.
  it.each(Object.entries(PRESETS))("spans %s with whole elements", (_name, p) => {
    expect(clampedAxis(p.nr, 0.55, 1).elements()).toHaveLength(p.nr - P);
    expect(periodicAxis(p.na).elements()).toHaveLength(p.na);
    // The two meshes are different, which is the whole reason both are offered.
    expect(p.nr - P).not.toBe(p.gnr - 1);
  });

  it("starts with both overlays off", () => {
    const s = defaultState();
    expect(MESH[s.mesh]).toBe(0);
    expect(s.contours).toBe(0);
    // …but with a usable width behind them, since neither draws without one.
    expect(s.lineWidth).toBeGreaterThan(0);
  });
});

describe("speed table", () => {
  it("spans both directions and is bounded", () => {
    const v = Object.values(SPEEDS);
    expect(Math.min(...v)).toBeLessThan(1);   // can be slowed below one step/frame
    expect(Math.max(...v)).toBeGreaterThan(1);
    // The frame loop's accumulator relies on this: at most `max` steps fall due
    // in any one frame, so a slow frame cannot cascade into a longer one.
    expect(Math.max(...v)).toBeLessThanOrEqual(16);
    expect(v.every((x) => x > 0)).toBe(true);
  });

  it("labels each rate with what it actually does", () => {
    for (const [label, v] of Object.entries(SPEEDS)) {
      const n = v >= 1 ? v : 1 / v;
      expect(label).toContain(String(n));
    }
    expect(Object.values(SPEEDS)).toContain(defaultState().speed);
  });
});
