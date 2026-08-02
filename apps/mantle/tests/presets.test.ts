/**
 * The control tables. No GPU and no DOM: these are the invariants the
 * solver assumes and the pane silently depends on, and every one of them fails
 * at *selection* time rather than at start-up, so a bad entry could ship
 * unnoticed behind a list item nobody clicked.
 */

import { describe, it, expect } from "vitest";
import {
  BOX_LENGTH, GEOMETRY, MESH, NU_WINDOWS, PRESETS, RADIUS_INNER, SPEEDS,
  VISCOSITY, WALLS, DEFAULT_PRESET, defaultState, geometryFor,
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

  // dtMax is an accuracy ceiling, so it has to fall with the grid — otherwise a
  // refinement silently buys resolution and spends it on a worse Courant number.
  it("reduces dtMax monotonically along the ladder", () => {
    const dr = entries.map(([, p]) => 1 / (p.gnr - 1));
    const dtMax = entries.map(([, p]) => p.dtMax);
    for (let i = 1; i < entries.length; i++) {
      expect(dr[i]).toBeLessThan(dr[i - 1]);
      expect(dtMax[i]).toBeLessThan(dtMax[i - 1]);
      // and roughly in proportion, so the Courant number stays comparable
      expect(dtMax[i] / dtMax[i - 1] / (dr[i] / dr[i - 1])).toBeCloseTo(1, 0);
    }
  });

  it("starts from a preset that exists", () => {
    expect(PRESETS[DEFAULT_PRESET]).toBeDefined();
    expect(defaultState().dtMax).toBe(PRESETS[DEFAULT_PRESET].dtMax);
    expect(PRESETS[defaultState().resolution]).toBeDefined();
  });
});

/**
 * The geometry table. Every entry here is consumed by something that cannot
 * report its own mistake: the metric is *compiled into* the WGSL, the box length
 * goes into the azimuthal knot vector, and the depth limits go into the radial
 * one — so an entry that named a domain the solver could not build would fail at
 * pipeline creation with a shader error, several layers from its cause.
 */
describe("geometry table", () => {
  const names = Object.keys(GEOMETRY) as (keyof typeof GEOMETRY)[];
  const wallNames = Object.keys(WALLS) as (keyof typeof WALLS)[];
  // Every combination the pane can produce, including the ones the annulus
  // ignores — `geometryFor` must not care which is which.
  const combos = names.flatMap((geometry) =>
    wallNames.map((walls) => [geometry, walls] as const));

  it.each(combos)("builds a solvable %s with %s edges", (name, walls) => {
    const g = geometryFor({ geometry: name, boxLength: BOX_LENGTH.default, walls });
    expect(g.hi).toBeGreaterThan(g.lo);
    expect(g.span).toBeGreaterThan(0);
    expect(g.lo).toBeGreaterThanOrEqual(0);
    // The width is what is drawn and the span is what is solved; the first is
    // never the larger, and they part company only for a mirrored box.
    expect(g.width).toBeLessThanOrEqual(g.span);
  });

  it.each(combos)("puts %s's hot boundary at its low end (%s)", (name, walls) => {
    // `lo` is hot and `hi` is cold everywhere in the solver — the buoyancy load,
    // the Dirichlet rows and the colour map all assume it — so the conduction
    // profile must run 1 → 0 in that direction, whichever domain it is.
    const g = geometryFor({ geometry: name, boxLength: BOX_LENGTH.default, walls });
    expect(g.conduction(g.lo)).toBeCloseTo(1, 12);
    expect(g.conduction(g.hi)).toBeCloseTo(0, 12);
  });

  it("gives the box the depth the literature states its benchmarks in", () => {
    const g = geometryFor({ geometry: "Cartesian box", boxLength: 4, walls: "periodic" });
    expect(g.lo).toBe(0);
    expect(g.hi).toBe(1);
    expect(g.span).toBe(4);
    // h = 1 and h′ = 0 is the entire difference from the annulus, and it is what
    // lets the lower boundary sit at z = 0 without anything dividing by it.
    expect(g.h(0)).toBe(1);
    expect(g.dh).toBe(0);
  });

  it("keeps the annulus at unit thickness, the dimensional clock's assumption", () => {
    const g = geometryFor({ geometry: "spherical annulus", boxLength: 4, walls: "periodic" });
    expect(g.lo).toBe(RADIUS_INNER);
    expect(g.hi - g.lo).toBeCloseTo(1, 12);
    expect(g.span).toBeCloseTo(2 * Math.PI, 12);
  });

  /**
   * Both boundaries of both domains must read Nu = 1 for pure conduction — that
   * is what `nuScale` is *for*, and it is the one relation in `geometry.ts` that
   * a plausible-looking wrong constant would leave silently intact everywhere
   * else. Integrating the conductive flux by hand here rather than through
   * `Temperature` keeps this a statement about the geometry alone.
   */
  it.each(combos)(
    "normalises %s (%s) so that conduction is Nu = 1 at both ends", (name, walls) => {
      const g = geometryFor({ geometry: name, boxLength: 3, walls });
      const na = 32, dphi = g.span / na, eps = (g.hi - g.lo) * 1e-6;
      // −dT_c/dr at each boundary, by a central difference just inside it.
      const flux = (r: number) =>
        (g.conduction(r - eps) - g.conduction(r + eps)) / (2 * eps);
      for (const r of [g.lo, g.hi])
        expect(g.nuScale(r) * flux(r) * na * dphi).toBeCloseTo(1, 5);
    });

  it("offers box lengths that span an aspect ratio worth having", () => {
    expect(BOX_LENGTH.min).toBeGreaterThan(0);
    expect(BOX_LENGTH.max).toBeGreaterThan(BOX_LENGTH.min);
    // The default has to land on a step of the slider, or the pane opens on a
    // value it cannot return to.
    const steps = (BOX_LENGTH.default - BOX_LENGTH.min) / BOX_LENGTH.step;
    expect(steps).toBeCloseTo(Math.round(steps), 9);
    expect(BOX_LENGTH.default).toBeGreaterThanOrEqual(BOX_LENGTH.min);
    expect(BOX_LENGTH.default).toBeLessThanOrEqual(BOX_LENGTH.max);
  });

  /**
   * The wall setting reaches the solver as one number — the *period* the knot
   * vector and every transform are built on — and reaches the renderer as
   * another, the width it draws. Getting the factor of two backwards would put
   * a walled box at half the resolution it claims, or draw a periodic one twice.
   */
  it("solves a walled box on twice the width it draws", () => {
    const wide = geometryFor({
      geometry: "Cartesian box", boxLength: 3, walls: "free-slip walls",
    });
    expect(wide.walls).toBe("free-slip");
    expect(wide.width).toBe(3);
    expect(wide.span).toBe(6);

    // …and a periodic one on exactly the width it draws.
    const wrap = geometryFor({
      geometry: "Cartesian box", boxLength: 3, walls: "periodic",
    });
    expect(wrap.walls).toBe("periodic");
    expect(wrap.width).toBe(3);
    expect(wrap.span).toBe(3);
  });

  /**
   * Nu is the horizontal *mean* flux, and for a mirrored box that mean is taken
   * over the doubled period — which is why `nuScale = 1/span` needs no case of
   * its own. The two boxes below have the same physical width and must therefore
   * report the same Nu from the same per-point flux, at half the sample spacing.
   */
  it("means the Nusselt flux over the period, not the drawn width", () => {
    const wrap = geometryFor({
      geometry: "Cartesian box", boxLength: 3, walls: "periodic",
    });
    const wall = geometryFor({
      geometry: "Cartesian box", boxLength: 3, walls: "free-slip walls",
    });
    const q = 1.7;   // any uniform boundary flux
    const na = 64;
    expect(wall.nuScale(0) * q * na * (wall.span / na))
      .toBeCloseTo(wrap.nuScale(0) * q * na * (wrap.span / na), 12);
  });

  it("has no walls to close on the annulus", () => {
    // Selecting a wall setting must not reach a domain with no ends: the
    // annulus is periodic whatever the list says, or the FFT would be solving
    // something the pane never showed.
    for (const walls of wallNames)
      expect(geometryFor({
        geometry: "spherical annulus", boxLength: 4, walls,
      }).walls).toBe("periodic");
  });

  it("starts from a geometry the list offers", () => {
    const s = defaultState();
    expect(GEOMETRY[s.geometry]).toBeDefined();
    expect(s.boxLength).toBe(BOX_LENGTH.default);
    expect(WALLS[s.walls]).toBeDefined();
  });
});

/**
 * The viscosity table decides which changes cost a rebuild. `main.ts`
 * rebuilds on a change of `variable`, `tackley` or `tosi`, and writes a
 * uniform otherwise, so an entry claiming strain-rate dependence without the
 * Krylov tier would send `n` to a solver with no μ buffer to put it in.
 */
describe("viscosity table", () => {
  it("only offers strain-rate dependence inside the Krylov tier", () => {
    const laws = Object.values(VISCOSITY);
    for (const { variable, strainRate } of laws)
      expect(strainRate && !variable).toBe(false);
    // One law per tier-and-law combination that exists, and a constant one: a
    // duplicate would make the list's rebuild classes ambiguous. Tackley,
    // Tosi and μ(T, d, ε̇) all share (variable, strainRate) — all three are
    // nonlinear, Krylov-tier laws — so `tackley` and `tosi` are part of the
    // key too; they are what actually decide the rebuild between those
    // three, since each compiles a different kernel.
    expect(laws.filter((l) => !l.variable)).toHaveLength(1);
    expect(new Set(laws.map((l) => `${l.variable}${l.strainRate}${l.tackley}${l.tosi}`)).size)
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

/**
 * The Nusselt plot's display window. Its values are step counts the plot compares
 * against the trace's own step axis, so a label and a value that disagree would
 * show a window other than the one named — and `Infinity` is load-bearing, not a
 * placeholder: it is what `NuTrace.first` reads as "no cutoff".
 */
describe("Nu window table", () => {
  it("offers an unbounded entry, and finite ones that are whole steps", () => {
    const v = Object.values(NU_WINDOWS);
    expect(v.filter((x) => !Number.isFinite(x))).toEqual([Infinity]);
    for (const x of v.filter(Number.isFinite)) {
      expect(Number.isInteger(x)).toBe(true);
      expect(x).toBeGreaterThan(0);
    }
    expect(new Set(v).size).toBe(v.length);
  });

  it("labels each window with the step count it applies", () => {
    for (const [label, v] of Object.entries(NU_WINDOWS)) {
      if (!Number.isFinite(v)) continue;
      // Thin spaces or not, the digits in the label must be the value's.
      expect(label.replace(/\D/g, "")).toBe(String(v));
    }
  });

  it("widens monotonically, so the list reads as a range", () => {
    const v = Object.values(NU_WINDOWS);
    for (let i = 1; i < v.length; i++) expect(v[i]).toBeGreaterThan(v[i - 1]);
  });

  // The default has to be an entry the list contains, or the pane opens showing a
  // selection that is not in its own options.
  it("starts from a window the list offers", () => {
    expect(Object.values(NU_WINDOWS)).toContain(defaultState().nuWindow);
  });

  /**
   * The narrowest window must still hold enough samples to be a curve. Polls
   * arrive every 15 frames, so the sample count over a window is
   * `steps / (15 · speed)` — at the fastest playback and the smallest window,
   * that has to stay above a handful of points.
   */
  it("keeps the narrowest window plottable at the fastest playback", () => {
    const narrowest = Math.min(...Object.values(NU_WINDOWS).filter(Number.isFinite));
    const fastest = Math.max(...Object.values(SPEEDS));
    expect(narrowest / (15 * fastest)).toBeGreaterThan(1);
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
