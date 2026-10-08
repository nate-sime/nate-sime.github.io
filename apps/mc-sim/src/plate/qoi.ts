// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The plate's quantities of interest and its hierarchy — the beam's five,
 * carried to two dimensions (`hierarchy.ts` says what each one is):
 *
 *   ω₁          √λ₁ of the plate,
 *   deflection  w at the centre — at the middle of the free edge, x = a, for
 *               the cantilever plate CFFF — where the point load also acts,
 *   compliance  ℓ(w),
 *   field       ‖w‖ over the plate,
 *   response    |w| at the same point under the load applied harmonically,
 *               with Rayleigh damping fitted to the uniform plate's first two
 *               modes.
 *
 * Level ℓ is ne₀·2^ℓ elements along each side, so dofs grow as h⁻² and the
 * work of a banded solve as h⁻⁴: γ = 4, against the beam's 1.
 */

import { gaussLegendre } from "../quad";
import { complexSolve } from "../band";
import {
  DEFAULT_FORCING, climb, type Forcing, type Harmonic, type Hierarchy, type LoadCase, type QoI,
} from "../hierarchy";
import { exactPlateEigenvalues, navierResponse, navierStatic } from "./exact";
import { Plate, gridValues, solveWork, type EdgeName, type PlateLoad, type PlateSpec } from "./plate";

export interface PlateCase {
  readonly edges: EdgeName;
  readonly load: LoadCase;
  /** Side along x over side along y. */
  readonly aspect: number;
  readonly nu: number;
  readonly forcing?: Forcing;
}

/** Where deflection is read and the point load acts. */
export const plateQoiPoint = (c: { edges: string; aspect: number }): { x: number; y: number } =>
  c.edges === "CFFF" ? { x: c.aspect, y: 0.5 } : { x: c.aspect / 2, y: 0.5 };

/** The load of a case: uniform pressure q (1, or a sampled field), or a point force P at the QoI point. */
export function plateLoadOf(c: PlateCase, q: Float64Array | null = null, P = 1): PlateLoad {
  return c.load === "uniform" ? { q: q ?? (() => 1) } : { forces: [{ ...plateQoiPoint(c), P }] };
}

const references = new Map<string, Float64Array>();

/**
 * The uniform plate's two lowest eigenvalues: exact where the edges allow,
 * else from quartic splines on 16 × 16 elements, solved once per case — good
 * to 10⁻⁶ or better, and all a damping fit needs.
 */
export function plateEigenvalues2(c: { edges: string; aspect: number; nu: number }): [number, number] {
  const ex = exactPlateEigenvalues(c.edges, c.aspect, 2);
  if (ex) return [ex[0], ex[1]];
  const key = `${c.edges}/${c.aspect}/${c.nu}`;
  let v = references.get(key);
  if (!v) {
    const modes = new Plate({ p: 4, ne: 16, edges: c.edges, aspect: c.aspect, nu: c.nu }).modes(c.edges === "FFFF" ? 5 : 2).values;
    v = c.edges === "FFFF" ? modes.slice(3) : modes;
    references.set(key, v);
  }
  return [v[0], v[1]];
}

/** Ω and the Rayleigh coefficients, fitted as for the beam: ζ on the uniform plate's first two modes. */
export function plateHarmonicOf(c: PlateCase): Harmonic {
  const { ratio, zeta } = c.forcing ?? DEFAULT_FORCING;
  const [l1, l2] = plateEigenvalues2(c);
  const w1 = Math.sqrt(l1), w2 = Math.sqrt(l2);
  return { Omega: ratio * w1, a: (2 * zeta * w1 * w2) / (w1 + w2), b: (2 * zeta) / (w1 + w2) };
}

/** Q for one solved plate; `c` the full coefficients of the static solution (or of Re u), `ci` of Im u. */
export function evaluatePlateQoI(
  plate: Plate, qoi: QoI, load: PlateLoad, at: { x: number; y: number }, h?: Harmonic,
): { Q: number; c: Float64Array; ci?: Float64Array } {
  if (qoi === "omega1") return { Q: Math.sqrt(plate.modes(1).values[0]), c: new Float64Array(0) };
  const F = plate.load(load);
  if (qoi === "response") {
    if (!h) throw new Error("the forced response needs a forcing frequency and damping");
    const { Omega: W, a, b } = h;
    const u = complexSolve(plate.K, plate.M, [1, W * b], [-W * W, W * a], F);
    const c = plate.expand(u.re), ci = plate.expand(u.im);
    return { Q: Math.hypot(plate.evaluate(c, at.x, at.y), plate.evaluate(ci, at.x, at.y)), c, ci };
  }
  const c = plate.solve(F);
  if (qoi === "deflection") return { Q: plate.evaluate(c, at.x, at.y), c };
  if (qoi === "compliance") {
    // ℓ(w) = Fᵀc over the free coefficients.
    let C = 0;
    for (let j = plate.loy; j < plate.hiy; j++)
      for (let i = plate.lox; i < plate.hix; i++) C += F[plate.dof(i, j)] * c[j * plate.n + i];
    return { Q: C, c };
  }
  return { Q: plateFieldNorm(plate, c, null), c };
}

/**
 * ‖w‖ — or ‖w − w_c‖ with w_c on a coarser plate of the same family — by Gauss
 * with p + 3 points per direction on every element: the coarse mesh nests in
 * the fine, so the integrand is a polynomial there and the rule exact.
 */
export function plateFieldNorm(
  plate: Plate, c: Float64Array, coarse: { plate: Plate; c: Float64Array } | null, ref?: (xs: Float64Array, ys: Float64Array) => Float64Array,
): number {
  const { xs, ys, wx, wy } = fineRule(plate);
  const w = plate.grid(c, xs, ys);
  if (coarse) {
    const wc = gridValues(coarse.plate.sx, coarse.plate.sy, coarse.c, xs, ys);
    for (let i = 0; i < w.length; i++) w[i] -= wc[i];
  }
  if (ref) {
    const we = ref(xs, ys);
    for (let i = 0; i < w.length; i++) w[i] -= we[i];
  }
  let s = 0;
  for (let iy = 0; iy < ys.length; iy++)
    for (let ix = 0; ix < xs.length; ix++) s += wx[ix] * wy[iy] * w[iy * xs.length + ix] ** 2;
  return Math.sqrt(s);
}

function fineRule(plate: Plate): { xs: Float64Array; ys: Float64Array; wx: Float64Array; wy: Float64Array } {
  const g = gaussLegendre(plate.spec.p + 3);
  const along = (breaks: Float64Array) => {
    const ne = breaks.length - 1, x = new Float64Array(ne * g.x.length), w = new Float64Array(x.length);
    for (let e = 0; e < ne; e++) {
      const a = breaks[e], h = breaks[e + 1] - a;
      g.x.forEach((xi, q) => { x[e * g.x.length + q] = a + h * xi; w[e * g.x.length + q] = h * g.w[q]; });
    }
    return { x, w };
  };
  const X = along(plate.sx.breaks), Y = along(plate.sy.breaks);
  return { xs: X.x, ys: Y.x, wx: X.w, wy: Y.w };
}

export interface PlateHierarchySpec {
  readonly p: number;
  readonly k: number;
  readonly plate: PlateCase;
  readonly qoi: QoI;
  readonly ne0: number;
  readonly levels: number;
}

/** Q of the continuum plate, where a closed form gives it. */
export function plateExact(c: PlateCase, qoi: QoI): number | null {
  if (qoi === "omega1") {
    const ev = exactPlateEigenvalues(c.edges, c.aspect, 1);
    return ev ? Math.sqrt(ev[0]) : null;
  }
  // The field's error needs w pointwise, which the series gives cheaply on a grid for a uniform load only.
  if (c.edges !== "SSSS" || (qoi === "field" && c.load === "point")) return null;
  const at = plateQoiPoint(c);
  if (qoi === "response") return navierResponse(c.aspect, c.load, at, plateHarmonicOf(c));
  const ns = navierStatic(c.aspect, c.load === "uniform" ? "uniform" : at);
  if (qoi === "deflection") return ns.w(at.x, at.y);
  if (qoi === "compliance") return ns.compliance;
  return ns.rms;
}

/**
 * What theory predicts for α where the solution is smooth: the simply
 * supported plate, whose modes and uniform-load deflection are trigonometric.
 * Elsewhere a corner where a clamped or free edge meets another carries a
 * singular term r^s that caps the rate below these, and the hierarchy measures
 * what it is. A point load makes w ~ r² log r under it: no smooth-data rate.
 */
export function plateTheoryRate(qoi: QoI, p: number, c: PlateCase): number | null {
  if (c.edges !== "SSSS") return null;
  if (qoi === "omega1") return 2 * (p - 1);
  if (c.load === "point") return null;
  if (qoi === "compliance") return 2 * (p - 1);
  if (qoi === "field") return Math.min(p + 1, 2 * (p - 1));
  return null;
}

export function runPlateHierarchy(spec: PlateHierarchySpec): Hierarchy {
  const { p, k, plate: pc, qoi } = spec;
  const at = plateQoiPoint(pc), load = plateLoadOf(pc);
  const harmonic = qoi === "response" ? plateHarmonicOf(pc) : undefined;
  const exact = plateExact(pc, qoi);
  const ns = qoi === "field" && pc.edges === "SSSS" && pc.load === "uniform" ? navierStatic(pc.aspect, "uniform") : null;
  type S = { plate: Plate; c: Float64Array };
  return climb<S>({
    p, exact,
    solve: (ne, perturb) => {
      const ps: PlateSpec = { p, k, ne, edges: pc.edges, aspect: pc.aspect, nu: pc.nu, perturb };
      const plate = new Plate(ps);
      const { Q, c } = evaluatePlateQoI(plate, qoi, load, at, harmonic);
      return { Q, dofs: plate.dofs, work: solveWork(ps), state: { plate, c } };
    },
    diff: (f, c) => (qoi === "field" ? plateFieldNorm(f.state.plate, f.state.c, c.state) : Math.abs(f.Q - c.Q)),
    // The field's error is ‖w_h − w‖, with w the Navier series on the same Gauss grid.
    error: (r) => (ns ? plateFieldNorm(r.state.plate, r.state.c, null, ns.grid!) : Math.abs(r.Q - exact!)),
  }, spec.ne0, spec.levels);
}
