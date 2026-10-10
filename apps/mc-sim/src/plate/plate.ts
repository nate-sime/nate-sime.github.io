// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The Kirchhoff plate, nondimensional, in a tensor-product spline space.
 *
 *   ∇²(D ∇²w) − (1 − ν)◇⁴(D, w) + μ ẅ = q + Σ P_j δ(x − a_j)   on (0, a) × (0, 1),
 *
 * lengths in units of the side b along y (so the plate is a × 1, a the aspect
 * ratio), bending stiffness D = Et³/12(1 − ν²) and mass μ = ρt relative to their
 * reference values. A uniform plate has D = μ = 1, and its frequencies ω̂ carry
 * the familiar Leissa scaling ω b² √(ρt/D). The weak form is
 *
 *   a(w, v) = ∫ D [ w_xx v_xx + w_yy v_yy + ν (w_xx v_yy + w_yy v_xx) + 2(1 − ν) w_xy v_xy ],
 *   m(w, v) = ∫ μ w v,   ℓ(v) = ∫ q v + Σ P_j v(a_j),
 *
 * which, like the beam's, asks for w ∈ H² — and the plate is where that hurts:
 * a C¹ Lagrange-type element on triangles is a famously awkward object (Argyris,
 * 21 dofs), while a tensor product of C¹ splines is just two beams' bases
 * multiplied. Degree p ≥ 2 and continuity k ≥ 1 in each direction, as for the
 * beam.
 *
 * Edges. Each edge is clamped (C), simply supported (S) or free (F); a name
 * lists them as left (x = 0), bottom (y = 0), right (x = a), top (y = 1), the
 * order Leissa uses. The essential conditions are imposed strongly, exactly as
 * on the beam: on an open knot vector w on x = 0 is the coefficient row i = 0
 * alone, and ∂w/∂x there involves rows 0 and 1 only. So
 *
 *   C ⇒ rows 0, 1 removed,   S ⇒ row 0,   F ⇒ nothing,
 *
 * and the natural conditions — zero normal moment on S and F, the Kirchhoff
 * shear and the corner forces on F, all of them ν-dependent — come from the
 * weak form with nothing more to do. The free coefficients form a rectangle
 * [lox, hix) × [loy, hiy) of the index grid.
 *
 * Numbering and bandwidth. Free coefficient (i, j) is unknown
 * (i − lox) + nfx (j − loy), x running fastest. Two functions interact when
 * both indices are within p, so the matrix is banded with half-width
 * p·nfx + p: N ≈ (ne + p)² unknowns, a band of ≈ p(ne + p), and a banded
 * Cholesky of O(N·bw²) = O(p² ne⁴) — h⁻⁴ against the beam's h⁻¹. That is the
 * cost the multilevel method has to beat in two dimensions, and the reason the
 * levels here stop at a few thousand unknowns.
 *
 * Coefficients arrive sampled at the quadrature points, as for the beam: the
 * grid of 1D Gauss points in x times that in y, flat with x fastest,
 * [gy·NX + gx] — the layout a separable random field fills cheaply
 * (`random/field2d.ts`).
 */

import { SymBand, cholSolve, cholesky } from "../band";
import { jitter } from "../beam/beam";
import { lowestModes, type Modes } from "../eig";
import { gaussLegendre } from "../quad";
import {
  basisOnElement, firstDof, locate, splineSpace, tabulate, type SplineSpace, type Tabulation,
} from "../spline";

export type Edge = "C" | "S" | "F";

/** Edge conditions, in Leissa's order: x = 0, y = 0, x = a, y = 1. */
export const EDGES = {
  SSSS: "simply supported",
  CCCC: "clamped",
  CSCS: "clamped on x = 0, a; simply supported on y = 0, 1",
  CFFF: "cantilever: clamped on x = 0, free elsewhere",
  FFFF: "free (Chladni's plate)",
} as const;

export type EdgeName = keyof typeof EDGES;

const edgesOf = (name: string): [Edge, Edge, Edge, Edge] => {
  if (!/^[CSF]{4}$/.test(name)) throw new Error(`edges are four of C, S, F; got "${name}"`);
  return name.split("") as [Edge, Edge, Edge, Edge];
};

const REMOVED: Record<Edge, number> = { C: 2, S: 1, F: 0 };

export interface PointLoad { readonly x: number; readonly y: number; readonly P: number; }

/** A plate coefficient: a function of (x, y), or its values at the quadrature grid. */
export type Coefficient2 = ((x: number, y: number) => number) | Float64Array;

export interface PlateLoad {
  readonly q?: Coefficient2;
  readonly forces?: readonly PointLoad[];
}

export interface PlateSpec {
  readonly p: number;
  /** Elements along each side. */
  readonly ne: number;
  /** Continuity; defaults to the maximal C^{p−1}. */
  readonly k?: number;
  readonly edges: string;
  /** Side along x, in units of the side along y. */
  readonly aspect?: number;
  /** Poisson's ratio. */
  readonly nu?: number;
  readonly stiffness?: Coefficient2;
  readonly mass?: Coefficient2;
  /** As for the beam: every assembled entry times (1 + δ), δ uniform in ±size — the round-off probe. */
  readonly perturb?: { readonly size: number; readonly seed: number };
}

export const DEFAULT_NU = 0.3;

/**
 * The two 1D tabulations a plate is assembled from — x and y — to second
 * derivatives at p + 1 Gauss points. One per level of a Monte Carlo run.
 */
export interface PlateTables {
  readonly tx: Tabulation;
  readonly ty: Tabulation;
}

export function plateTables(p: number, ne: number, k = p - 1, aspect = 1): PlateTables {
  const rule = gaussLegendre(p + 1);
  return { tx: tabulate(splineSpace(p, ne, k, 0, aspect), rule, 2), ty: tabulate(splineSpace(p, ne, k, 0, 1), rule, 2) };
}

/** The quadrature points of a plate's tables, flat [gy·NX + gx], as two coordinate arrays. */
export function gridPoints(t: PlateTables): { x: Float64Array; y: Float64Array } {
  const NX = t.tx.x.length, NY = t.ty.x.length, x = new Float64Array(NX * NY), y = new Float64Array(NX * NY);
  for (let gy = 0; gy < NY; gy++)
    for (let gx = 0; gx < NX; gx++) { x[gy * NX + gx] = t.tx.x[gx]; y[gy * NX + gx] = t.ty.x[gy]; }
  return { x, y };
}

/** Why a spec cannot be solved statically, or null if it can. */
export function admissible(spec: PlateSpec, modal = false): string | null {
  const k = spec.k ?? spec.p - 1;
  if (spec.p < 2 || k < 1)
    return "The plate's weak form needs second derivatives to be square-integrable, so the space must be at least C¹ (p ≥ 2, k ≥ 1).";
  const e = edgesOf(spec.edges);
  // Rigid motions are w = c₀ + c₁x + c₂y. A clamped edge stops all three, and
  // so do any two supported edges (two lines on which a plane vanishes force
  // it to zero); one simply supported edge leaves the rotation about it.
  const held = e.filter((v) => v !== "F").length;
  if (!modal && (held === 0 || (held === 1 && !e.includes("C"))))
    return held === 0
      ? "A free plate is a mechanism: it moves rigidly under any load, so it has modes (three rigid ones at ω = 0 among them) but no static solution."
      : "One simply supported edge leaves the plate free to rotate about it: a mechanism.";
  if (freeDofs(spec) < 1) return "The mesh is too coarse to leave any free coefficient.";
  return null;
}

/** The free index range along each direction, and the counts. */
function layout(spec: PlateSpec): { n: number; lox: number; hix: number; loy: number; hiy: number; nfx: number; nfy: number } {
  const n = splineSpace(spec.p, spec.ne, spec.k ?? spec.p - 1).n;
  const [l, b, r, t] = edgesOf(spec.edges);
  const lox = REMOVED[l], hix = n - REMOVED[r], loy = REMOVED[b], hiy = n - REMOVED[t];
  return { n, lox, hix, loy, hiy, nfx: Math.max(0, hix - lox), nfy: Math.max(0, hiy - loy) };
}

/** Free coefficients, without building the plate: the size of every solve. */
export function freeDofs(spec: PlateSpec): number {
  const { nfx, nfy } = layout(spec);
  return nfx * nfy;
}

/** Half-bandwidth of the free system, p·nfx + p. */
export function bandwidth(spec: PlateSpec): number {
  return spec.p * layout(spec).nfx + spec.p;
}

/**
 * Work of one banded factorisation, N·bw² multiply–adds up to a constant: the
 * unit multilevel Monte Carlo counts a plate sample's cost in.
 */
export const solveWork = (spec: PlateSpec): number => freeDofs(spec) * bandwidth(spec) ** 2;

export class Plate {
  readonly sx: SplineSpace;
  readonly sy: SplineSpace;
  readonly tables: PlateTables;
  readonly aspect: number;
  readonly nu: number;
  /** Coefficients per direction; the full vector is n × n, x fastest. */
  readonly n: number;
  readonly lox: number;
  readonly hix: number;
  readonly loy: number;
  readonly hiy: number;
  readonly nfx: number;
  readonly nfy: number;
  /** Free blocks, numbered (i − lox) + nfx (j − loy). */
  readonly K: SymBand;
  readonly M: SymBand;
  private LK?: SymBand;
  private readonly jitter?: () => number;

  constructor(readonly spec: PlateSpec, tables?: PlateTables) {
    const why = admissible(spec, true);
    if (why) throw new Error(why);
    const k = spec.k ?? spec.p - 1;
    this.aspect = spec.aspect ?? 1;
    this.nu = spec.nu ?? DEFAULT_NU;
    if (tables) {
      const { tx, ty } = tables;
      if (tx.space.p !== spec.p || tx.space.k !== k || tx.space.ne !== spec.ne || tx.space.b !== this.aspect ||
        ty.space.ne !== spec.ne || tx.d < 2 || ty.d < 2)
        throw new Error("the tables given are not of this plate's space");
    }
    this.tables = tables ?? plateTables(spec.p, spec.ne, k, this.aspect);
    this.sx = this.tables.tx.space;
    this.sy = this.tables.ty.space;
    const lay = layout(spec);
    this.n = lay.n;
    this.lox = lay.lox; this.hix = lay.hix; this.loy = lay.loy; this.hiy = lay.hiy;
    this.nfx = lay.nfx; this.nfy = lay.nfy;
    const D = sampled(this.tables, spec.stiffness ?? (() => 1));
    const mu = sampled(this.tables, spec.mass ?? (() => 1));
    const bw = spec.p * this.nfx + spec.p, N = this.nfx * this.nfy;
    this.K = new SymBand(N, Math.min(bw, Math.max(0, N - 1)));
    this.M = new SymBand(N, this.K.bw);
    this.assemble(D, mu);
    if (spec.perturb) {
      this.jitter = jitter(spec.perturb.size, spec.perturb.seed);
      this.K.data.forEach((v, i, a) => (a[i] = v * this.jitter!()));
      this.M.data.forEach((v, i, a) => (a[i] = v * this.jitter!()));
    }
  }

  /** Number of free coefficients. */
  get dofs(): number { return this.nfx * this.nfy; }

  /** Free unknown of coefficient (i, j), or −1 if it is constrained. */
  dof(i: number, j: number): number {
    return i < this.lox || i >= this.hix || j < this.loy || j >= this.hiy ? -1 : i - this.lox + this.nfx * (j - this.loy);
  }

  /** The full n × n coefficient vector (x fastest), zeros on the constrained ones. */
  expand(free: ArrayLike<number>): Float64Array {
    const c = new Float64Array(this.n * this.n);
    for (let j = this.loy; j < this.hiy; j++)
      for (let i = this.lox; i < this.hix; i++) c[j * this.n + i] = free[this.dof(i, j)];
    return c;
  }

  /**
   * Both matrices, element by element. At each quadrature point the 2D second
   * derivatives are products of 1D tables — B_xx = N″(x) N(y), B_xy = N′(x) N′(y),
   * B_yy = N(x) N″(y) — and the bending energy density is
   * D [B_xx B_xx + B_yy B_yy + ν (B_xx B_yy + B_yy B_xx) + 2(1 − ν) B_xy B_xy].
   */
  private assemble(D: Float64Array, mu: Float64Array): void {
    const { tx, ty } = this.tables, p = this.spec.p, P1 = p + 1, L = P1 * P1, nu = this.nu;
    const nq = tx.nq, NX = tx.x.length, sx = (tx.d + 1) * P1, sy = (ty.d + 1) * P1;
    const bxx = new Float64Array(L), byy = new Float64Array(L), bxy = new Float64Array(L), b0 = new Float64Array(L);
    const ke = new Float64Array(L * L), me = new Float64Array(L * L), ids = new Int32Array(L);
    for (let ey = 0; ey < this.sy.ne; ey++)
      for (let ex = 0; ex < this.sx.ne; ex++) {
        ke.fill(0);
        me.fill(0);
        for (let qy = 0; qy < nq; qy++)
          for (let qx = 0; qx < nq; qx++) {
            const gx = ex * nq + qx, gy = ey * nq + qy, g = gy * NX + gx;
            const w = tx.w[gx] * ty.w[gy], wd = w * D[g], wm = w * mu[g];
            const X = gx * sx, Y = gy * sy;
            for (let b = 0; b < P1; b++)
              for (let a = 0; a < P1; a++) {
                const l = b * P1 + a;
                b0[l] = tx.B[X + a] * ty.B[Y + b];
                bxx[l] = tx.B[X + 2 * P1 + a] * ty.B[Y + b];
                byy[l] = tx.B[X + a] * ty.B[Y + 2 * P1 + b];
                bxy[l] = tx.B[X + P1 + a] * ty.B[Y + P1 + b];
              }
            for (let r = 0; r < L; r++) {
              const xr = bxx[r], yr = byy[r], cr = bxy[r], mr = wm * b0[r];
              for (let s = 0; s <= r; s++) {
                ke[r * L + s] += wd * (xr * bxx[s] + yr * byy[s] + nu * (xr * byy[s] + yr * bxx[s]) + 2 * (1 - nu) * cr * bxy[s]);
                me[r * L + s] += mr * b0[s];
              }
            }
          }
        const fx = firstDof(this.sx, ex), fy = firstDof(this.sy, ey);
        for (let b = 0; b < P1; b++) for (let a = 0; a < P1; a++) ids[b * P1 + a] = this.dof(fx + a, fy + b);
        for (let r = 0; r < L; r++) {
          const I = ids[r];
          if (I < 0) continue;
          for (let s = 0; s <= r; s++) {
            const J = ids[s];
            if (J < 0) continue;
            this.K.add(I, J, ke[r * L + s]);
            this.M.add(I, J, me[r * L + s]);
          }
        }
      }
  }

  /** Load vector over the free coefficients. */
  load(load: PlateLoad): Float64Array {
    const F = new Float64Array(this.dofs), { tx, ty } = this.tables, P1 = this.spec.p + 1;
    if (load.q) {
      const q = sampled(this.tables, load.q), nq = tx.nq, NX = tx.x.length, sx = (tx.d + 1) * P1, sy = (ty.d + 1) * P1;
      for (let ey = 0; ey < this.sy.ne; ey++)
        for (let ex = 0; ex < this.sx.ne; ex++) {
          const fx = firstDof(this.sx, ex), fy = firstDof(this.sy, ey);
          for (let qy = 0; qy < nq; qy++)
            for (let qx = 0; qx < nq; qx++) {
              const gx = ex * nq + qx, gy = ey * nq + qy, wq = tx.w[gx] * ty.w[gy] * q[gy * NX + gx];
              for (let b = 0; b < P1; b++)
                for (let a = 0; a < P1; a++) {
                  const I = this.dof(fx + a, fy + b);
                  if (I >= 0) F[I] += wq * tx.B[gx * sx + a] * ty.B[gy * sy + b];
                }
            }
        }
    }
    for (const { x, y, P } of load.forces ?? []) {
      const ex = locate(this.sx, x), ey = locate(this.sy, y);
      const Nx = basisOnElement(this.sx, ex, x, 0), Ny = basisOnElement(this.sy, ey, y, 0);
      const fx = firstDof(this.sx, ex), fy = firstDof(this.sy, ey);
      for (let b = 0; b < P1; b++)
        for (let a = 0; a < P1; a++) {
          const I = this.dof(fx + a, fy + b);
          if (I >= 0) F[I] += P * Nx[a] * Ny[b];
        }
    }
    if (this.jitter) F.forEach((v, i) => (F[i] = v * this.jitter!()));
    return F;
  }

  /** Static deflection: free coefficients in, full coefficient vector out. */
  solve(load: PlateLoad | Float64Array): Float64Array {
    const F = load instanceof Float64Array ? load : this.load(load);
    this.LK ??= cholesky(this.K);
    return this.expand(cholSolve(this.LK, F));
  }

  /**
   * The `count` lowest modes, vectors full and M-orthonormal. With no edge held
   * the plate has three rigid modes at λ = 0 and K is singular, so the iteration
   * runs on K + σM, which shifts every eigenvalue by σ and none of the vectors.
   */
  modes(count: number): Modes {
    const free = !edgesOf(this.spec.edges).some((e) => e !== "F");
    const sigma = free ? 1 : 0;
    let K = this.K;
    if (sigma) {
      K = new SymBand(this.K.n, this.K.bw, this.K.data.slice());
      K.data.forEach((v, i, a) => (a[i] = v + sigma * this.M.data[i]));
    }
    const r = lowestModes(K, this.M, count);
    return { values: r.values.map((v) => v - sigma), vectors: r.vectors.map((v) => this.expand(v)), iterations: r.iterations };
  }

  /** w at (x, y) from a full coefficient vector. */
  evaluate(c: ArrayLike<number>, x: number, y: number): number {
    const ex = locate(this.sx, x), ey = locate(this.sy, y), P1 = this.spec.p + 1;
    const Nx = basisOnElement(this.sx, ex, x, 0), Ny = basisOnElement(this.sy, ey, y, 0);
    const fx = firstDof(this.sx, ex), fy = firstDof(this.sy, ey);
    let w = 0;
    for (let b = 0; b < P1; b++)
      for (let a = 0; a < P1; a++) w += c[(fy + b) * this.n + fx + a] * Nx[a] * Ny[b];
    return w;
  }

  /** w on the grid xs × ys, flat [iy·xs.length + ix]: separable, so O(n p) per point, not per point and function. */
  grid(c: ArrayLike<number>, xs: ArrayLike<number>, ys: ArrayLike<number>): Float64Array {
    return gridValues(this.sx, this.sy, c, xs, ys);
  }
}

/** Σ c_ij B_i(x) B_j(y) on a grid, for any tensor space — the coarse level of a pair too. */
export function gridValues(
  sx: SplineSpace, sy: SplineSpace, c: ArrayLike<number>, xs: ArrayLike<number>, ys: ArrayLike<number>,
): Float64Array {
  const n = sx.n, P1 = sx.p + 1, nx = xs.length;
  // T[j·nx + ix] = Σ_i c_ij B_i(x_ix), then w = Σ_j B_j(y) T[j].
  const T = new Float64Array(sy.n * nx);
  for (let ix = 0; ix < nx; ix++) {
    const e = locate(sx, xs[ix]), N = basisOnElement(sx, e, xs[ix], 0), f0 = firstDof(sx, e);
    for (let j = 0; j < sy.n; j++) {
      let s = 0;
      for (let a = 0; a < P1; a++) s += c[j * n + f0 + a] * N[a];
      T[j * nx + ix] = s;
    }
  }
  const out = new Float64Array(nx * ys.length);
  for (let iy = 0; iy < ys.length; iy++) {
    const e = locate(sy, ys[iy]), N = basisOnElement(sy, e, ys[iy], 0), f0 = firstDof(sy, e);
    for (let b = 0; b < P1; b++) {
      const nb = N[b], o = (f0 + b) * nx;
      for (let ix = 0; ix < nx; ix++) out[iy * nx + ix] += nb * T[o + ix];
    }
  }
  return out;
}

function sampled(t: PlateTables, c: Coefficient2): Float64Array {
  if (typeof c !== "function") return c;
  const { x, y } = gridPoints(t);
  return Float64Array.from(x, (xi, i) => c(xi, y[i]));
}
