/**
 * Passive-marker particles, f64: the reference the GPU particle push and the
 * chemical-composition projection are checked against, exactly as
 * `Temperature` is the reference `advectSource` is checked against
 * (`wgsl.ts`'s own header). A tracer answers two different questions with the
 * one mechanism: run with `species` never referenced, it draws the *pathline*
 * a parcel of mantle actually follows — the thing ψ's instantaneous
 * streamlines cannot show in a time-dependent flow. Read its `c`, and the
 * same cloud is a marker-in-cell discretisation of a chemically distinct
 * layer, carried with no numerical diffusion at all, which is the point of
 * using markers rather than advecting a field.
 *
 * Particles live in (r, φ), the one representation decision the plan makes —
 * see `geometry.ts` and `particles.ts` for why.
 */

import type { Geometry } from "../geometry";
import { mat } from "../linalg";
import { cubic, type Velocity } from "./temperature";
import {
  seedParticles, DEFAULT_LAYER_DEPTH, TARGET_PARTICLES_PER_CELL,
  SPECIES_CONDITIONS, type SpeciesConditionName,
} from "../particles";

/**
 * φ folded into one period. Periodic domains (the annulus, and a box with
 * periodic side walls) wrap, the way the domain itself repeats. A free-slip
 * box does not: x = 0 and x = width are solid, impermeable walls, so nothing
 * in continuous flow ever reaches them with nonzero speed, and a particle
 * that numerically overshoots one by the tracer path's own truncation error
 * belongs reflected back into the domain — the same treatment the radial
 * clamp below gives an analogous overshoot at the mantle's own top and
 * bottom, not a wraparound into what would be its own mirror image.
 */
export function closePhi(geom: Geometry, phi: number): number {
  if (geom.walls !== "free-slip") {
    const span = geom.span;
    return ((phi % span) + span) % span;
  }
  const w = geom.width, period = 2 * w;
  const p = ((phi % period) + period) % period;
  return p > w ? period - p : p;
}

/**
 * Forward explicit-midpoint RK2 trace of one tracer through the flow —
 * forward in the literal sense: unlike `Temperature.departure`, which traces
 * *backward* from a grid point to the place its value came from, this traces
 * a parcel of mantle *forward* from where it is to where the flow carries it
 * next, which is what a tracer's path (a pathline) actually is.
 *
 * The one line this changes from `departure` run at a negated dt is the
 * final transverse step: it is taken using h at the RK2 *midpoint* radius,
 * not the starting radius. On the annulus, where h = r genuinely changes
 * over a step, using the wrong h leaves an O(dt²) error every step — first
 * order over a long-run trajectory. A single semi-Lagrangian hop, one cell
 * long, never accumulates enough of that to matter next to its own
 * interpolation error; a tracer pushed for thousands of steps has no
 * interpolation error to hide behind, so the same shortcut would show up as
 * particles slowly drifting off the streamlines they should be following. In
 * a box h is constant and the distinction is silently absent — which is why
 * `Temperature.departure` itself is left exactly as it is: fixing it there
 * would move every Nusselt number in the test suite for a term the
 * semi-Lagrangian scheme cannot see past anyway.
 */
export function advanceRK2(
  geom: Geometry, u: Velocity, r: number, phi: number, dt: number,
): [number, number] {
  const a = u(r, phi);
  const rm = Math.min(geom.hi, Math.max(geom.lo, r + 0.5 * dt * a.ur));
  const pm = phi + 0.5 * dt * (a.up / geom.h(r));
  const b = u(rm, pm);
  return [r + dt * b.ur, phi + dt * (b.up / geom.h(rm))];
}

export interface ParticleOptions {
  /** Tracer count. Defaults to `TARGET_PARTICLES_PER_CELL` per composition-grid cell. */
  count?: number;
  /** Names the draw — the same seed reproduces the same cloud. See `particles.ts`. */
  seed?: number;
  /** Initial composition profile; see `SPECIES_CONDITIONS`. */
  species?: SpeciesConditionName;
  /** Thickness of the dense layer, as a fraction of the mantle's depth. */
  layerDepth?: number;
  /** Composition-grid resolution. Default: half the temperature grid, each direction. */
  cnr?: number;
  cna?: number;
}

/**
 * A cloud of tracers plus the composition field they project onto.
 *
 * The composition grid is coarser than, and independent of, the temperature
 * grid: `C` enters the buoyancy load only through the load's own quadrature,
 * which projects onto the ψ spline space regardless, so a grid a quarter the
 * size resolves everything the flow can respond to while asking a quarter as
 * many tracers to reach the noise floor a marker-in-cell estimate needs. See
 * the plan's particle-budget argument for the arithmetic.
 */
export class Particles {
  r!: Float64Array;
  phi!: Float64Array;
  c!: Float64Array;
  /** Composition field, marker-in-cell projected. Unpopulated cells hold their last value — see `project`. */
  readonly C: Float64Array[];
  readonly cnr: number;
  readonly cna: number;
  readonly cdr: number;
  readonly cdphi: number;
  private readonly species: SpeciesConditionName;
  private readonly layerDepth: number;

  constructor(
    readonly geom: Geometry, gnr: number, gna: number, o: ParticleOptions = {},
  ) {
    this.cnr = o.cnr ?? Math.floor((gnr + 1) / 2);
    this.cna = o.cna ?? gna / 2;
    this.cdr = (geom.hi - geom.lo) / (this.cnr - 1);
    this.cdphi = geom.span / this.cna;
    this.C = mat(this.cnr, this.cna);
    this.species = o.species ?? "dense basal layer";
    this.layerDepth = o.layerDepth ?? DEFAULT_LAYER_DEPTH;
    const count = o.count ?? TARGET_PARTICLES_PER_CELL * this.cnr * this.cna;
    this.seed(count, o.seed ?? 1);
  }

  /**
   * Draw a fresh cloud, area-uniform, and paint the species profile onto it —
   * what pressing "reseed" means. Projects immediately, so a caller reading
   * `C` right after `seed` — including this class's own constructor, ahead of
   * the simulation's first Stokes solve — sees this cloud's composition
   * rather than a stale or zeroed field.
   */
  seed(count: number, seedValue: number): void {
    const drawn = seedParticles(this.geom, count, seedValue);
    // Widened, not reinterpreted: seedParticles draws in f32 so the GPU cloud
    // starts bit-identical, and every f32 value is exactly representable in
    // f64, so this loses nothing while giving the f64 push room to diverge
    // from the GPU's only at the rate its own arithmetic does.
    this.r = Float64Array.from(drawn.r);
    this.phi = Float64Array.from(drawn.phi);
    const cond = SPECIES_CONDITIONS[this.species];
    this.c = new Float64Array(count);
    for (let p = 0; p < count; p++)
      this.c[p] = cond.composition(this.geom, this.r[p], this.layerDepth);
    this.project();
  }

  /**
   * Advance every tracer one step along the flow, in place. A particle push
   * is a pointwise read-modify-write — unlike the temperature grid, a
   * tracer's new position depends on no neighbour's, so there is nothing to
   * double-buffer and nothing to race.
   */
  push(u: Velocity, dt: number): void {
    const { geom } = this;
    for (let p = 0; p < this.r.length; p++) {
      const [rn, pn] = advanceRK2(geom, u, this.r[p], this.phi[p], dt);
      this.r[p] = Math.min(geom.hi, Math.max(geom.lo, rn));
      this.phi[p] = closePhi(geom, pn);
    }
  }

  /**
   * Marker-in-cell projection of the tracer composition onto `C`: each
   * particle deposits its `c` onto the four composition-grid nodes bracketing
   * it, weighted by the usual bilinear (area) weights, and a node's value is
   * the weighted mean of everything that landed near it — self-normalising
   * against however many tracers happen to be nearby, so a cell that
   * momentarily holds twice the average count is not read as twice the
   * composition.
   *
   * A cell no tracer touched this step keeps whatever composition it last
   * held rather than being reset to zero. Zero is not "no data" here, it is
   * "no dense material" — snapping an unvisited cell to it and back as
   * tracers wander past would be a spike in the buoyancy load with no
   * physical cause, worse than the mild staleness of leaving it alone.
   *
   * The free-slip mirror is imposed exactly as `Temperature.diffuse` imposes
   * it for T: tracers only ever occupy `[0, width]`, so the upper half of the
   * period is written as the reflection of the lower rather than left to
   * whatever nothing scattered there.
   */
  project(): void {
    const { cnr, cna, cdr, cdphi, geom } = this;
    const num = mat(cnr, cna), den = mat(cnr, cna);
    for (let p = 0; p < this.r.length; p++) {
      const x = Math.min(cnr - 1, Math.max(0, (this.r[p] - geom.lo) / cdr));
      const i = Math.min(cnr - 2, Math.floor(x)), tr = x - i;
      let y = this.phi[p] / cdphi;
      y -= Math.floor(y / cna) * cna;
      const j = Math.floor(y), tp = y - j, jn = (j + 1) % cna;
      const c = this.c[p];
      const w00 = (1 - tr) * (1 - tp), w01 = (1 - tr) * tp;
      const w10 = tr * (1 - tp), w11 = tr * tp;
      num[i][j] += w00 * c; den[i][j] += w00;
      num[i][jn] += w01 * c; den[i][jn] += w01;
      num[i + 1][j] += w10 * c; den[i + 1][j] += w10;
      num[i + 1][jn] += w11 * c; den[i + 1][jn] += w11;
    }
    for (let i = 0; i < cnr; i++)
      for (let j = 0; j < cna; j++)
        if (den[i][j] > 0) this.C[i][j] = num[i][j] / den[i][j];

    if (geom.walls === "free-slip")
      for (let i = 0; i < cnr; i++)
        for (let j = cna / 2 + 1; j < cna; j++)
          this.C[i][j] = this.C[i][cna - j];
  }

  /**
   * Monotone bicubic read-back of the projected composition at an arbitrary
   * (r, φ) — `Temperature.sample`, generalised over the composition grid's
   * own (coarser) spacing rather than the temperature grid's. This is what
   * lets the buoyancy load (`buoyancyLoad` in `step.ts`) read C at exactly
   * the quadrature points it already samples T at, so `T − B·C` is one
   * effective field rather than two fields assembled on different rules.
   */
  sample(r: number, phi: number): number {
    const { cnr, cna, cdr, cdphi, geom, C } = this;
    const x = Math.min(cnr - 1, Math.max(0, (r - geom.lo) / cdr));
    const i = Math.min(cnr - 2, Math.floor(x)), tr = x - i;
    let y = phi / cdphi;
    y -= Math.floor(y / cna) * cna;
    const j = Math.floor(y), tp = y - j;

    const col = new Float64Array(4);
    for (let m = -1; m <= 2; m++) {
      const row = C[Math.min(cnr - 1, Math.max(0, i + m))];
      col[m + 1] = cubic(row[(j - 1 + cna) % cna], row[j % cna],
        row[(j + 1) % cna], row[(j + 2) % cna], tp);
    }
    return cubic(col[0], col[1], col[2], col[3], tr);
  }
}
