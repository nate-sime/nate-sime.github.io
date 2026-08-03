/**
 * Shared ground for the particle feature: the tracer/species colour registry,
 * area-uniform seeding, and the dense-layer initial condition. Nothing here
 * touches the GPU or the f64 solver — it is read by both (`gpu/particles.ts`
 * and `solver/particles.ts`, once they exist) and by the pane, the way
 * `colormaps.ts` is read by the fragment shader and the colour-bar legend
 * alike. One table per decision, so a mode or a species profile is a row
 * added here rather than a code path added somewhere else.
 *
 * Particles are stored and pushed in (r, φ) — see `geometry.ts` for why that
 * is the one representation decision this feature makes — so every physical
 * quantity below is stated on that pair rather than on world (x, y).
 */

import type { ColormapName } from "./colormaps";
import type { Geometry } from "./geometry";

// ---- reproducible seeding -----------------------------------------------------

/** A stream of draws in [0, 1), advanced one call at a time. */
export type Rng = () => number;

/**
 * mulberry32: a small, fast generator, not a cryptographic one — the property
 * that matters here is that the same seed reproduces the same cloud of
 * tracers on every run, not that the draws are unpredictable. That
 * reproducibility is what lets "reseed" be a named, sharable initial state
 * rather than a fresh scatter every time, and it is what lets the CPU and GPU
 * particle sets start bit-identical for the parity tests in §8 of the plan:
 * both draw from the same stream, so `Math.random` — seedless and
 * unreproducible — is the one thing this cannot use.
 */
export const mulberry32 = (seed: number): Rng => {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) | 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
};

/**
 * Fraction of the mantle's thickness a radius has descended, 0 at the hot
 * boundary and 1 at the cold one — the same convention `Geometry.lo`/`hi`
 * fix, just linear in r rather than the conductive profile `conduction(r)` is.
 * This is the depth a "seeded near the core–mantle boundary" or "top quarter
 * of the layer" statement means, and it is what both the *initial depth*
 * colour mode and the dense-layer species profile below are stated in terms
 * of.
 */
export const depthOf = (g: Geometry, r: number): number => (r - g.lo) / (g.hi - g.lo);

/**
 * Draw a radius so that tracers are uniform in *area*, not in r. On the
 * annulus, area grows linearly with r, so a naive uniform-in-r draw crowds
 * 24% more tracers per unit area against the core–mantle boundary than
 * against the surface — the seeding would still average to the right
 * composition, but the noise on that average would be depth-dependent, and a
 * quiet reading near one boundary next to a noisy one near the other reads as
 * a physical signal it is not.
 *
 * `xi` is one uniform draw in [0, 1); the map below is the inverse of the
 * cumulative area fraction out to r, i.e. the one draw that lands a particle
 * at r with probability proportional to the area of that radial shell. In a
 * box that area is independent of r (`dh = 0`) and the map is linear; on the
 * annulus it grows with r itself (`dh = 1`) and the map is the square root
 * that undoes it. Both are the same statement — invert the area measure —
 * specialised the way `geometry.ts`'s own `h` specialises everything else.
 */
export const seedRadius = (g: Geometry, xi: number): number =>
  g.dh === 0
    ? g.lo + xi * (g.hi - g.lo)
    : Math.sqrt(g.lo * g.lo + xi * (g.hi * g.hi - g.lo * g.lo));

/** One seeded cloud: parallel arrays, ready to upload as the GPU particle buffer's r and φ columns. */
export interface ParticleSeed {
  r: Float32Array;
  phi: Float32Array;
}

/**
 * Seed `count` tracers, area-uniform in (r, φ), from a named `seed` so the
 * result is reproducible. The transverse draw is uniform over `g.width`
 * rather than `g.span`: on a free-slip box the solved period is folded about
 * x = 0, and the far half is a mirror image, not new ground — seeding it too
 * would place every particle twice. On the annulus and a periodic box `width
 * === span` and the distinction is silent.
 */
export const seedParticles = (g: Geometry, count: number, seed: number): ParticleSeed => {
  const rng = mulberry32(seed);
  const r = new Float32Array(count);
  const phi = new Float32Array(count);
  for (let i = 0; i < count; i++) {
    r[i] = seedRadius(g, rng());
    phi[i] = rng() * g.width;
  }
  return { r, phi };
};

// ---- the chemical species: initial conditions ---------------------------------

/** Default thickness of the dense basal layer, as a fraction of the mantle's depth. */
export const DEFAULT_LAYER_DEPTH = 0.2;

/** One way to paint the two-species field onto a freshly seeded cloud. */
export interface SpeciesCondition {
  label: string;
  composition(g: Geometry, r: number, layerDepth: number): number;
}

/**
 * Registry of composition initial conditions. One entry to start — a dense,
 * chemically distinct layer sitting on the hot (basal) boundary, the setup
 * the thermochemical-plume and LLSVP literature runs, whose fate under a
 * given buoyancy ratio (stays a stable blanket, domes up, or gets stirred
 * into the overturn and entrained) is the question §5 of the plan exists to
 * let a particle answer without numerical diffusion answering it instead. A
 * second initial condition is a second row, not a new code path.
 */
export const SPECIES_CONDITIONS = {
  "dense basal layer": {
    label: "dense basal layer",
    composition: (g, r, layerDepth) => (depthOf(g, r) < layerDepth ? 1 : 0),
  },
} as const satisfies Record<string, SpeciesCondition>;

export type SpeciesConditionName = keyof typeof SPECIES_CONDITIONS;

// ---- colour: one attribute, a registry of how to fill it ----------------------

/**
 * One entry in the tracer-colour registry. `wgsl` is `null` for a mode that
 * needs no push-time computation — either because the colour is fixed
 * (`uniform`) or because it was written into the particle's `a` once, at
 * seeding, on the host (`initial depth`, `initial φ`); otherwise it is the
 * WGSL expression for `a`, evaluated by the push kernel every step from the
 * particle's *current* (r, φ) — the mechanism that makes `temperature` and
 * `speed` track a tracer's present surroundings rather than the place it
 * started. `colormap` is which of `colormaps.ts`'s maps the render pipeline
 * compiles in for that mode; `discrete` marks a mode whose map should read as
 * distinct bands (the two species) rather than a continuous scale.
 */
export interface TintEntry {
  label: string;
  wgsl: string | null;
  colormap: ColormapName;
  discrete?: boolean;
}

/**
 * How a tracer is coloured, as one table read by the push kernel, the render
 * pipeline and the pane alike. The default worth reaching for is *initial
 * depth*: colouring a parcel by where it started, rather than by anything it
 * currently is, is what makes stirring visible — two tracers seeded side by
 * side that end up on opposite sides of the cell say something a field plot
 * on a fixed grid cannot.
 *
 * `species` is the coupled mode, and deliberately not a separate rendering
 * path: it is one more row whose attribute is the composition a particle
 * already carries (`p.c`) and whose map is a two-stop discrete one, which is
 * the sense in which "coloured by chemical species" and "coloured by depth"
 * are the same mechanism.
 */
export const PARTICLE_TINT = {
  uniform: { label: "uniform", wgsl: null, colormap: "viridis" },
  "initial depth": { label: "initial depth", wgsl: null, colormap: "viridis" },
  "initial φ": { label: "initial φ", wgsl: null, colormap: "turbo" },
  temperature: { label: "temperature", wgsl: "sample_T(r, phi)", colormap: "inferno" },
  speed: { label: "speed", wgsl: "length(velocity(r, phi)) / stat[3]", colormap: "plasma" },
  age: { label: "age", wgsl: "min(1.0, pt.age / pt.tau)", colormap: "magma" },
  species: { label: "species", wgsl: "p.c", colormap: "coolwarm", discrete: true },
} as const satisfies Record<string, TintEntry>;

export type TintMode = keyof typeof PARTICLE_TINT;

// ---- constants shared with the GPU scatter and the CPU twin -------------------

/**
 * Fixed-point scale for the composition projection's atomic accumulators
 * (plan §5): WGSL has integer atomics only, so the cloud-in-cell weights are
 * scaled by this and rounded before `atomicAdd`. Chosen so a cell would need
 * upward of 4096 tracers before the accumulator could overflow, against a
 * design target of 8–32 tracers per cell — headroom, not a limit that is ever
 * meant to be approached. The CPU twin uses ordinary floating-point sums and
 * has no need of this scale itself, but shares the constant so a convergence
 * test can hold the two projections to the same rounding budget.
 */
export const CIC_FIXED_POINT_SCALE = 1 << 20;

/**
 * Tracers per composition-grid cell the ratio method wants for a tolerable
 * noise floor on the projected composition — the number both the CPU twin's
 * default tracer count and (later) the pane's particle-count ladder are sized
 * from, so the two never quote different budgets for the same picture.
 */
export const TARGET_PARTICLES_PER_CELL = 16;
