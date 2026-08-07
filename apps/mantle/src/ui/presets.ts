/**
 * Control state and the tables behind the list controls.
 *
 * Kept separate from `controls.ts` so the invariants can be regression-tested
 * without a DOM: an `N_φ` that is not a power of two, or a dt that does not fall
 * with the grid, are both mistakes that would otherwise surface only when a user
 * happened to pick that entry.
 */

import type { ColormapName } from "../colormaps";
import { annulus, box, type Geometry, type Walls } from "../geometry";
import {
  DEFAULT_LAYER_DEPTH, PARTICLE_TINT, type SpeciesConditionName, type TintMode,
} from "../particles";

/**
 * Resolution ladder. `N_φ` (both `na` and `gna`) must be a power of two — the
 * azimuthal FFT is radix-2 — and 1024 is the ceiling, where the
 * transform's shared-memory footprint reaches the 16 KB workgroup limit.
 *
 * `dtMax` falls with the radial grid for the reason a flat `dt` used to: a
 * coarse step spends a fine grid's headroom on a worse Courant number. It is
 * now a *ceiling* rather than the step itself — `main.ts` sizes the actual
 * step from the advective CFL every poll (see `adaptiveDt`) and never exceeds
 * this. Measured init and per-step costs on a discrete GPU are noted — the
 * growth is almost entirely in the f64 factorisation at start-up, not in the
 * frame loop.
 */
export const PRESETS = {
  "coarse · ψ 48×128": { nr: 48, na: 128, gnr: 97, gna: 128, dtMax: 2e-4 },
  "standard · ψ 96×256": { nr: 96, na: 256, gnr: 193, gna: 256, dtMax: 1e-4 },
  "fine · ψ 128×256": { nr: 128, na: 256, gnr: 257, gna: 512, dtMax: 7e-5 },
  "finer · ψ 160×512": { nr: 160, na: 512, gnr: 321, gna: 512, dtMax: 6e-5 },
  "finest · ψ 192×512": { nr: 192, na: 512, gnr: 385, gna: 512, dtMax: 5e-5 },
} as const;

export type PresetName = keyof typeof PRESETS;

/**
 * Steps advanced per animation frame. Fractional rates are the point: a coarse
 * mesh solves a step in well under a millisecond, so at one step per frame the
 * dynamics run far too fast to watch. Slowing down *must not* be done by
 * shrinking dt — that changes the trajectory being solved, not the rate it is
 * shown at — so the frame loop carries an accumulator instead and simply steps
 * less often.
 */
export const SPEEDS = {
  "1 step / 16 frames": 1 / 16,
  "1 step / 8 frames": 1 / 8,
  "1 step / 4 frames": 1 / 4,
  "1 step / 2 frames": 1 / 2,
  "1 step / frame": 1,
  "2 steps / frame": 2,
  "4 steps / frame": 4,
  "8 steps / frame": 8,
  "16 steps / frame": 16,
} as const;

/**
 * How much of the run the two corner plots — Nusselt number and RMS velocity —
 * show, as a span of solver **steps**. One list rather than two: both are
 * polls of the same frame loop at the same cadence, so a second "how much of
 * the run" control would offer nothing the first does not already set.
 *
 * Steps rather than samples, even though samples are what the trace stores: the
 * diagnostic poll rate is set by the frame loop, so a window of "the last 400
 * polls" would mean a different stretch of the simulation at every playback
 * speed, and would silently change under the user when they moved that list. A
 * step count is the simulation's own clock and the readout already displays it.
 *
 * Steps rather than nondimensional time for a related reason: `dt` moves with the
 * resolution preset, so a time window would show a different number of points on
 * a coarse mesh than a fine one. Either choice trades one invariance for another;
 * this one keeps the *sample density* of the plot fixed, which is what governs
 * whether the curve is a curve or a smear.
 *
 * `Infinity` shows everything the trace still holds — bounded by `NU_CAPACITY`,
 * not by this list, so a long enough run eventually rolls its own beginning off
 * the left edge whichever entry is selected.
 */
export const NU_WINDOWS = {
  "last 500 steps": 500,
  "last 2 000 steps": 2_000,
  "last 10 000 steps": 10_000,
  "last 50 000 steps": 50_000,
  "all": Infinity,
} as const;

/**
 * Viscosity law. `variable` picks the tier — whether the FFT
 * radial solve is the entire solve or the inner preconditioner — and
 * `strainRate` picks whether the power-law index `n` is the user's or pinned
 * to 1.
 *
 * **Only `variable` is a rebuild.** The Krylov tier allocates the
 * quadrature-point stress and viscosity buffers — tens of megabytes at the top
 * of the resolution ladder — so a constant-μ run must not carry them, and that
 * is a decision about allocation, which cannot be made per frame. The two
 * variable laws share every buffer and every pipeline: `n = 1` collapses the
 * power law to the identity exactly (see `rheology.ts`), so switching between
 * them is one uniform write.
 *
 * **Depth is not a fourth entry**, and the asymmetry with ε̇ is the point. The
 * strain rate is what makes the law *nonlinear* — it needs the strain kernels,
 * the Picard lag and a sweep count, so the two sides of it are genuinely
 * different solves and the list has to name which one is running. Depth is one
 * more factor in a coefficient that is still linear in ψ, costs one buffer read,
 * and is carried by the preconditioner exactly. So it is a slider inside both
 * variable laws — set it to zero contrast and each is the μ(T) law it was —
 * and the names below say so rather than the list doubling in length.
 *
 * **Tackley is a fifth exception to "only `variable` is a rebuild."** It sits
 * in the Krylov tier like the other two variable laws, but its pointwise law
 * is a different formula on a different set of parameters (σ_Y, σ_b, η*, no
 * γ/c/n), evaluated by a different GPU kernel — so `tackley` also gates a
 * rebuild, not a uniform write.
 *
 * **Tosi is a sixth exception, of the same shape as Tackley's.** Its
 * pointwise law is also a different GPU kernel, so `tosi` gates a rebuild
 * too — but unlike Tackley it *does* read γ and c (the paper's own γ_T, γ_z),
 * so it keeps the contrast/depth sliders rather than hiding them, and adds
 * σ_Y, σ_b, η* (shared with Tackley — the two laws state the identical
 * yielding branch).
 *
 * **Blankenbach is a seventh exception, of the same shape as Tosi's minus the
 * yielding branch.** It reads γ and c too — the paper's own b, c, referenced
 * at the cold surface rather than this app's T = ½, d = ½ (see
 * `rheology.ts`) — so it keeps the contrast/depth sliders, but its pointwise
 * law has no ε̇ dependence at all, so it hides n, the Picard sweeps and σ_Y,
 * σ_b, η* the same way `μ(T, d)` already does. The point of having it at all
 * rather than routing case 2a/2b through `μ(T, d)` at a rescaled Ra is that a
 * reader can then enter the paper's own Ra, b and c directly and get the
 * paper's own problem, with nothing to work out first.
 */
export const VISCOSITY = {
  "constant": {
    variable: false, strainRate: false, tackley: false, tosi: false,
    blankenbach: false, vanKeken: false,
  },
  "μ(T, d)": {
    variable: true, strainRate: false, tackley: false, tosi: false,
    blankenbach: false, vanKeken: false,
  },
  "μ(T, d, ε̇)": {
    variable: true, strainRate: true, tackley: false, tosi: false,
    blankenbach: false, vanKeken: false,
  },
  // Tackley (2000): Arrhenius creep and Bingham yielding in parallel — a
  // structurally different law, not a further tier of the power law above, so
  // it carries its own parameters (σ_Y, σ_b, η*) rather than γ, c, n.
  "Tackley": {
    variable: true, strainRate: true, tackley: true, tosi: false,
    blankenbach: false, vanKeken: false,
  },
  // Tosi et al. (2015): the community benchmark's harmonic combination of the
  // μ(T, d) exponential (as the paper's own γ_T, γ_z) and the same Bingham
  // yielding branch Tackley states. n means nothing to it (no power law), but
  // it is still ε̇-dependent and nonlinear, so it keeps the Picard sweeps.
  "Tosi": {
    variable: true, strainRate: true, tackley: false, tosi: true,
    blankenbach: false, vanKeken: false,
  },
  // Blankenbach et al. (1989): μ(T, d) = exp(−bT + cd), the paper's own
  // uncentred reference — see the note above and `rheology.ts`.
  "Blankenbach": {
    variable: true, strainRate: false, tackley: false, tosi: false,
    blankenbach: true, vanKeken: false,
  },
  // van Keken et al. (1997): μ(φ) = η_light + φ(η_dense − η_light) — the one
  // law here with no T dependence at all, so no ε̇ and no contrast sliders
  // either; it reads a tracer cloud's own composition instead (see
  // `rheology.ts`), which is what makes it the only law that stays correct
  // under a nonzero `Rb` (`Rb`'s own header in `gpu/wgsl.ts` explains why
  // every other law here would silently corrupt under that combination).
  "van Keken": {
    variable: true, strainRate: false, tackley: false, tosi: false,
    blankenbach: false, vanKeken: true,
  },
} as const;

export type ViscosityName = keyof typeof VISCOSITY;

/**
 * Three of the seven laws above, under names that owe nothing to the
 * literature — what the pane's "how the rock behaves" control offers by
 * default. The other four (Tackley, Tosi, Blankenbach, van Keken) are each
 * a specific paper's own formula and stay under "law" in the advanced
 * viscosity folder, named the way that paper names them; nothing here
 * removes them; it only decides what a first-time reader sees before they
 * ask for the rest.
 */
export const SIMPLE_VISCOSITY = {
  "same everywhere": "constant",
  "stiffer when cold": "μ(T, d)",
  "stiffer when cold, weaker where it's moving": "μ(T, d, ε̇)",
} as const satisfies Record<string, ViscosityName>;

export type SimpleViscosityName = keyof typeof SIMPLE_VISCOSITY;

/**
 * Mesh overlay, and the mode the render shader reads from the uniform.
 *
 * There are *two* meshes and they are not the same size, so one toggle would
 * have had to pick one and be wrong about the other half the time: ψ lives in a
 * spline space of `nr − 3` × `na` elements, T on a grid of `gnr − 1` × `gna`
 * cells — the readout names both. Neither is a rebuild; the overlay is drawn
 * from the element counts in the uniform, so this is one write.
 */
export const MESH = {
  "off": 0,
  "ψ elements": 1,
  "T grid": 2,
} as const;

export type MeshName = keyof typeof MESH;

/**
 * Domain. The annulus is the geometry the write-up derives and the one the
 * dimensional clock is scaled against; the box is the cell the mantle
 * convection literature states its benchmarks in — depth 0 → 1, hot below, and
 * an adjustable length.
 *
 * **A geometry is a rebuild**, and so is the length. The metric is emitted into
 * the WGSL rather than branched on a uniform (a box would otherwise evaluate
 * `1/r` at z = 0), and the length reaches the azimuthal knot vector — so the
 * quadrature tables, the discrete symbols and the per-mode inverses are all
 * built against it. There is no version of this that is a slider tracking the
 * pointer, which is why `boxLength` fires on release and says what it is doing.
 */
export const GEOMETRY = {
  "spherical annulus": "annulus",
  "Cartesian box": "box",
} as const;

export type GeometryName = keyof typeof GEOMETRY;

/**
 * Inner radius of the annulus, in units of the mantle's own thickness: Earth's
 * core–mantle boundary, 3486 km, divided by the 6371 − 3486 = 2885 km between
 * it and the surface. The outer radius is one unit further out by
 * construction, so the shell has unit thickness the same way the box has unit
 * depth. `ui/dimensional.ts` reads the same choice from the other side when it
 * puts years on the clock.
 */
export const RADIUS_INNER = 1.208318891;

/**
 * van Keken et al. (1997)'s own domain width for their Rayleigh–Taylor
 * benchmark suite, in units of the (unit) depth — chosen in the original
 * paper so that the fastest-growing linear mode of the instability is a
 * single cell across the box. It is the floor `BOX_LENGTH` opens at below,
 * not a round number: this is the one case in the benchmark table that
 * needs a domain narrower than it is deep.
 */
export const VAN_KEKEN_WIDTH = 0.9142;

/**
 * Bounds on the box length, in units of its depth.
 *
 * The floor used to be a flat 1 — a domain narrower than it is deep, where a
 * single cell cannot fit — until the van Keken Rayleigh–Taylor case needed
 * exactly that. `VAN_KEKEN_WIDTH` is now the floor itself rather than a
 * value the old floor excluded, so the single-cell argument still holds
 * everywhere above it. The ceiling is set by the *azimuthal* resolution,
 * which the preset ladder fixes — at L = 8 and `na = 256` a spline element
 * is 0.031 across against 0.008 in depth, four to one, and past that the
 * transverse direction is the accuracy bottleneck rather than the radial one
 * the ladder is sized on. 4 is the default because it is the aspect ratio
 * the box benchmarks are usually run at, and because it is close to the
 * annulus it sits beside in the list: at the mid radius that domain is
 * 2π·0.775 around by 0.45 deep.
 *
 * **`step` is far finer than the pane ever displays, on purpose** — the same
 * reason `CONTRAST`'s is: `pane.refresh()` snaps whatever `state` already
 * holds to the nearest step multiple the instant a benchmark is selected,
 * before the reader has touched anything, and a step of 0.5 would round
 * `VAN_KEKEN_WIDTH` itself to 1 on the spot. `format` is what actually keeps
 * the slider's own display at two decimals.
 */
export const BOX_LENGTH = {
  min: VAN_KEKEN_WIDTH, max: 8, step: 0.0001, default: 4,
  format: (v: number) => v.toFixed(2),
} as const;

/**
 * What closes the box left and right.
 *
 * *periodic* is what the transform gives natively — x wraps, and the domain has
 * no ends. *free-slip walls* is impermeable, stress-free and insulating at both,
 * reached by mirroring: the solve runs on twice the width and holds the state
 * even about x = 0, which `geometry.ts` explains is the walled problem itself
 * rather than a stand-in for it.
 *
 * The trade is **resolution**: a walled box of a given width is solved on twice
 * the period, so at a fixed `na` it resolves half as finely as a periodic one.
 * That is the whole price, and it is why periodic remains the default.
 *
 * The list is offered on the box alone. The annulus has no ends to close.
 */
export const WALLS = {
  "periodic": "periodic",
  "free-slip walls": "free-slip",
} as const satisfies Record<string, Walls>;

export type WallsName = keyof typeof WALLS;

/** The `Geometry` a `State` selects. */
export const geometryFor = (s: {
  geometry: GeometryName; boxLength: number; walls: WallsName;
}): Geometry =>
  GEOMETRY[s.geometry] === "annulus"
    ? annulus(RADIUS_INNER, RADIUS_INNER + 1)
    : box(s.boxLength, WALLS[s.walls]);

/**
 * Bounds on the two contrast sliders — log₁₀ of the viscosity ratio across
 * temperature, and across depth. Named constants, not inline numbers in
 * `controls.ts`, for the same reason `BOX_LENGTH` is: `BENCHMARKS` entries
 * are checked against them, so a slider range that moved without the check
 * following would fail at selection time rather than silently clipping a
 * benchmark's own contrast.
 *
 * **`step` is far finer than the two decimals `format` displays, on purpose.**
 * Tweakpane's own step constraint does not just gate the drag — `refresh()`
 * re-applies it to whatever is already in `state` and writes the rounded
 * result back (`ComplexValue.setRawValue`, via `InputBindingValue.fetch` →
 * `push`), which a benchmark preset hits immediately, before the reader ever
 * touches the slider. At the old 0.25 step, that rounded Blankenbach 2b's own
 * b = ln 16384 (log₁₀ 4.214…) to 4.25 on `pane.refresh()` alone — a contrast
 * 8.5% off the paper's own. `step: 0.0001` shrinks that silent rounding to a
 * few thousandths of a percent, well past where it could matter to the
 * physics; `format` is what actually keeps the display at two decimals, the
 * same split the Courant slider already draws between the two.
 */
export const CONTRAST = {
  min: 0, max: 5, step: 0.0001, format: (v: number) => v.toFixed(2),
} as const;
export const DEPTH_CONTRAST = {
  min: 0, max: 3, step: 0.0001, format: (v: number) => v.toFixed(2),
} as const;

/**
 * Labels of the rheology sliders, named once because two places must agree
 * on them: the pane, and the legend under the equation that tells the reader
 * which slider sets γ, which sets c and which sets n. Renaming a slider without
 * the legend following would leave it pointing at a control that is not there.
 *
 * The thermal one keeps the bare name it has always had. Adding "thermal" to it
 * would be more symmetric and less useful: it is the contrast the app opens with
 * and the one the write-up means, so it is the depth slider that has to say what
 * it is a contrast *across*.
 */
export const LABELS = {
  contrast: "log₁₀ contrast",
  depth: "log₁₀ depth contrast",
  n: "power-law n",
  sigmaY: "yield stress σ_Y",
  sigmaB: "yield-stress gradient σ_b",
  etaStar: "min. plastic viscosity η*",
  etaLight: "η light",
  etaDense: "η dense",
} as const;

/**
 * The tracer overlay. `off` never touches the GPU; `visual` and `chemical`
 * both keep a live cloud of tracers pushed through the flow, differing only
 * in whether the composition they carry is allowed to push back on the
 * buoyancy driving that same flow.
 *
 * The three-way list is a convenience over two independent facts, kept
 * distinct because they cost different things (see `gpu/particles.ts` and
 * `main.ts`'s `attachParticles`): `attached` decides whether a `GpuParticles`
 * object exists at all — its own buffers and pipelines, tens of milliseconds
 * to build or tear down — while `coupled` only decides whether the buoyancy
 * load's `T − B·C` term reads a live composition or is held at `B = 0`, a
 * single uniform write regardless of `attached`. So `visual ↔ chemical` is
 * free, and `off ↔` either of the other two is the one transition that pays
 * the construction cost.
 */
export const PARTICLES = {
  "off": { attached: false, coupled: false },
  "visual": { attached: true, coupled: false },
  "chemical": { attached: true, coupled: true },
} as const;

export type ParticlesName = keyof typeof PARTICLES;

/**
 * Tracer-count ladder. The ratio method (`particles.ts`, `gpu/particles.ts`)
 * wants on the order of 16 tracers per composition-grid cell for a tolerable
 * noise floor on the projected field; at the standard preset's composition
 * grid (97×128) that is ≈200 000, which is why it sits in the middle of the
 * ladder rather than at either end — the entries either side are for a
 * coarser or finer run than the default, not a knowingly under- or
 * over-sampled one.
 */
export const PARTICLE_COUNTS = {
  "50 000": 50_000,
  "100 000": 100_000,
  "200 000": 200_000,
  "400 000": 400_000,
  "800 000": 800_000,
} as const;

/**
 * Dot radius, in screen pixels. Large enough to read as a mark against the
 * field it is drawn over, small enough that a few hundred thousand of them
 * do not paint the canvas solid.
 */
export const PARTICLE_SIZE = { min: 0.5, max: 6, step: 0.1, default: 2.5 } as const;

/**
 * Dot opacity. Draw order is buffer order, which is arbitrary (see
 * `gpu/particles.ts`) — at a 2–3 px dot that never matters for one tracer
 * against another, but it is what makes a single dot's own alpha, not the
 * order it happens to be drawn in, the only handle on how a *dense* cloud
 * reads once many dots overlap.
 */
export const PARTICLE_OPACITY = { min: 0.05, max: 1, step: 0.05, default: 0.85 } as const;

/**
 * log₁₀ Rb — the compositional Rayleigh number's own slider coordinate,
 * spanning decades for the same reason `logRa`'s does. Rb weighs the
 * compositional term against the thermal one in the buoyancy load, which
 * reads `Ra·T − Rb·C` rather than `Ra·(T − B·C)` (see `particles.ts` and
 * `tqSource`/`bcSource` in `gpu/wgsl.ts`) — deliberately not a ratio of the
 * two, since a ratio cannot survive `Ra = 0`, the purely compositional
 * (isothermal) limit `State.isothermal` reaches. Rb = 1 is van Keken et
 * al.'s own isothermal Rayleigh–Taylor scale; Rb ~ 10⁵ is a dense layer's
 * compositional Rayleigh number entrained by an ordinary mantle-like Ra —
 * both physically reasonable settings, and the range covers both.
 */
export const LOG_RB = { min: -2, max: 6, step: 0.05, default: 0 } as const;

/**
 * Bounds on η_light and η_dense, the two endpoints of van Keken's own
 * composition-linear viscosity law (`rheology.ts`). Both default to 1 — the
 * paper's own isoviscous case 1a — and giving them a contrast reaches cases
 * 1b/1c without changing anything else about the law's shape.
 */
export const ETA_VAN_KEKEN = { min: 0.01, max: 100, step: 0.01, default: 1 } as const;

/**
 * Thickness of the dense basal layer, as a fraction of the mantle's depth —
 * the slider's own default mirrors `DEFAULT_LAYER_DEPTH` in `particles.ts`.
 */
export const LAYER_DEPTH = {
  min: 0.02, max: 0.5, step: 0.01, default: DEFAULT_LAYER_DEPTH,
} as const;

export interface State {
  geometry: GeometryName;
  /** Width of the Cartesian box, in units of its depth. Ignored by the annulus. */
  boxLength: number;
  /** What closes the box left and right. Ignored by the annulus. */
  walls: WallsName;
  /** log₁₀ Ra — the slider's coordinate, and the one the physics is smooth in. */
  logRa: number;
  /**
   * Force `Ra = 0` regardless of `logRa` — the purely compositional
   * ("isothermal", in the sense of no thermal expansion at all) limit van
   * Keken et al. (1997)'s own Rayleigh–Taylor case states. A separate
   * checkbox rather than widening `logRa` down towards it, for the same
   * reason `PARTICLES`'s `coupled` flag is separate from `Rb`'s own value:
   * `logRa` stays where it is genuinely useful (spanning the decades where
   * convective onset and plume count actually change) and this is a clean
   * override on top, not a slider stretched to cover an edge case it would
   * otherwise spend most of its travel on.
   */
  isothermal: boolean;
  /** Ceiling on the adaptive step — see `adaptiveDt` and `PRESETS`. */
  dtMax: number;
  /**
   * Target Courant number the adaptive step is sized to hold: `dt ≈ courant /
   * maxSpeed`. Matches `SimOptions.cfl`'s default in `solver/step.ts`.
   */
  courant: number;
  /** Steps per frame; may be < 1. See SPEEDS. */
  speed: number;
  paused: boolean;
  contours: number;
  lineWidth: number;
  mesh: MeshName;
  /** Temperature colour map — a render-pipeline rebuild, not a uniform. */
  colormap: ColormapName;
  /**
   * Span of the two corner plots (Nusselt number, RMS velocity), in solver
   * steps; `Infinity` for the whole trace.
   */
  nuWindow: number;
  wavenumber: number;
  resolution: PresetName;
  viscosity: ViscosityName;
  /** log₁₀ of the viscosity contrast μ(T=0)/μ(T=1); γ = ln of it. */
  logContrast: number;
  /**
   * log₁₀ of the viscosity contrast across the layer's depth, μ at the hot
   * boundary over μ at the cold one *at fixed T*; c = ln of it. Positive
   * stiffens the deep interior, which is the sign the mantle has.
   */
  logDepthContrast: number;
  /** Krylov iterations per solve. Fixed budget ⇒ predictable frame time. */
  iters: number;
  /** Power-law index. 1 is Newtonian; ≈3 is dislocation creep. */
  n: number;
  /** Rheology updates per solve. Each one costs a full Krylov budget. */
  picard: number;
  /** Constant ductile yield stress, Tackley and Tosi laws. */
  sigmaY: number;
  /** Gradient of brittle yield stress with depth, Tackley and Tosi laws. */
  sigmaB: number;
  /** Minimum plastic viscosity, Tackley and Tosi laws. */
  etaStar: number;
  /** Viscosity of the light material, van Keken's own composition-linear law. */
  etaLight: number;
  /** Viscosity of the dense material, van Keken's own composition-linear law. */
  etaDense: number;
  /** Show the text readout (domain, Ra, law, resolution, Nu, …) over the canvas. */
  debug: boolean;
  /** Tracer overlay mode — see `PARTICLES`. */
  particles: ParticlesName;
  /** Tracer count; see `PARTICLE_COUNTS`. Only ever allocated while `particles !== "off"`. */
  particleCount: number;
  /** How a tracer is coloured; see `PARTICLE_TINT` in `particles.ts`. */
  particleTint: TintMode;
  /** Initial composition profile; see `SPECIES_CONDITIONS` in `particles.ts`. Only read while `particles` is `"chemical"`. */
  particleSpecies: SpeciesConditionName;
  /**
   * The colour map `particleTint`'s row compiles in — tracked here for the
   * pane's legend, not independently chosen: the tint registry fixes one
   * map per mode, the same way `PARTICLE_TINT` documents.
   */
  particleColormap: ColormapName;
  /** Dot radius, screen pixels. */
  particleSize: number;
  /** Dot opacity. */
  particleOpacity: number;
  /** log₁₀ Rb. Only reaches the buoyancy load while `particles` is `"chemical"` — see `PARTICLES`, `LOG_RB`. */
  logRb: number;
  /** Thickness (or, for "van Keken interface", height) of the seeded composition boundary — see `SPECIES_CONDITIONS`. */
  layerDepth: number;
}

export const DEFAULT_PRESET: PresetName = "standard · ψ 96×256";

/**
 * Krylov budget, fixed so that frame time does not depend on the flow.
 *
 * Measured in f64 with the μ̄(r) preconditioner at a 10³ contrast, a *cold*
 * solve reduces the residual by 7e-2 in 8 iterations and 2e-3 in 16. The running
 * loop is never cold — it warm-starts from the previous frame — and there four
 * iterations already reach the floor f32 storage of ψ imposes. Twelve is
 * therefore margin: for the one cold solve at start-up or reseed, and for
 * raising the contrast a decade or two without revisiting this.
 */
export const DEFAULT_ITERS = 12;

/**
 * Illustrative starting points for a reader who does not yet know what a
 * Rayleigh number is — three named pictures rather than three numbers, so
 * "try an example" has something to offer before the literature benchmarks
 * below it do. Chosen for what they look like on screen, not for
 * reproducing a published figure — that difference is what keeps this table
 * separate from `BENCHMARKS`, whose every field is instead traceable to a
 * paper.
 *
 * All three leave `geometry` alone (the annulus the app opens on), so
 * picking one from the default state or from another quick start lands on
 * the same domain either way.
 */
export const QUICK_STARTS = {
  // Just past onset: one or two lazy cells, most of the layer doing very
  // little — the picture that says "convection" needs no more heat than
  // this to happen at all.
  "Sluggish mantle": {
    viscosity: "constant",
    logRa: Math.log10(2e3),
    wavenumber: 2,
  },
  // Three decades above onset: many small, fast plumes, the isoviscous
  // regime pushed as far as it goes before the two other quick starts add
  // temperature-dependent stiffening back in.
  "Vigorous convection": {
    viscosity: "constant",
    logRa: Math.log10(1e6),
    wavenumber: 8,
  },
  // A 10³ contrast between hot and cold rock is what turns isoviscous cells
  // into a few broad, persistent upwellings — the picture "plume" usually
  // means, and the reason `μ(T, d)` exists at all.
  "Whole-mantle plumes": {
    viscosity: "μ(T, d)",
    logRa: Math.log10(1e5),
    logContrast: Math.log10(1e3),
    logDepthContrast: 0,
    wavenumber: 4,
  },
} as const satisfies Record<string, Partial<State>>;

export type QuickStartName = keyof typeof QUICK_STARTS;

/**
 * Named parameter sets reproducing benchmarks from the literature. Each entry
 * is a partial `State` — only the fields the benchmark's definition actually
 * constrains — applied over whatever the pane currently holds.
 *
 * Deliberately not itself a `State` field: there is nothing here for the
 * solver to track once loaded (see the pane's own "snap back to custom" note
 * in `controls.ts`), so this table is consulted once, on selection, rather
 * than carried in the object `build()` reads every rebuild.
 *
 * Resolution and dt are never part of an entry: the ladder in `PRESETS`
 * governs accuracy independently of which physics problem is loaded, exactly
 * as `onResolution` and `onGeometry` already act as independent knobs — a
 * benchmark should not reach past the user's current choice there.
 */
export const BENCHMARKS = {
  // Blankenbach et al. (1989), cases 1a–1c: the same unit square, free-slip
  // on all four sides, isoviscous problem at three Rayleigh numbers a decade
  // apart — 10⁴, 10⁵, 10⁶ — so only `logRa` differs between them. All three
  // share the paper's own single-cell initial condition: at this aspect ratio
  // a higher seed mode settles onto a multi-cell steady state instead, which
  // is a different solution than the one each case reports.
  "Blankenbach 1a": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 4,
    viscosity: "constant",
    wavenumber: 1,
  },
  "Blankenbach 1b": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 5,
    viscosity: "constant",
    wavenumber: 1,
  },
  // The paper's steadiest-reported case; some later studies find it weakly
  // time-dependent, which is not a discrepancy this table can resolve — it
  // states the same problem 1a and 1b do, not a claim about its long-run
  // behaviour.
  "Blankenbach 1c": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 6,
    viscosity: "constant",
    wavenumber: 1,
  },
  // Blankenbach et al. (1989), case 2a: the same unit square and free-slip
  // walls, now temperature-dependent viscosity μ = exp(−bT) at the paper's
  // own Ra = 10⁴, b = ln 1000 — no depth term, that is case 2b's addition.
  //
  // Entered directly, with no rescaling: the `Blankenbach` law states the
  // paper's own uncentred formula (see `rheology.ts`), unlike this app's own
  // `μ(T, d)`, which centres its thermal exponent on T = ½ and so would need
  // Ra rescaled to reproduce this case. `logContrast: 3` is b = ln(10³) =
  // ln 1000, and `logDepthContrast: 0` is c = 0 — no depth term.
  "Blankenbach 2a": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 4,
    viscosity: "Blankenbach",
    logContrast: 3,
    logDepthContrast: 0,
    wavenumber: 1,
  },
  // Blankenbach et al. (1989), case 2b: 2a's temperature dependence plus a
  // depth term, at the paper's own Ra = 10⁴, b = ln 16384, c = ln 64 — *not*
  // 2a's b = ln 1000 with a depth term added on top, which is the guess the
  // shared name invites and would be wrong: 2b states its own, larger,
  // thermal contrast, independent of 2a's.
  //
  // The box is 2.5 wide, not the unit square 1a/2a use — the paper's own
  // cell geometry is wider for this case, presumably because the 64-fold
  // stiffening of the deep interior shifts the preferred cell width.
  // `logContrast`/`logDepthContrast` are computed rather than rounded
  // decimals so the value entered is exactly the paper's, not a
  // transcription of finitely many digits of it.
  "Blankenbach 2b": {
    geometry: "Cartesian box",
    boxLength: 2.5,
    walls: "free-slip walls",
    logRa: 4,
    viscosity: "Blankenbach",
    logContrast: Math.log10(16384),
    logDepthContrast: Math.log10(64),
    wavenumber: 1,
  },
  // van Keken et al. (1997), case 1: the isoviscous isothermal
  // Rayleigh–Taylor instability at the head of that paper's thermochemical
  // benchmark suite — a lighter fluid underlying a heavier one across a
  // perturbed interface, the unstable arrangement, with *no* thermal
  // buoyancy at all ("isothermal" in the sense that the thermal expansion
  // coefficient is zero, not that T is undefined — T is still advected and
  // diffused, it just never reaches the momentum balance). `isothermal:
  // true` forces Ra = 0 regardless of `logRa`, and the momentum source
  // reduces to the paper's own f = φ(0, −1)ᵀ, a unit compositional Rayleigh
  // number acting alone — see `isothermal`'s own header on why a checkbox
  // rather than a widened `logRa` floor.
  //
  // Domain width (0.9142), interface height (0.2) and the perturbation
  // itself (amplitude 1/50, one cosine half-wavelength across the width) are
  // the paper's own — see `VAN_KEKEN_WIDTH` and the "van Keken interface"
  // row of `SPECIES_CONDITIONS`. η_light = η_dense = 1 is case 1a, the
  // isoviscous one; giving the two a contrast reaches cases 1b/1c without
  // anything else here changing. `logRa`/`wavenumber` are entered anyway,
  // inert while `isothermal` is checked, so unchecking it lands on an
  // ordinary thermal run rather than Ra = 1 (log₁₀ 0) by accident.
  //
  // `dtMax` is the one field every other benchmark leaves at the resolution
  // preset's own value and this one cannot: that value is an accuracy
  // ceiling tuned to the O(10¹–10²) speeds a Ra = 10⁴ thermal run reaches
  // (`PRESETS`'s own header), and Rb = 1 acting alone tops out three to four
  // orders of magnitude slower — a GPU probe of this exact case measures
  // max|u| ≈ 0.07–0.12 through the growth phase. Left at the resolution's
  // own ceiling (2×10⁻⁴ at the default resolution), `adaptiveDt` never finds
  // a reason to raise it — `courant/maxSpeed` is ~10, nowhere near binding —
  // so the run advances at a fixed, tiny dt and looks stopped: reaching the
  // displacement the probe saw by t = 20 would take on the order of 10⁵
  // steps, minutes of wall clock even at the pane's top speed. 10⁻² is
  // still far below where `courant/maxSpeed` would start binding at these
  // speeds, so `adaptiveDt` stays free to shrink it again if the instability
  // goes on to accelerate past this phase.
  "van Keken 1997": {
    geometry: "Cartesian box",
    boxLength: VAN_KEKEN_WIDTH,
    walls: "free-slip walls",
    isothermal: true,
    logRa: Math.log10(1e4),
    wavenumber: 4,
    viscosity: "van Keken",
    etaLight: 1,
    etaDense: 1,
    particles: "chemical",
    particleSpecies: "van Keken interface",
    layerDepth: 0.2,
    logRb: 0,   // Rb = 1
    dtMax: 1e-2,
  },
  // Tosi et al. (2015), case 1: unit square, free-slip on all four sides, at
  // the paper's own Ra = 100, with a purely temperature-dependent linear
  // viscosity η_lin(T) = exp(−γ_T T), γ_T = ln 10⁵ — no depth term, and no
  // plastic branch at all in this case.
  //
  // η_lin is `Blankenbach`'s own exp(−bT + cd) with c = 0 (see
  // `rheology.ts`), so that law reproduces it exactly at the paper's own Ra,
  // with nothing to rescale first. `Tosi`'s harmonic combination is the
  // wrong tool for a case with no yielding: it does not recover η_lin as
  // σ_Y, η* → ∞, only 2η_lin — a harmonic *mean* of a finite viscosity
  // against an infinite one is twice the finite one, not the finite one.
  "Tosi 1": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 2,
    viscosity: "Blankenbach",
    logContrast: 5,
    logDepthContrast: 0,
    wavenumber: 1,
  },
  // Tosi et al. (2015), case 2: case 1's linear viscosity plus the paper's
  // own depth term, γ_z = ln 10 — still no yielding, so `Blankenbach` again,
  // now with its own c = γ_z.
  "Tosi 2": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 2,
    viscosity: "Tosi",
    logContrast: 5,
    logDepthContrast: 0,
    wavenumber: 1,
    sigmaY: 1,
    sigmaB: 0,
    etaStar: 1e-3,
  },
  // Tosi et al. (2015), case 3: case 1's temperature-dependent linear branch
  // plus the paper's own Bingham yielding branch — constant yield stress
  // σ_Y = 1, plastic-branch floor η* = 10⁻³ — combined the way the paper's
  // Eq. 7 states, which is what the `Tosi` law's harmonic combination is
  // for. `sigmaB: 0` states the paper's own *constant* yield stress
  // explicitly, rather than leaving it to whatever the pane currently holds:
  // σ_b is this app's own extension, and the paper's Cases 1–4 fix it to
  // zero (see `rheology.ts`).
  "Tosi 3": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 2,
    viscosity: "Blankenbach",
    logContrast: 5,
    logDepthContrast: 1,
    wavenumber: 1,
  },
  // Tosi et al. (2015), case 4: case 2's depth-dependent linear branch plus
  // case 3's yielding — the paper's most demanding case, reported bistable
  // between a mobile-lid and a stagnant-lid regime at these same parameters
  // depending on what the run is seeded with. This entry fixes the
  // parameters the paper states; which basin the wavenumber-1 initial
  // condition settles into is a property of the run, not of this table.
  "Tosi 4": {
    geometry: "Cartesian box",
    boxLength: 1,
    walls: "free-slip walls",
    logRa: 2,
    viscosity: "Tosi",
    logContrast: 5,
    logDepthContrast: 1,
    sigmaY: 1,
    sigmaB: 0,
    etaStar: 1e-3,
    wavenumber: 1,
  },
} as const satisfies Record<string, Partial<State>>;

export type BenchmarkName = keyof typeof BENCHMARKS;

export const defaultState = (): State => ({
  // n = 3 is dislocation creep, the mantle's dominant deformation mechanism.
  // picard = 1 is pure time-lagging: a second sweep
  // doubles the solve, and Stokes being quasi-static, the previous frame's
  // strain rate is already an O(dt) guess at this one's.
  //
  // The annulus opens the app because it is the geometry the write-up derives,
  // and the one whose free-slip condition is the point being made.
  geometry: "spherical annulus",
  boxLength: BOX_LENGTH.default,
  // Periodic by default: it is what the transform gives natively, and it spends
  // the azimuthal resolution on the domain rather than on its mirror image.
  walls: "periodic",
  logRa: Math.log10(1e4),
  isothermal: false,
  dtMax: PRESETS[DEFAULT_PRESET].dtMax,
  courant: 1.0,
  speed: 2,
  paused: false,
  // Both overlays start off: the temperature field is the subject, and the first
  // thing on screen should be it rather than a lattice drawn over it. Each is a
  // uniform write away.
  contours: 0,
  lineWidth: 1.1,
  mesh: "off",
  colormap: "inferno",
  // The whole trace by default. The plot's first job is the initial transient
  // settling, which is the one thing a fixed window would cut off — narrowing it
  // is what the reader does once they want to see the settled state resolved.
  nuWindow: NU_WINDOWS.all,
  wavenumber: 4,
  resolution: DEFAULT_PRESET,
  viscosity: "constant",
  logContrast: 3,
  // No depth dependence to start: it is the term the reader has to *ask* for,
  // and at zero each variable law is the μ(T) law the write-up derives.
  logDepthContrast: 0,
  iters: DEFAULT_ITERS,
  n: 3,
  picard: 1,
  sigmaY: 1,
  sigmaB: 1,
  etaStar: 1e-3,
  // Isoviscous by default — van Keken et al.'s own case 1a.
  etaLight: ETA_VAN_KEKEN.default,
  etaDense: ETA_VAN_KEKEN.default,
  // Off by default: the readout is a wall of numbers for anyone not
  // debugging, and sits over the top-left corner of the one thing the app is
  // actually showing.
  debug: false,
  // Off by default, for the same reason the two field overlays are: the
  // temperature field is the subject and the first thing on screen should be
  // it, not a cloud of dots drawn over it. `initial depth` is the tint worth
  // switching to once a reader turns the overlay on — colouring a parcel by
  // where it started, rather than by anything it currently is, is what makes
  // stirring visible (see `particles.ts`).
  particles: "off",
  particleCount: PARTICLE_COUNTS["200 000"],
  particleTint: "initial depth",
  particleColormap: PARTICLE_TINT["initial depth"].colormap,
  particleSpecies: "dense basal layer",
  particleSize: PARTICLE_SIZE.default,
  particleOpacity: PARTICLE_OPACITY.default,
  logRb: LOG_RB.default,
  layerDepth: LAYER_DEPTH.default,
});
