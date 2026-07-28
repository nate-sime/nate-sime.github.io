/**
 * Control state and the tables behind the list controls.
 *
 * Kept separate from `controls.ts` so the invariants can be regression-tested
 * without a DOM: an `N_φ` that is not a power of two, or a dt that does not fall
 * with the grid, are both mistakes that would otherwise surface only when a user
 * happened to pick that entry.
 */

import { annulus, box, type Geometry } from "../geometry";

/**
 * Resolution ladder. `N_φ` (both `na` and `gna`) must be a power of two — the
 * azimuthal FFT is radix-2 — and 1024 is the ceiling, where the
 * transform's shared-memory footprint reaches the 16 KB workgroup limit.
 *
 * `dt` falls with the radial grid so the advective Courant number stays roughly
 * constant across the ladder; it is an accuracy parameter, not a speed knob
 * (that is `speed`, below). Measured init and per-step costs on a discrete GPU
 * are noted — the growth is almost entirely in the f64 factorisation at start-up,
 * not in the frame loop.
 */
export const PRESETS = {
  "coarse · ψ 48×128": { nr: 48, na: 128, gnr: 97, gna: 128, dt: 2e-4 },
  "standard · ψ 96×256": { nr: 96, na: 256, gnr: 193, gna: 256, dt: 1e-4 },
  "fine · ψ 128×256": { nr: 128, na: 256, gnr: 257, gna: 512, dt: 7e-5 },
  "finer · ψ 160×512": { nr: 160, na: 512, gnr: 321, gna: 512, dt: 6e-5 },
  "finest · ψ 192×512": { nr: 192, na: 512, gnr: 385, gna: 512, dt: 5e-5 },
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
 * How much of the run the Nusselt plot shows, as a span of solver **steps**.
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
 */
export const VISCOSITY = {
  "constant": { variable: false, strainRate: false },
  "μ(T)": { variable: true, strainRate: false },
  "μ(T, ε̇)": { variable: true, strainRate: true },
} as const;

export type ViscosityName = keyof typeof VISCOSITY;

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
 * Radius ratio of the annulus: Earth's core–mantle boundary against the
 * surface, 3480/6371 = 0.546. `ui/dimensional.ts` reads the same choice from the
 * other side when it puts years on the clock.
 */
export const RADIUS_RATIO = 0.55;

/**
 * Bounds on the box length, in units of its depth.
 *
 * The floor is a domain narrower than it is deep, where a single cell cannot
 * fit; the ceiling is set by the *azimuthal* resolution, which the preset ladder
 * fixes — at L = 8 and `na = 256` a spline element is 0.031 across against 0.008
 * in depth, four to one, and past that the transverse direction is the accuracy
 * bottleneck rather than the radial one the ladder is sized on. 4 is the default
 * because it is the aspect ratio the box benchmarks are usually run at, and
 * because it is close to the annulus it sits beside in the list: at the mid
 * radius that domain is 2π·0.775 around by 0.45 deep.
 */
export const BOX_LENGTH = { min: 1, max: 8, step: 0.5, default: 4 } as const;

/** The `Geometry` a `State` selects. */
export const geometryFor = (s: {
  geometry: GeometryName; boxLength: number;
}): Geometry =>
  GEOMETRY[s.geometry] === "annulus" ? annulus(RADIUS_RATIO, 1) : box(s.boxLength);

/**
 * Labels of the two rheology sliders, named once because two places must agree
 * on them: the pane, and the legend under the equation that tells the reader
 * which slider sets γ and which sets n. Renaming a slider without the legend
 * following would leave it pointing at a control that is not there.
 */
export const LABELS = {
  contrast: "log₁₀ contrast",
  n: "power-law n",
} as const;

export interface State {
  geometry: GeometryName;
  /** Length of the Cartesian box, in units of its depth. Ignored by the annulus. */
  boxLength: number;
  /** log₁₀ Ra — the slider's coordinate, and the one the physics is smooth in. */
  logRa: number;
  dt: number;
  /** Steps per frame; may be < 1. See SPEEDS. */
  speed: number;
  paused: boolean;
  contours: number;
  lineWidth: number;
  mesh: MeshName;
  /** Span of the Nusselt plot, in solver steps; `Infinity` for the whole trace. */
  nuWindow: number;
  wavenumber: number;
  resolution: PresetName;
  viscosity: ViscosityName;
  /** log₁₀ of the viscosity contrast μ(T=0)/μ(T=1); γ = ln of it. */
  logContrast: number;
  /** Krylov iterations per solve. Fixed budget ⇒ predictable frame time. */
  iters: number;
  /** Power-law index. 1 is Newtonian; ≈3 is dislocation creep. */
  n: number;
  /** Rheology updates per solve. Each one costs a full Krylov budget. */
  picard: number;
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
  logRa: Math.log10(2e4),
  dt: PRESETS[DEFAULT_PRESET].dt,
  speed: 2,
  paused: false,
  // Both overlays start off: the temperature field is the subject, and the first
  // thing on screen should be it rather than a lattice drawn over it. Each is a
  // uniform write away.
  contours: 0,
  lineWidth: 1.1,
  mesh: "off",
  // The whole trace by default. The plot's first job is the initial transient
  // settling, which is the one thing a fixed window would cut off — narrowing it
  // is what the reader does once they want to see the settled state resolved.
  nuWindow: NU_WINDOWS.all,
  wavenumber: 4,
  resolution: DEFAULT_PRESET,
  viscosity: "constant",
  logContrast: 3,
  iters: DEFAULT_ITERS,
  n: 3,
  picard: 1,
});
