/**
 * Control state and the tables behind the two list controls.
 *
 * Kept separate from `controls.ts` so the invariants can be regression-tested
 * without a DOM: an `N_φ` that is not a power of two, or a dt that does not fall
 * with the grid, are both mistakes that would otherwise surface only when a user
 * happened to pick that entry.
 */

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

export interface State {
  /** log₁₀ Ra — the slider's coordinate, and the one the physics is smooth in. */
  logRa: number;
  dt: number;
  /** Steps per frame; may be < 1. See SPEEDS. */
  speed: number;
  paused: boolean;
  contours: number;
  lineWidth: number;
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
  logRa: Math.log10(2e4),
  dt: PRESETS[DEFAULT_PRESET].dt,
  speed: 2,
  paused: false,
  contours: 24,
  lineWidth: 1.1,
  wavenumber: 4,
  resolution: DEFAULT_PRESET,
  viscosity: "constant",
  logContrast: 3,
  iters: DEFAULT_ITERS,
  n: 3,
  picard: 1,
});
