/**
 * Dimensional time, for a simulation that is otherwise entirely nondimensional.
 *
 * The solver has no physical units in it anywhere and should not acquire any: the
 * whole point of the Boussinesq scaling is that Ra is the only parameter, and one
 * run at Ra = 2×10⁴ stands for every combination of gravity, viscosity, layer
 * depth and thermal expansivity that produces it. Attaching seconds to the clock
 * is therefore a *display* choice, and it is one that needs its assumption stated
 * rather than buried — which is why it lives in one file, is exported as data,
 * and is printed on screen beside the numbers it produces.
 *
 * The governing equations are written with lengths scaled by whichever length is
 * 1 in code units and temperature diffusion of unit strength,
 *
 *     ∂T/∂t + u·∇T = ∇²T,
 *
 * so one unit of nondimensional time is one thermal diffusion time across that
 * length, ℓ²/κ. **Which length that is differs between the two geometries**, and
 * the difference is a factor of five in the clock, so it is not a detail:
 *
 *   annulus — `r_o = 1`, the outer radius. ℓ is Earth's surface radius, and the
 *             radius ratio r_i/r_o = 0.55 the geometry already uses is the
 *             core–mantle boundary against it, 3480 km / 6371 km = 0.546.
 *   box     — depth `= 1`, the convecting layer. ℓ is the mantle's thickness,
 *             6371 − 3480 = 2891 km, which is the same physical layer measured
 *             the way a box measures it.
 *
 * Nothing else in the model fixes ℓ or κ, so the reference below does.
 *
 * **Expect very large numbers, and read them as a result rather than a fault.**
 * R_o²/κ is about 1.3×10¹² years — thermal diffusion across the mantle is some
 * two orders of magnitude slower than the age of the Earth, which is exactly why
 * the real mantle convects instead of conducting. A run that settles by t ≈ 0.1
 * has therefore run for ~130 Gyr of diffusion time, and that is not a
 * contradiction: at Ra = 2×10⁴ this model is thousands of times more viscous than
 * the mantle, whose Ra is 10⁷–10⁹, so it needs many diffusion times to do what
 * the Earth does in a fraction of one. The dimensional clock is honest about the
 * scaling; it is not a claim that the simulation is Earth.
 */

import type { GeometryKind } from "../geometry";

/**
 * The scales that turn nondimensional time into seconds. Exported so the readout
 * can print them: a reader given "133 Gyr" and no reference has been told a
 * number, not a quantity.
 */
export const REFERENCE = {
  /** Outer radius in metres — Earth's surface, the length the annulus' r_o = 1 is. */
  Ro: 6.371e6,
  /** Layer thickness in metres — the mantle, the length the box's depth 1 is. */
  depth: 2.891e6,
  /** Thermal diffusivity in m²/s — a standard mantle value. */
  kappa: 1e-6,
} as const;

/** Seconds in a Julian year, the definition the "yr" below means. */
export const YEAR = 3.15576e7;

/** The length that is 1 in code units, in metres. */
export const lengthScale = (kind: GeometryKind): number =>
  kind === "annulus" ? REFERENCE.Ro : REFERENCE.depth;

/** Seconds per unit of nondimensional time: ℓ²/κ. */
export const timeUnit = (kind: GeometryKind): number =>
  lengthScale(kind) ** 2 / REFERENCE.kappa;

/** Years per unit of nondimensional time — 1.3×10¹² on the annulus, 2.6×10¹¹ in a box. */
export const timeUnitYears = (kind: GeometryKind): number =>
  timeUnit(kind) / YEAR;

/**
 * Ladder of year multiples, largest first. It runs past Gyr because it has to:
 * one nondimensional time unit is already 1.3 Tyr, so a run of order t = 1 leaves
 * the geological units behind entirely, and clamping at Gyr would print four
 * digits where a unit change says it better.
 */
const LADDER: readonly (readonly [number, string])[] = [
  [1e12, "Tyr"], [1e9, "Gyr"], [1e6, "Myr"], [1e3, "kyr"], [1, "yr"],
];

/** Three significant figures, without an exponent inside the ladder's range. */
const sig3 = (v: number): string =>
  v < 10 ? v.toFixed(2) : v < 100 ? v.toFixed(1) : v.toFixed(0);

/**
 * A nondimensional time as a dimensional one, with the unit that keeps it
 * readable. `—` for a clock that has not started, matching the readout's
 * convention for a diagnostic that is not there yet.
 */
export function dimensionalTime(kind: GeometryKind, t: number): string {
  if (!Number.isFinite(t)) return "—";
  const yr = (t * timeUnit(kind)) / YEAR;
  const a = Math.abs(yr);
  if (a === 0) return "0 yr";
  // Below a year the ladder has nothing left to divide by; an exponent is the
  // honest fallback. Unreachable at these scales, and cheaper than pretending.
  if (a < 1) return `${yr.toExponential(2)} yr`;
  const [scale, unit] = LADDER.find(([s]) => a >= s) ?? [1, "yr"];
  return `${sig3(yr / scale)} ${unit}`;
}

/**
 * The scaling itself, for the readout's header — the assumption, stated. The
 * symbol names the length that is 1 in code units, so a reader can see at a
 * glance that the box and the annulus are not on the same clock.
 */
export const referenceNote = (kind: GeometryKind): string => {
  const sym = kind === "annulus" ? "R_o" : "d";
  return `t in ${sym}²/κ = ${timeUnitYears(kind).toExponential(2)} yr` +
    // κ in exponential form: the default stringification of 1e-6 is "0.000001",
    // six leading zeros a reader has to count to recognise a number they know.
    `   (${sym} = ${(lengthScale(kind) / 1e3).toFixed(0)} km, ` +
    `κ = ${REFERENCE.kappa.toExponential(0)} m²/s)`;
};
