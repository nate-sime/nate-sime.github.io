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
 * The governing equations are written with lengths scaled by the outer radius
 * (r_o = 1 in code units) and temperature diffusion of unit strength,
 *
 *     ∂T/∂t + u·∇T = ∇²T,
 *
 * so one unit of nondimensional time is one thermal diffusion time across the
 * outer radius, R_o²/κ. Nothing else in the model fixes R_o or κ, so the
 * reference below does — at Earth's mantle, which the geometry was already chosen
 * to match: the radius ratio r_i/r_o = 0.55 is the core–mantle boundary against
 * the surface, 3480 km / 6371 km = 0.546.
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

/**
 * The scales that turn nondimensional time into seconds. Exported so the readout
 * can print them: a reader given "133 Gyr" and no reference has been told a
 * number, not a quantity.
 */
export const REFERENCE = {
  /** Outer radius in metres — Earth's surface, the length the code's r_o = 1 is. */
  Ro: 6.371e6,
  /** Thermal diffusivity in m²/s — a standard mantle value. */
  kappa: 1e-6,
} as const;

/** Seconds in a Julian year, the definition the "yr" below means. */
export const YEAR = 3.15576e7;

/** Seconds per unit of nondimensional time: R_o²/κ. */
export const TIME_UNIT = REFERENCE.Ro ** 2 / REFERENCE.kappa;

/** Years per unit of nondimensional time — about 1.3×10¹². */
export const TIME_UNIT_YEARS = TIME_UNIT / YEAR;

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
export function dimensionalTime(t: number): string {
  if (!Number.isFinite(t)) return "—";
  const yr = (t * TIME_UNIT) / YEAR;
  const a = Math.abs(yr);
  if (a === 0) return "0 yr";
  // Below a year the ladder has nothing left to divide by; an exponent is the
  // honest fallback. Unreachable at these scales, and cheaper than pretending.
  if (a < 1) return `${yr.toExponential(2)} yr`;
  const [scale, unit] = LADDER.find(([s]) => a >= s) ?? [1, "yr"];
  return `${sig3(yr / scale)} ${unit}`;
}

/** The scaling itself, for the readout's header — the assumption, stated. */
export const referenceNote = (): string =>
  `t in R_o²/κ = ${TIME_UNIT_YEARS.toExponential(2)} yr` +
  // κ in exponential form: the default stringification of 1e-6 is "0.000001",
  // six leading zeros a reader has to count to recognise a number they know.
  `   (R_o = ${(REFERENCE.Ro / 1e3).toFixed(0)} km, ` +
  `κ = ${REFERENCE.kappa.toExponential(0)} m²/s)`;
