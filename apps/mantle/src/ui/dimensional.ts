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
 * length, ℓ²/κ. **Both geometries use the same ℓ**: the mantle's own thickness.
 * The box has depth 1 by construction; the annulus is nondimensionalised the
 * same way, with `r_o − r_i = 1` (see `geometry.ts`), so a code-unit length is
 * already a multiple of that one physical length in either domain. There is
 * therefore one clock, not two.
 *
 * Nothing else in the model fixes ℓ or κ, so the reference below does.
 *
 * **Expect very large numbers, and read them as a result rather than a fault.**
 * d²/κ is about 2.6×10¹¹ years — thermal diffusion across the mantle is some two
 * orders of magnitude slower than the age of the Earth, which is exactly why the
 * real mantle convects instead of conducting. A run that settles by t ≈ 0.1 has
 * therefore run for ~26 Gyr of diffusion time, and that is not a contradiction:
 * at Ra = 2×10⁴ this model is thousands of times more viscous than the mantle,
 * whose Ra is 10⁷–10⁹, so it needs many diffusion times to do what the Earth
 * does in a fraction of one. The dimensional clock is honest about the scaling;
 * it is not a claim that the simulation is Earth.
 */

/**
 * The scale that turns nondimensional time into seconds. Exported so the readout
 * can print it: a reader given "133 Gyr" and no reference has been told a
 * number, not a quantity.
 */
export const REFERENCE = {
  /** Mantle thickness in metres — the length that is 1 in code units. */
  depth: 2.885e6,
  /** Thermal diffusivity in m²/s — a standard mantle value. */
  kappa: 1e-6,
} as const;

/** Seconds in a Julian year, the definition the "yr" below means. */
export const YEAR = 3.15576e7;

/** The length that is 1 in code units, in metres. */
export const lengthScale = REFERENCE.depth;

/**
 * The same length in km — the unit a mantle depth is actually read in.
 * Dividing a dimensional depth by this gives the nondimensional fraction of
 * the layer's thickness the solver wants; multiplying goes the other way. See
 * `ui/controls.ts`'s Brandenburg A₀-transition-depth slider for the one place
 * that round-trip currently happens.
 */
export const MANTLE_THICKNESS_KM = lengthScale / 1e3;

/** Seconds per unit of nondimensional time: ℓ²/κ. */
export const timeUnit = lengthScale ** 2 / REFERENCE.kappa;

/** Years per unit of nondimensional time. */
export const timeUnitYears = timeUnit / YEAR;

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
  const yr = (t * timeUnit) / YEAR;
  const a = Math.abs(yr);
  if (a === 0) return "0 yr";
  // Below a year the ladder has nothing left to divide by; an exponent is the
  // honest fallback. Unreachable at these scales, and cheaper than pretending.
  if (a < 1) return `${yr.toExponential(2)} yr`;
  const [scale, unit] = LADDER.find(([s]) => a >= s) ?? [1, "yr"];
  return `${sig3(yr / scale)} ${unit}`;
}

/**
 * The scaling itself, for the readout's header — the assumption, stated. `d`
 * names the length that is 1 in code units: the mantle's own thickness, the
 * same in both geometries.
 */
export const referenceNote = (): string =>
  `t in d²/κ = ${timeUnitYears.toExponential(2)} yr` +
  // κ in exponential form: the default stringification of 1e-6 is "0.000001",
  // six leading zeros a reader has to count to recognise a number they know.
  `   (d = ${(lengthScale / 1e3).toFixed(0)} km, ` +
  `κ = ${REFERENCE.kappa.toExponential(0)} m²/s)`;
