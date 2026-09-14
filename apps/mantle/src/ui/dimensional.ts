/**
 * Dimensional display scales for the otherwise nondimensional solver.
 *
 * A scale is deliberately passed to formatting functions rather than stored in
 * solver state. A profile can change what code-unit length means without
 * teaching the numerical kernels about metres or seconds.
 */

export interface DimensionalScale {
  /** Physical length represented by one code-unit depth, in metres. */
  readonly depth: number;
  /** Thermal diffusivity in m²/s. */
  readonly kappa: number;
}

export const YEAR = 3.15576e7;

export const createDimensionalScale = (depth: number, kappa: number): DimensionalScale => {
  if (!(depth > 0) || !(kappa > 0)) throw new Error("Dimensional scales require positive depth and diffusivity.");
  return { depth, kappa };
};

/** The legacy Earth reference, retained as the default for box benchmarks. */
export const REFERENCE = createDimensionalScale(2.885e6, 1e-6);

/** Seconds per unit of nondimensional time: ℓ²/κ. */
export function timeUnitFor(scale: DimensionalScale): number {
  return scale.depth ** 2 / scale.kappa;
}

/** Years per unit of nondimensional time. */
export function timeUnitYearsFor(scale: DimensionalScale): number {
  return timeUnitFor(scale) / YEAR;
}

/** Metres per second per unit of nondimensional velocity: κ/ℓ. */
export function velocityUnitFor(scale: DimensionalScale): number {
  return scale.kappa / scale.depth;
}

/** cm/yr per unit of nondimensional velocity. */
export function velocityUnitCmPerYearFor(scale: DimensionalScale): number {
  return velocityUnitFor(scale) * 100 * YEAR;
}

/** Compatibility aliases for existing consumers and the Earth default. */
export const lengthScale = REFERENCE.depth;
export const MANTLE_THICKNESS_KM = lengthScale / 1e3;
export const timeUnit = timeUnitFor(REFERENCE);
export const timeUnitYears = timeUnitYearsFor(REFERENCE);
export const velocityUnit = velocityUnitFor(REFERENCE);
export const velocityUnitCmPerYear = velocityUnitCmPerYearFor(REFERENCE);

const LADDER: readonly (readonly [number, string])[] = [
  [1e12, "Tyr"], [1e9, "Gyr"], [1e6, "Myr"], [1e3, "kyr"], [1, "yr"],
];

const sig3 = (v: number): string =>
  v < 10 ? v.toFixed(2) : v < 100 ? v.toFixed(1) : v.toFixed(0);

/** A nondimensional time as a readable dimensional value. */
export function dimensionalTime(t: number, scale: DimensionalScale = REFERENCE): string {
  if (!Number.isFinite(t)) return "—";
  const yr = t * timeUnitYearsFor(scale);
  const a = Math.abs(yr);
  if (a === 0) return "0 yr";
  if (a < 1) return `${yr.toExponential(2)} yr`;
  const [unitScale, unit] = LADDER.find(([s]) => a >= s) ?? [1, "yr"];
  return `${sig3(yr / unitScale)} ${unit}`;
}

/** State the display assumption beside dimensional readouts. */
export function referenceNote(scale: DimensionalScale = REFERENCE): string {
  return `t in d²/κ = ${timeUnitYearsFor(scale).toExponential(2)} yr` +
    `   (d = ${(scale.depth / 1e3).toFixed(0)} km, κ = ${scale.kappa.toExponential(0)} m²/s)`;
}

/** A nondimensional velocity as cm/yr. */
export function dimensionalVelocity(v: number, scale: DimensionalScale = REFERENCE): string {
  if (!Number.isFinite(v)) return "—";
  return `${(v * velocityUnitCmPerYearFor(scale)).toFixed(2)} cm/yr`;
}
