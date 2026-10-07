// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The one place the beam acquires physical units.
 *
 * The solver works in the span L, the reference stiffness EI and the reference
 * mass per length ρA, and never sees a metre. Every result therefore converts
 * with one scale, and the conversion is a display choice made here:
 *
 *   x  = L x̂
 *   w  = W ŵ,   W = q₀L⁴/EI  (distributed load q₀)   or   P₀L³/EI  (point load P₀)
 *   ω  = Ω ω̂,   Ω = (EI / (ρA L⁴))^½,   f = ω / 2π
 *   ℓ(w) = W·(q₀L or P₀) ℓ̂   — the work done, in joules
 *
 * so a toggle between the two never re-solves anything. Relative errors are the
 * same number in both, which is the point of plotting them.
 */

import type { LoadCase } from "../hierarchy";

export interface ReferenceBeam {
  /** Span, m. */
  readonly L: number;
  /** Young's modulus, Pa. */
  readonly E: number;
  /** Section width and depth, m (rectangular, at the root for a tapered beam). */
  readonly b: number;
  readonly h: number;
  /** Density, kg/m³. */
  readonly rho: number;
  /** Load magnitudes: distributed (N/m) and point (N). */
  readonly q0: number;
  readonly P0: number;
}

/** A 1 m steel bar, 20 mm square — a ruler you could hold and pluck. */
export const STEEL_BAR: ReferenceBeam = { L: 1, E: 210e9, b: 0.02, h: 0.02, rho: 7850, q0: 100, P0: 100 };

export const flexuralRigidity = (r: ReferenceBeam) => (r.E * r.b * r.h ** 3) / 12;
export const massPerLength = (r: ReferenceBeam) => r.rho * r.b * r.h;

/** rad/s per unit ω̂. */
export const frequencyScale = (r: ReferenceBeam) => Math.sqrt(flexuralRigidity(r) / (massPerLength(r) * r.L ** 4));

/** m per unit ŵ. */
export const deflectionScale = (r: ReferenceBeam, load: LoadCase) =>
  load === "uniform" ? (r.q0 * r.L ** 4) / flexuralRigidity(r) : (r.P0 * r.L ** 3) / flexuralRigidity(r);

/** J per unit ℓ̂. */
export const complianceScale = (r: ReferenceBeam, load: LoadCase) =>
  deflectionScale(r, load) * (load === "uniform" ? r.q0 * r.L : r.P0);

const PREFIX: readonly [number, string][] = [
  [1e9, "G"], [1e6, "M"], [1e3, "k"], [1, ""], [1e-3, "m"], [1e-6, "µ"], [1e-9, "n"], [1e-12, "p"],
];

/** A value in SI with an engineering prefix: 0.00453 m → "4.53 mm". */
export function si(v: number, unit: string, digits = 4): string {
  if (!Number.isFinite(v)) return "—";
  if (v === 0) return `0 ${unit}`;
  const a = Math.abs(v);
  const [s, pre] = PREFIX.find(([s]) => a >= s * 0.9995) ?? PREFIX[PREFIX.length - 1];
  // Fixed decimals, not toPrecision: 100 at two digits must read "100", not "1.0e+2".
  const m = v / s, decimals = Math.max(0, digits - 1 - Math.floor(Math.log10(Math.abs(m)) + 1e-9));
  return `${m.toFixed(decimals)} ${pre}${unit}`;
}

/** A nondimensional value: plain when moderate, exponent when not. */
export function plain(v: number, digits = 6): string {
  if (!Number.isFinite(v)) return "—";
  const a = Math.abs(v);
  return a !== 0 && (a < 1e-3 || a >= 1e5) ? v.toExponential(digits - 1) : v.toPrecision(digits);
}

/** How results are shown: the toggle, and the beam the dimensional view assumes. */
export interface Display {
  readonly dimensional: boolean;
  readonly ref: ReferenceBeam;
}

export const fmt = {
  length: (d: Display, x: number) => (d.dimensional ? si(x * d.ref.L, "m") : plain(x, 4)),
  deflection: (d: Display, w: number, load: LoadCase) =>
    d.dimensional ? si(w * deflectionScale(d.ref, load), "m") : plain(w),
  /** ω̂ ↦ f in Hz, or ω̂ itself. */
  frequency: (d: Display, omega: number) =>
    d.dimensional ? si((omega * frequencyScale(d.ref)) / (2 * Math.PI), "Hz") : plain(omega),
  compliance: (d: Display, c: number, load: LoadCase) =>
    d.dimensional ? si(c * complianceScale(d.ref, load), "J") : plain(c),
};

/** The assumption behind every dimensional number, printed beside them. */
export function referenceNote(r: ReferenceBeam): string {
  return `reference: L = ${si(r.L, "m", 3)}, ${si(r.b, "m", 3)} × ${si(r.h, "m", 3)}, ` +
    `E = ${si(r.E, "Pa", 3)}, ρ = ${r.rho} kg/m³  →  EI = ${si(flexuralRigidity(r), "N·m²", 3)}, ` +
    `ρA = ${massPerLength(r).toPrecision(3)} kg/m; q₀ = ${si(r.q0, "N/m", 3)}, P₀ = ${si(r.P0, "N", 3)}`;
}

/** A prefix for an axis whose largest magnitude is `maxAbs` (in base units): the factor to multiply by, and the unit. */
export function axisUnit(maxAbs: number, unit: string): { factor: number; label: string } {
  const [s, pre] = PREFIX.find(([s]) => maxAbs >= s) ?? PREFIX[PREFIX.length - 1];
  return { factor: 1 / s, label: `${pre}${unit}` };
}
