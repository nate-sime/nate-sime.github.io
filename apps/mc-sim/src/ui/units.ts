// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The one place the beam — and the plate — acquire physical units.
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
 *
 * A plate is the same with L the side along y, the bending stiffness
 * D = Et³/12(1 − ν²) for EI, the mass per area ρt for ρA, and a pressure q₀:
 * W = q₀L⁴/D or P₀L²/D, Ω = (D / (ρt L⁴))^½, ℓ(w) = W·(q₀L² or P₀).
 */

import type { LoadCase, QoI } from "../hierarchy";

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

export interface ReferencePlate {
  /** Side along y, m; the side along x is the aspect ratio times it. */
  readonly L: number;
  readonly E: number;
  readonly nu: number;
  /** Thickness, m. */
  readonly t: number;
  readonly rho: number;
  /** Pressure (Pa) and point force (N). */
  readonly q0: number;
  readonly P0: number;
}

export type Reference = ReferenceBeam | ReferencePlate;

const isPlate = (r: Reference): r is ReferencePlate => "t" in r;

/** A 1 m steel bar, 20 mm square — a ruler you could hold and pluck. */
export const STEEL_BAR: ReferenceBeam = { L: 1, E: 210e9, b: 0.02, h: 0.02, rho: 7850, q0: 100, P0: 100 };

/** A 1 m square steel plate, 5 mm thick, under a 1 kPa pressure. */
export const STEEL_PLATE: ReferencePlate = { L: 1, E: 210e9, nu: 0.3, t: 0.005, rho: 7850, q0: 1000, P0: 100 };

export const flexuralRigidity = (r: ReferenceBeam) => (r.E * r.b * r.h ** 3) / 12;
export const massPerLength = (r: ReferenceBeam) => r.rho * r.b * r.h;
export const plateRigidity = (r: ReferencePlate) => (r.E * r.t ** 3) / (12 * (1 - r.nu ** 2));
export const massPerArea = (r: ReferencePlate) => r.rho * r.t;

/** The stiffness and mass the nondimensional ones are relative to: EI and ρA, or D and ρt. */
const stiffnessOf = (r: Reference) => (isPlate(r) ? plateRigidity(r) : flexuralRigidity(r));
const massOf = (r: Reference) => (isPlate(r) ? massPerArea(r) : massPerLength(r));

/** rad/s per unit ω̂. */
export const frequencyScale = (r: Reference) => Math.sqrt(stiffnessOf(r) / (massOf(r) * r.L ** 4));

/** m per unit ŵ. */
export const deflectionScale = (r: Reference, load: LoadCase) =>
  load === "uniform" ? (r.q0 * r.L ** 4) / stiffnessOf(r) : (r.P0 * r.L ** (isPlate(r) ? 2 : 3)) / stiffnessOf(r);

/** J per unit ℓ̂. */
export const complianceScale = (r: Reference, load: LoadCase) =>
  deflectionScale(r, load) * (load === "uniform" ? r.q0 * r.L ** (isPlate(r) ? 2 : 1) : r.P0);

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

/** How results are shown: the toggle, and the beam or plate the dimensional view assumes. */
export interface Display {
  readonly dimensional: boolean;
  readonly ref: Reference;
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
export function referenceNote(r: Reference): string {
  if (isPlate(r))
    return `reference: side L = ${si(r.L, "m", 3)}, thickness ${si(r.t, "m", 3)}, E = ${si(r.E, "Pa", 3)}, ν = ${r.nu}, ` +
      `ρ = ${r.rho} kg/m³  →  D = ${si(plateRigidity(r), "N·m", 3)}, ρt = ${massPerArea(r).toPrecision(3)} kg/m²; ` +
      `q₀ = ${si(r.q0, "Pa", 3)}, P₀ = ${si(r.P0, "N", 3)}`;
  return `reference: L = ${si(r.L, "m", 3)}, ${si(r.b, "m", 3)} × ${si(r.h, "m", 3)}, ` +
    `E = ${si(r.E, "Pa", 3)}, ρ = ${r.rho} kg/m³  →  EI = ${si(flexuralRigidity(r), "N·m²", 3)}, ` +
    `ρA = ${massPerLength(r).toPrecision(3)} kg/m; q₀ = ${si(r.q0, "N/m", 3)}, P₀ = ${si(r.P0, "N", 3)}`;
}

/** A prefix for an axis whose largest magnitude is `maxAbs` (in base units): the factor to multiply by, and the unit. */
export function axisUnit(maxAbs: number, unit: string): { factor: number; label: string } {
  const [s, pre] = PREFIX.find(([s]) => maxAbs >= s) ?? PREFIX[PREFIX.length - 1];
  return { factor: 1 / s, label: `${pre}${unit}` };
}

/**
 * How a quantity of interest is drawn on an axis: the factor from its
 * nondimensional value, and the unit — hertz for ω₁ (as f = ω/2π), metres for
 * a deflection or the RMS field, joules for compliance; bare when nondimensional.
 */
export function qoiScale(d: Display, qoi: QoI, load: LoadCase): { factor: number; unit: string } {
  if (!d.dimensional) return { factor: 1, unit: "" };
  if (qoi === "omega1") return { factor: frequencyScale(d.ref) / (2 * Math.PI), unit: "Hz" };
  if (qoi === "compliance") return { factor: complianceScale(d.ref, load), unit: "J" };
  return { factor: deflectionScale(d.ref, load), unit: "m" };
}
