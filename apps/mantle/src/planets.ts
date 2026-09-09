/**
 * Planetary examples are parameter profiles for one solver, plus independent
 * presentation data. This first registry entry reproduces the application's
 * existing Earth assumptions exactly; later bodies add entries, not solver
 * branches.
 */

import type { SurfaceMaterialId } from "./gpu/surfaceAssets";
import type { State } from "./ui/presets";

export type PlanetId = "earth";

export interface PlanetDefinition {
  readonly id: PlanetId;
  readonly label: string;
  readonly model: {
    readonly name: string;
    readonly version: string;
    readonly summary: string;
    readonly caveats: readonly string[];
    readonly sources: readonly { label: string; url: string }[];
  };
  readonly physical: {
    readonly surfaceRadiusKm: number;
    readonly mantleBottomRadiusKm: number;
    readonly thermalDiffusivityM2s: number;
  };
  readonly solver: {
    readonly state: Pick<State,
      "geometry" | "radialWalls" | "logRa" | "viscosity" |
      "logContrast" | "logDepthContrast" | "n" | "picard">;
    readonly initialWavenumber: number;
  };
  readonly visual: {
    readonly surface: SurfaceMaterialId;
    readonly atmosphere?: { readonly color: readonly [number, number, number]; readonly strength: number };
    readonly axialTiltDeg: number;
    readonly overviewOrbit: number;
    readonly overviewRadius: number;
  };
}

export interface PlanetRadii {
  /** Mantle thickness, the physical length represented by one code unit. */
  readonly depthKm: number;
  /** CMB radius in mantle-thickness units. */
  readonly ri: number;
  /** Surface radius in mantle-thickness units; always `ri + 1`. */
  readonly ro: number;
}

export const radiiFor = (planet: PlanetDefinition): PlanetRadii => {
  const depthKm = planet.physical.surfaceRadiusKm - planet.physical.mantleBottomRadiusKm;
  if (!(depthKm > 0) || !(planet.physical.thermalDiffusivityM2s > 0))
    throw new Error(`Invalid physical scale for ${planet.id}.`);
  // The existing app stored the Earth ratio to nine decimal places. Keeping
  // that derived precision makes profile extraction behavior-preserving while
  // still deriving it from the source radii rather than carrying a second value.
  const ri = Math.round((planet.physical.mantleBottomRadiusKm / depthKm) * 1e9) / 1e9;
  return { depthKm, ri, ro: ri + 1 };
};

/** Earth's previous default: 3486 km CMB, 6371 km surface, 1e-6 m²/s κ. */
export const EARTH: PlanetDefinition = {
  id: "earth",
  label: "Earth",
  model: {
    name: "Earth reference mantle",
    version: "PREM geometry / current app reference",
    summary: "The existing Earth-like annulus and display scale, expressed as data.",
    caveats: [
      "This is a 2-D Boussinesq teaching model, not a full Earth simulation.",
      "The thermal diffusivity is a display reference; it does not alter the nondimensional solver.",
    ],
    sources: [
      { label: "Dziewonski & Anderson (1981), Preliminary Reference Earth Model", url: "https://doi.org/10.1016/0031-9201(81)90046-7" },
      { label: "NASA Visible Earth, Blue Marble", url: "https://visibleearth.nasa.gov/images/57752" },
    ],
  },
  physical: { surfaceRadiusKm: 6371, mantleBottomRadiusKm: 3486, thermalDiffusivityM2s: 1e-6 },
  solver: {
    state: {
      geometry: "spherical annulus", radialWalls: "free-slip", logRa: Math.log10(1e6),
      viscosity: "constant", logContrast: 3, logDepthContrast: 0, n: 3, picard: 1,
    },
    initialWavenumber: 4,
  },
  visual: {
    surface: "earth-daymap", axialTiltDeg: 0, overviewOrbit: 1, overviewRadius: 1,
  },
};

export const PLANETS: Record<PlanetId, PlanetDefinition> = { earth: EARTH };

export const planetFor = (id: PlanetId | null | undefined): PlanetDefinition => PLANETS[id ?? "earth"];

/** The small physics subset a named profile owns; presentation controls stay out. */
export type PlanetOwnedState = Pick<State,
  "activePlanet" | "geometry" | "radialWalls" | "logRa" | "viscosity" |
  "logContrast" | "logDepthContrast" | "n" | "picard" | "wavenumber">;

/**
 * Derive the modified badge instead of storing a flag that can become stale.
 * A box benchmark and a missing profile are necessarily custom configurations.
 */
export const isPlanetProfileModified = (state: PlanetOwnedState): boolean => {
  if (state.activePlanet === null) return true;
  const planet = planetFor(state.activePlanet);
  const expected = planet.solver.state;
  return state.geometry !== expected.geometry ||
    state.radialWalls !== expected.radialWalls ||
    state.logRa !== expected.logRa ||
    state.viscosity !== expected.viscosity ||
    state.logContrast !== expected.logContrast ||
    state.logDepthContrast !== expected.logDepthContrast ||
    state.n !== expected.n ||
    state.picard !== expected.picard ||
    state.wavenumber !== planet.solver.initialWavenumber;
};
