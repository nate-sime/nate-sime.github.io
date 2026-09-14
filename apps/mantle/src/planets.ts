/**
 * Planetary examples are parameter profiles for one solver, plus independent
 * presentation data. This first registry entry reproduces the application's
 * existing Earth assumptions exactly; later bodies add entries, not solver
 * branches.
 */

import type { SurfaceMaterialId } from "./gpu/surfaceAssets";
import type { State } from "./ui/presets";

export type PlanetId = "earth" | "venus" | "mars";

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

/**
 * Reduced profile of King (2018)'s reference Venus mantle calculation.
 *
 * The radii, diffusivity, and boundary condition reproduce that published
 * spherical-shell setup. Its published Rayleigh number is far above this
 * teaching grid's resolvable range, so the profile deliberately uses 10^6.5
 * for a legible, resolved transient. This app cannot reproduce its full
 * temperature/depth-dependent viscosity, plastic yielding, internal heating,
 * or spherical-harmonic initial condition, so those are intentionally called
 * out rather than hidden behind a decorative surface.
 */
export const VENUS: PlanetDefinition = {
  id: "venus",
  label: "Venus",
  model: {
    name: "King (2018) reduced Venus mantle profile",
    version: "reference geometry / resolved teaching Ra",
    summary: "A reduced annulus mapping of a published 3-D Venus stagnant-lid calculation.",
    caveats: [
      "The published calculation uses temperature- and depth-dependent viscosity, plastic yielding, internal heating, and 3-D spherical-shell flow; this app maps only its geometry, free-slip boundaries, diffusivity, and Rayleigh number.",
      "The published reference Rayleigh number (3.18 × 10^8) is not resolvable on this interactive grid; this profile uses log₁₀ Ra = 6.5 for a resolved illustrative transient.",
      "Venus' core radius is model-dependent because direct seismic constraints are unavailable; this profile uses 3,110 km from the selected model, not a measured boundary.",
      "The exterior uses a NASA/JPL-Caltech Magellan radar-derived texture. Its colour treatment is not a natural-colour view, and neither the exterior nor atmosphere is a surface-temperature or atmospheric simulation.",
    ],
    sources: [
      { label: "King (2018), Venus resurfacing constrained by geoid and topography", url: "https://doi.org/10.1002/2017JE005475" },
      { label: "NASA/JPL planetary physical parameters", url: "https://ssd.jpl.nasa.gov/planets/phys_par.html" },
    ],
  },
  physical: { surfaceRadiusKm: 6052, mantleBottomRadiusKm: 3110, thermalDiffusivityM2s: 1e-6 },
  solver: {
    state: {
      geometry: "spherical annulus", radialWalls: "free-slip", logRa: 6.5,
      viscosity: "constant", logContrast: 3, logDepthContrast: 0, n: 3, picard: 1,
    },
    initialWavenumber: 1,
  },
  visual: {
    surface: "venus-magellan",
    atmosphere: { color: [0.94, 0.63, 0.18], strength: 0.42 },
    axialTiltDeg: 177.36,
    overviewOrbit: 0.72,
    overviewRadius: 0.95,
  },
};

/**
 * Reduced Mars profile based on the geometry used by Roberts (2006). The
 * chosen Ra is a resolved point inside the active-convection range discussed
 * by Li et al. (2007), rather than an attempt to reproduce a stagnant-lid,
 * internally heated 3-D Mars calculation with this 2-D solver.
 */
export const MARS: PlanetDefinition = {
  id: "mars",
  label: "Mars",
  model: {
    name: "Roberts (2006) reduced Mars mantle profile",
    version: "layered-viscosity geometry / resolved teaching Ra",
    summary: "A reduced annulus mapping of published Mars mantle geometry and active-convection estimates.",
    caveats: [
      "The reference calculations use 3-D spherical shells, internal and basal heating, and depth-dependent material properties; this app maps only a core radius, surface radius, free-slip boundaries, a reference diffusivity, and constant viscosity.",
      "logâ‚â‚€ Ra = 7 is a resolved teaching midpoint of the 2 Ã— 10^6 to 3 Ã— 10^7 active-convection range discussed by Li et al. (2007); it is not a reconstruction of Mars' present thermal state.",
      "The core radius is model-dependent because Mars lacks direct seismic constraints. This profile uses the 1,650 km layered-viscosity case in Roberts (2006), alongside its 3,400 km planetary radius.",
      "The exterior is a Viking-image-derived NASA/JPL-Caltech map for orientation only; it does not represent topography, albedo physics, surface temperature, or the atmosphere.",
    ],
    sources: [
      { label: "Roberts (2006), Geoid and topography of Mars from a dynamic mantle", url: "https://doi.org/10.1029/2005JE002668" },
      { label: "Li et al. (2007), Could thermal mantle plumes have caused magnetic anomalies on Mars?", url: "https://doi.org/10.1029/2007GL030544" },
      { label: "NASA 3D Resources, Mars image texture", url: "https://science.nasa.gov/3d-resources/mars/" },
    ],
  },
  physical: { surfaceRadiusKm: 3400, mantleBottomRadiusKm: 1650, thermalDiffusivityM2s: 2.132e-6 },
  solver: {
    state: {
      geometry: "spherical annulus", radialWalls: "free-slip", logRa: 7,
      viscosity: "constant", logContrast: 3, logDepthContrast: 0, n: 3, picard: 1,
    },
    initialWavenumber: 1,
  },
  visual: {
    surface: "mars-viking",
    atmosphere: { color: [0.66, 0.25, 0.10], strength: 0.08 },
    axialTiltDeg: 25.19,
    overviewOrbit: 1.52,
    overviewRadius: 0.53,
  },
};

export const PLANETS: Record<PlanetId, PlanetDefinition> = { earth: EARTH, venus: VENUS, mars: MARS };

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
