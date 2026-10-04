// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The guided tours, as data.
 *
 * No DOM and no GPU here, for exactly the reason `presets.ts` has none: this
 * is a table `tour.ts` renders and `main.ts` drives, and every invariant it
 * relies on — that a target is one the pane actually exposes, that a preset
 * name exists, that a focus point is inside the domain the step is running in
 * — fails at *step* time, part-way through a tour most readers will click
 * through exactly once. `tests/tour.test.ts` checks them without a browser,
 * the same trade `tests/presets.test.ts` makes for `QUICK_STARTS`.
 *
 * A tour is a sequence of small, declarative changes rather than a script of
 * imperative calls: each step says which control it is about, what it wants
 * the model set to, and how long to wait before moving on. The driver decides
 * *how* to make each of those true, which is what lets a step be reordered,
 * skipped or stepped back into without the table knowing anything about the
 * pane, the camera or the solver.
 */

import { EARTH, MARS, PLANETS, radiiFor, type PlanetId } from "../planets";
import {
  DEFAULT_PRESET, PRESETS, type BenchmarkName, type QuickStartName, type State,
} from "./presets";

/**
 * Every element a step may point at. The first group are Tweakpane blades,
 * handed out by `buildPane` (`controls.ts`) — Tweakpane renders no IDs of its
 * own, so a blade's `element` is the only stable way to find one, and the
 * labels are no help: `vigour`'s flips between "convective vigour" and
 * "log₁₀ Ra" with the advanced toggle. The second group are the static
 * containers in index.html, which only `main.ts` can resolve.
 *
 * Declared as a list rather than a union type so `tests/tour.test.ts` can
 * check membership at run time, while `Record<TourTargetName, HTMLElement>`
 * on the resolver still makes the compiler insist every name is wired.
 */
export const TOUR_TARGETS = [
  // pane blades (ui/controls.ts)
  "tutorials", "planet", "preset", "vigour", "seedPattern", "seed", "rock", "paused", "speed", "flow", "tracers",
  "colormap", "restart", "resetView", "view3d", "advanced",
  // the "?" on the domain folder's title, standing in for all six
  "domainHelp",
  // advanced pane blades, one row per folder — what `SECTION_HELP` points at
  "geometry", "boxWidth", "walls", "radialWalls", "resolution",
  "isothermal", "courant", "dtMax", "dtInitial",
  "law", "equation", "contrast", "depthContrast", "powerLawN", "cgIterations",
  "picardSweeps", "yieldStress", "yieldGradient", "etaStar", "etaLight", "etaDense",
  "seedMode",
  "streamlineDensity", "meshOverlay", "lineWidth", "plotWindow",
  "tracerOverlay", "tracerCount", "tracerColour", "tracerSize", "tracerOpacity",
  "composition", "logRb", "layerDepth", "reseedTracers",
  // static containers (index.html)
  "canvas", "traces", "caption",
] as const;

export type TourTargetName = (typeof TOUR_TARGETS)[number];

/** A world-space point to fly the 2-D camera to, or the way back out. */
export type TourFocus = { zoom: number; x: number; y: number; ms: number } | "reset";

/**
 * How long to hold a step before offering the next one automatically. Solver
 * steps rather than wall-clock wherever an *effect* is what is being waited
 * for: the same number of steps is the same amount of physics on a fast
 * machine and a slow one, where the same number of seconds is not.
 * Wall-clock is for steps that are only asking the reader to look at
 * something already on screen.
 */
export type TourDwell = { steps: number } | { ms: number };

/** A small table in a card's body, for numbers a reader compares across rows. */
export interface TourTable {
  head: readonly string[];
  rows: readonly (readonly string[])[];
}

/** A bulleted list in a card's body, for points a reader takes away one at a time. */
export interface TourList {
  items: readonly string[];
}

export interface TourStep {
  /** Stable across edits to the prose — it is what a step is identified by in tests and in a log. */
  id: string;
  title: string;
  /** One paragraph, table or list per entry. */
  body: readonly (string | TourTable | TourList)[];
  /**
   * The control this step is about, lit while everything else is dimmed.
   * `null` dims the whole screen — the opening and closing cards are about
   * the app, not about any one control.
   */
  target: TourTargetName | null;
  /** An additional control to light with `target` when the experiment needs both. */
  highlight?: TourTargetName;
  /**
   * A second chrome panel left at full brightness and reachable, without a
   * spotlight of its own. `highlight` cannot do this when the two live in
   * different panels — the shades are tiled inside one host (see `targetOf`
   * in `tour.ts`) — and an experiment worked on the pane is usually read off
   * the corner traces.
   */
  companion?: TourTargetName;
  /**
   * A fine Rayleigh-number control on the card: the live value, typeable,
   * over a slider spanning only `min`–`max` in log₁₀ Ra. The simple view's
   * vigour slider carries no number and spans nine decades in one short
   * track, and a step that asks the reader to bracket a threshold has to
   * show where they are and let them move by a percent or two. Kept in step
   * with the pane's slider both ways — it writes through `setLogRa`, and
   * re-reads `State` every frame. Typed values may go outside the range;
   * the card's slider then rests at its end.
   */
  raControl?: { min: number; max: number };
  /**
   * Reconciled against the live `Globe3D.viewMode`, and only toggled when the
   * two differ — so stepping backwards through a tour doesn't flip the view
   * on a step that never asked to change it.
   */
  view?: "2d" | "3d";
  /**
   * Show the pane's advanced view for this step. Unlike `view`, applied on
   * every step of a guided tour, absent meaning the simple view: most steps
   * point at simple controls the advanced view hides ("how the rock
   * behaves", "show flow lines"), so a step stepped back into from an
   * advanced one has to get its own view back without having to say so.
   * Never applied by a section's help, which is opened from the advanced
   * view and is about it.
   */
  advanced?: boolean;
  /** Applied through `PaneHandle.applyPatch` — literally the path a preset selection takes. */
  patch?: Partial<State>;
  /** Same, by name, for the entries already in `QUICK_STARTS`/`BENCHMARKS`. */
  preset?: QuickStartName | BenchmarkName;
  /**
   * The planet this step runs on. Selected only when it is not already the
   * live one, with the step's `patch` written over its profile so the planet
   * is built once, as the step describes it — so a tour can state it on every
   * step that depends on it, and "back" lands on the right planet.
   */
  planet?: PlanetId;
  /**
   * Drag the convective-vigour slider for the reader over `ms`, rather than
   * snapping `logRa` to its destination: the point of the step is that the
   * picture responds continuously to a control, and a value that teleports
   * shows the ends without the middle.
   *
   * `from`, when given, is where the slider is put before it travels, so a
   * step about crossing a value always crosses it, whatever the reader left
   * the slider at on the step before.
   */
  ramp?: { from?: number; to: number; ms: number };
  /** Fly the 2-D camera. Ignored in the 3-D view, which has its own camera. */
  focus?: TourFocus;
  /** Draw a guide along the annulus' outer surface for this step. */
  surfaceGuide?: boolean;
  /** Start this step from the standard small thermal perturbation. */
  reseed?: boolean;
  /** Advance on its own once the effect has had time to develop. "Next" still skips ahead. */
  dwell?: TourDwell;
  /**
   * Offer a "replay" button that enters the step again from scratch — the
   * reseed, the ramp from its `from`, the dwell — for a step whose point is
   * something seen happening once, which a reader who looked away has missed.
   */
  replay?: boolean;
  /**
   * Fields of `State` that "finish" on this step puts back to what the reader
   * had before the tour, where the tour moved them. For settings a tour
   * needs while it runs but that the reader should not be left with, such as
   * a costly grid; unlike "restore the run I had", which puts back everything,
   * the run itself carries on.
   */
  restoreOnFinish?: readonly (keyof State)[];
  /**
   * The one line telling the reader what to actually look for, kept out of
   * `body` because the card styles it apart: the explanation is why the
   * control exists, this is what the canvas is about to do.
   */
  watch?: string;
}

/**
 * The outer boundary of the annulus, in the world coordinates `focus` is
 * written in — `RADIUS_INNER + 1`, since the shell is built with unit
 * thickness (see `presets.ts`). Named here rather than left as a literal in
 * the one step that flies there, because that step's whole point is *which*
 * boundary it is looking at.
 */
const OUTER_RADIUS = 2.208318891;

/**
 * Everything an onset test depends on, pinned by the tour's first step rather
 * than inherited: a reader arriving from a no-slip benchmark, from the box or
 * from a finer grid would otherwise be told a
 * threshold their run does not have. The planet is the one thing left alone —
 * switching profile is an animated rebuild of its own — so the thresholds
 * quoted are Earth's proportions, the ones the app opens on.
 *
 * The time-step cap is raised tenfold for the tour. Near onset the flow is
 * nearly still, so the Courant rule never binds and the cap alone sets how
 * fast simulated time passes; at the resolution's own cap a disturbance at
 * Ra = 100 takes most of a minute to fade. Growth rates here are at most a
 * few tens, so σ·dt stays below a few per cent — and the closing step puts
 * the resolution's own cap back.
 *
 * Stated on every experiment step, not just the first, together with the
 * domain that step runs in: "back" from the box must land in the annulus again,
 * and `tour.ts` drops the fields that already hold, so restating them costs
 * nothing — no rebuild unless the domain actually has to change.
 */
const ONSET_SETUP = {
  resolution: DEFAULT_PRESET, dtMax: 1e-3, courant: 1, speed: 2,
  particles: "off", viscosity: "constant", isothermal: false, paused: false,
} as const satisfies Partial<State>;

/** The annulus the app opens on, seeded with its lowest-threshold pattern (see below). */
const ONSET_ANNULUS = {
  ...ONSET_SETUP, geometry: "spherical annulus", radialWalls: "free-slip", wavenumber: 4,
} as const satisfies Partial<State>;

/**
 * One wavelength of the free-slip critical roll pair, k = π/√2 — the box
 * `tests/temperature.test.ts` checks the onset of convection in, seeded the
 * same way, so its threshold is 27π⁴/4 by construction.
 */
const ONSET_BOX = {
  ...ONSET_SETUP, geometry: "Cartesian box", boxLength: 2 * Math.SQRT2, walls: "periodic",
  wavenumber: 1,
} as const satisfies Partial<State>;

/**
 * Measured thresholds, rounded for prose: σ from the growth of v_rms either
 * side of onset, interpolated to zero (σ is linear in Ra for a fixed
 * pattern), on a 24 × 64 grid at the tour's own time-step cap. Seed mode
 * 4 is the lowest of the annulus' patterns: mode 2 ≈ 1,053, 3 ≈ 716,
 * four ≈ 673, five ≈ 751, six ≈ 916. The no-slip box, at the free-slip
 * box's width, ≈ 1,979: above the flat layer's 1,708 because that minimum
 * is for rolls about 2.0 deep-units wide, not this box's 2.83.
 * `tests/onset-tour.test.ts` re-checks the signs the prose relies on.
 */
/**
 * The card slider's span for every onset step, Ra = 100 to ≈ 3,160: from the
 * conductive step's Ra to past the no-slip box's threshold, and narrow
 * enough that a pixel of drag on the card is a percent or two of Ra.
 */
const ONSET_RA_RANGE = { min: 2, max: 3.5 } as const;

const ANNULUS_RA_C = "670";
/**
 * Every seed pattern's threshold in the annulus, and what it does at the
 * pattern step's Ra = 800 — the table on that step's card.
 *
 * Modes 1–6 are measured as above (mode 1 ≈ 3,164). Modes 7 and up are the
 * flat free-slip layer's (k² + π²)³/k² at each mode's mid-depth wavelength,
 * scaled by the 2–3% that formula falls short of the measured modes 1–6 by:
 * mode 10 ≈ 2,600, 11 ≈ 3,400 and 12 ≈ 4,400, the last two past the card
 * slider's end. Above threshold, a seed far from mode 4 grows but is
 * overtaken: measured on a 16 × 64 grid, mode 1 at Ra = 3,160 ends as five
 * cells, mode 2 at 1,500 as four, and mode 10 at 3,000 as six. The signs at
 * 800 are measured for every mode, on a 24 × 64 grid.
 */
const ANNULUS_MODES: TourTable = {
  head: ["Mode", "Needs Ra ≈", "At Ra = 800"],
  rows: [
    ["1", "3,160", "fades"],
    ["2", "1,050", "fades"],
    ["3", "720", "grows slowly"],
    ["4", ANNULUS_RA_C, "grows fastest"],
    ["5", "750", "grows slowly"],
    ["6", "920", "fades"],
    ["7", "1,200", "fades"],
    ["8", "1,500", "fades"],
    ["9", "2,000", "fades"],
    ["10+", "2,600+", "fades"],
  ],
};
const BOX_NO_SLIP_RA_C = "2,000";

/**
 * Surface gravity, m/s² (NASA planetary fact sheets). `planets.ts` has no use
 * for it — its profiles state Ra directly — so it lives here, with the one
 * tour that scales Ra from planet to planet.
 */
const SURFACE_GRAVITY: Record<PlanetId, number> = { mars: 3.71, venus: 8.87, earth: 9.81 };

/**
 * Ra ∝ g·d³, relative to Mars, with every other property in it held equal:
 * mantle thickness from the profiles' own radii. Venus ≈ 11.4, Earth ≈ 11.8.
 *
 * Deliberately not the profiles' own Ra, which order the planets the other
 * way (Mars 10⁷, Venus 10^6.5): each is a resolved teaching value chosen
 * against its own sources (see their caveats in `planets.ts`), not a
 * comparison between planets. The tour's claim *is* that comparison, so it
 * sets Ra itself.
 */
const raRelativeToMars = (id: PlanetId): number =>
  (SURFACE_GRAVITY[id] * radiiFor(PLANETS[id]).depthKm ** 3)
  / (SURFACE_GRAVITY.mars * radiiFor(MARS).depthKm ** 3);

/** Mars at 10⁵, the others scaled from it. */
const PLANET_LOG_RA: Record<PlanetId, number> = {
  mars: 5,
  venus: 5 + Math.log10(raRelativeToMars("venus")),
  earth: 5 + Math.log10(raRelativeToMars("earth")),
};

const formatRatio = (id: PlanetId): string => raRelativeToMars(id).toFixed(0);
const formatStress = (v: number): string => v.toLocaleString("en-US");

/** The comparison on the three-planet tour's "why Venus" card. */
const PLANET_SCALES: TourTable = {
  head: ["Planet", "Mantle depth d (km)", "g (m/s²)", "Ra relative to Mars"],
  rows: (["mars", "venus", "earth"] as const).map((id) => [
    PLANETS[id].label,
    radiiFor(PLANETS[id]).depthKm.toLocaleString("en-US"),
    SURFACE_GRAVITY[id].toFixed(2),
    id === "mars" ? "1" : `≈ ${raRelativeToMars(id).toFixed(1)}`,
  ]),
};

/**
 * The numerics every three-planet step runs at. The app's own default grid,
 * not a finer one: the finest is about four times slower a step, which makes
 * each lid take too long to form for the tour to be watchable, and the
 * default reproduces the same three regimes (see `STRONG_LID`). Courant 2 is
 * past the amber warning (accuracy, not stability — see the numerics help)
 * and roughly halves the wall-clock time a lid takes to form. The closing
 * step puts the reader's own values back on "finish".
 */
const PLANET_SETUP = {
  resolution: DEFAULT_PRESET, courant: 2, wavenumber: 5,
  isothermal: false, paused: false,
} as const satisfies Partial<State>;

/**
 * One rock for all three planets: Tosi et al. (2015)'s law, a 10⁵ thermal
 * contrast with no depth term and the benchmark's own η*, so that between
 * Venus and Earth only σ_Y differs.
 */
const TOSI = {
  viscosity: "Tosi", logContrast: 5, logDepthContrast: 0, sigmaB: 0, etaStar: 1e-3, picard: 1,
} as const satisfies Partial<State>;

/**
 * The two yield stresses, measured on the GPU solver on the default grid at
 * Courant 2, running the tour's own sequence — a fresh Mars seed, then Venus
 * and Earth each carrying the previous planet's field, as a planet change
 * does. The ratio of surface to interior v_rms, over the last 40% of 6,000
 * steps: Mars 0.02 and Venus 0.08 at 60,000 (stagnant lids; Venus' interior
 * about five times Mars'), Earth 0.94 at 10,000 (mobile). The finest grid
 * gives the same regimes (0.02, 0.12, 0.75 over 3,000 steps). At 20,000
 * Venus' lid is episodic (0.49), and Earth breaks anywhere from 3,000 to
 * 10,000. The Tosi law is bistable (see "Tosi 4" in `presets.ts`), and which
 * branch a run takes depends on the flow it starts from: on the finest grid,
 * Mars under this rock carrying the uniform-rock step's flow stays mobile at
 * 60,000, which is why its lid step reseeds — and a reseed in the Krylov
 * tier also clears ψ (`GpuSimulation.seedTemperatureDisturbance`).
 */
const STRONG_LID = 6e4;
const WEAK_LID = 1e4;

export const TOURS = {
  /**
   * The tour for someone who has just arrived: what the picture is, then the
   * two controls that change it most (how hard the layer is driven, and how
   * the rock responds), then the three overlays that explain what it is
   * doing. Ends in a live, configured run rather than back where it started —
   * see `tour.ts` on why that is the wanted ending and not a leak.
   */
  "First look": [
    {
      id: "welcome",
      title: "A planet's mantle, cut open",
      body: [
        "This globe is cut open to show one slice through the mantle: the "
        + "rocky shell between the planet's core and its surface. Colour is "
        + "temperature: hot at the core–mantle boundary on the inside, cold at "
        + "the surface on the outside.",
        "Rock here is solid, but over millions of years it flows. Hot rock "
        + "rises because it is less dense, cold rock sinks, and the whole "
        + "layer turns over. That is mantle convection, and it is what drives "
        + "plate tectonics.",
      ],
      target: "canvas",
      // Stated rather than inherited from the app's opening pose, so the
      // card's "this globe" is true however the reader left the view.
      view: "3d",
      watch: "Look for hot rock rising from the inner edge and cold rock "
        + "sinking from the outer one. This is running live, not a recording.",
    },
    {
      id: "flat-view",
      title: "The flat view the model runs in",
      body: [
        "The model doesn't solve the whole globe. It solves this one slice, "
        + "laid flat: the spherical annulus. Its inner ring is the core–mantle "
        + "boundary and its outer ring is the surface.",
        "You'll stay in this view for the rest of the tour. The \"3D view\" "
        + "button, lit here, takes you back to the globe at any time.",
      ],
      target: "view3d",
      view: "2d",
    },
    {
      id: "sluggish",
      title: "Just past the onset of convection",
      body: [
        "The \"try an example\" list loads a ready-made scene. This one, "
        + "\"Barely convecting\", has only just enough vigour to move: one or "
        + "two lazy cells, with most of the layer doing very little.",
        "Convection doesn't fade in gradually. Below a critical vigour nothing "
        + "moves, and heat crosses the layer by conduction alone. Above it, "
        + "the layer turns over on its own.",
      ],
      target: "preset",
      preset: "Barely convecting",
      // The example states a problem, not a playback state — so the two
      // fields that would otherwise leave this step watching a still frame
      // (the isothermal override forces Ra to 0; a paused solver takes no
      // steps for the dwell below to count) are corrected on top of it.
      patch: { isothermal: false, paused: false },
      dwell: { steps: 400 },
      watch: "One or two broad, slow cells, and long stretches where little changes.",
    },
    {
      id: "vigour",
      title: "Convective vigour: the Rayleigh number",
      body: [
        "How vigorously a layer convects is measured by one number, the "
        + "Rayleigh number, Ra. It is a ratio: buoyancy, which lifts hot rock "
        + "and sinks cold rock, divided by what resists it. Viscosity slows "
        + "the flow, and thermal diffusion lets a warm blob lose its heat "
        + "before it can rise. The larger Ra is, the more decisively buoyancy "
        + "wins.",
        "Ra has no units, so the same value means the same kind of flow "
        + "whatever the size of the layer. It grows with the cube of the "
        + "layer's thickness, which is how a mantle thousands of kilometres "
        + "deep reaches such large values.",
        "The \"convective vigour\" slider sets Ra, on a logarithmic scale: "
        + "each step along it multiplies Ra rather than adding to it. "
        + "\"Barely convecting\" ran at Ra = 2,000, just above the critical "
        + "value where convection first starts. Watch it rise to a million "
        + "(10⁶). Earth's mantle is estimated at 10⁷–10⁸, beyond the top of "
        + "this ramp.",
      ],
      target: "vigour",
      ramp: { to: 6, ms: 5000 },
      dwell: { steps: 300 },
      watch: "The cells break into many thin, fast upwellings and "
        + "downwellings, and the hot and cold layers along the inner and outer "
        + "rings grow thinner.",
    },
    {
      id: "rock",
      title: "How the rock behaves",
      body: [
        "So far every part of the layer has been equally stiff. Real rock is "
        + "not: it is far stiffer where it is cold. \"How the rock behaves\" "
        + "is now set to \"stiffer when cold\", which makes the coldest rock a "
        + "thousand times as viscous as the hottest. Vigour is back down to "
        + "Ra = 10⁵, so the change is easier to see.",
        "Stiffening the cold rock turns a field of interchangeable cells into "
        + "a few broad, long-lived upwellings. These are what geologists call "
        + "plumes.",
      ],
      target: "rock",
      patch: { viscosity: "μ(T, d)", logRa: 5, logContrast: 3, logDepthContrast: 0 },
      dwell: { steps: 500 },
      watch: "A stiff cold lid forming along the outer ring, and plume stems "
        + "that persist instead of drifting apart.",
    },
    {
      id: "boundary-layer",
      title: "The thermal boundary layer",
      body: [
        "We've zoomed in on the outer ring of the annulus. The thin cold band "
        + "along it is the thermal boundary layer, this model's version of a "
        + "planet's lithosphere.",
        "Nearly the entire temperature drop across the mantle happens inside "
        + "that band. The interior beneath it is close to uniform in "
        + "temperature, because convection stirs it faster than conduction "
        + "can build a gradient. The band thickens until it is too heavy, then "
        + "drips off and sinks.",
      ],
      target: "canvas",
      view: "2d",
      // Just inside the outer boundary, at the top of the annulus. The zoom
      // is limited by `clampPan` in main.ts, which keeps at least 1/zoom of
      // the domain on screen — at 5.5 this frames the cold lid together with
      // the upper half of the interior it is falling into, which is the
      // comparison the step is making.
      focus: { zoom: 5.5, x: 0, y: OUTER_RADIUS - 0.26, ms: 1500 },
      surfaceGuide: true,
      // Counted in solver steps like the other steps waiting on an effect: a
      // drip detaching is physics, and takes longer on a slower machine.
      dwell: { steps: 400 },
      watch: "A cold finger thickening, detaching, and sinking away from the surface.",
    },
    {
      id: "flow-lines",
      title: "Flow lines",
      body: [
        "\"Show flow lines\" is now switched on, and rock moves along these "
        + "lines. They are streamlines: contours of the stream function, which "
        + "is what the solver computes.",
        "Closed loops are convection cells. Where the lines crowd together the "
        + "flow is fast, and where they are far apart it is nearly still.",
      ],
      target: "flow",
      patch: { contours: 24 },
      focus: "reset",
      dwell: { steps: 300 },
      watch: "Loops merging and splitting as neighbouring plumes compete.",
    },
    {
      id: "diagnostics",
      title: "Reading the run",
      body: [
        "The lower panel plots the Nusselt number, Nu. Its two curves are the "
        + "heat flowing through the core–mantle boundary (inner) and through "
        + "the surface (outer), each divided by what conduction alone would "
        + "carry. Nu = 1 means no convection at all; Nu = 10 means convection "
        + "is carrying ten times the heat conduction alone would.",
        "The upper panel plots root-mean-square velocity: how fast the whole "
        + "layer is moving (v_rms), and how fast its surface is moving "
        + "(surface v_rms). Surface velocity is the model's closest stand-in "
        + "for the speed of tectonic plates, and the right-hand axis converts "
        + "it to cm/yr.",
        "For comparison, the Atlantic Ocean widens by about 2–4 cm/yr as the "
        + "plates on either side of it pull apart. North America and Europe "
        + "are separating at around 2 cm/yr, and South America and Africa at "
        + "around 3–4 cm/yr.",
      ],
      target: "traces",
      dwell: { ms: 7000 },
      watch: "A flat trace is a steady state result. A wobbling one is "
        + "time-dependent flow, not noise in the measurement.",
    },
    {
      id: "section-help",
      title: "Advanced controls, and where to learn more",
      body: [
        "Everything else is behind \"advanced controls\", which the tour has "
        + "just switched on: the domain and its boundary conditions, the "
        + "numerics, viscosity laws including yielding and power-law creep, "
        + "the initial condition, the view and the tracers, each in its own "
        + "section below.",
        "If you want to learn more about a section, click the ? beside its "
        + "title. It walks through that section's controls one at a time, "
        + "explaining what each does and the physics behind it, without "
        + "changing your run.",
      ],
      target: "domainHelp",
      highlight: "advanced",
      advanced: true,
      watch: "Every advanced section has its own ?, like the one lit here beside \"domain\".",
    },
    {
      id: "done",
      title: "Over to you",
      body: [
        "Scroll to zoom and drag to pan the canvas at any time, and press H "
        + "to hide the interface.",
        "The other guided tutorials pick up from here: \"convection onset\" "
        + "looks closely at the threshold you saw at the start, and \"three "
        + "planet tour\" follows Mars, Venus and Earth to show why only "
        + "Earth's surface moves as plates. Published benchmark "
        + "cases from Blankenbach, Tosi and van Keken are in the \"try an "
        + "example\" list, below the ready-made scenes.",
        "The run is left exactly where the tour finished, so you can carry "
        + "on from here, or use \"restore the run I had\" below to put it back.",
      ],
      target: "tutorials",
    },
  ],
  /**
   * A linear-stability experiment the reader runs by hand: seed, watch the
   * disturbance fade or grow, move the threshold's one control, repeat. First
   * in the annulus the app opens on, then in the flat box where the answer is
   * known exactly (27π⁴/4, which `tests/temperature.test.ts` pins the solver
   * to), then with the one change that moves it most — boundaries that grip.
   *
   * Every threshold quoted in the prose is a measured one, and
   * `tests/onset-tour.test.ts` re-measures each claim against the step's own
   * settings on a coarse grid, so a solver or preset change that moves one
   * fails there rather than in front of a reader.
   */
  "Convection onset": [
    {
      id: "conduction",
      title: "Conduction before convection",
      body: [
        "This experiment starts from rock of uniform viscosity and a "
        + "deliberately weak thermal drive, Ra = 100. The highlighted \"seed "
        + "disturbance\" button resets the temperature to the conductive "
        + "profile plus a small harmonic perturbation, here mode 4, chosen with "
        + "\"seed pattern\" just above it: four warm and "
        + "four cool lobes spaced evenly around the annulus, strongest mid-mantle "
        + "and vanishing at both boundaries, at 5% of the temperature "
        + "difference across the layer.",
        "The warm lobes are lighter than the rock around them and the cool "
        + "lobes heavier, so together they set the rock moving. The test is "
        + "whether that motion carries heat in a way that strengthens the "
        + "temperature disturbance, or whether it is smoothed away first.",
        "Here the temperature disturbance fades: at this Ra, convection cannot "
        + "start. Viscosity and thermal diffusion remove it faster than "
        + "buoyancy can feed it, and heat crosses the layer by conduction "
        + "alone: hot at the inner ring, cold at the outer, varying smoothly "
        + "in between.",
        "Both rings of the annulus are free-slip: rock cannot cross them, but "
        + "slides along them without friction. That is a fair match for a "
        + "planet. At the bottom, the liquid iron of the outer core is far too "
        + "runny to grip the mantle; at the top, only ocean or air lies above.",
      ],
      target: "seed",
      highlight: "seedPattern",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      view: "2d",
      patch: { ...ONSET_ANNULUS, logRa: 2 },
      reseed: true,
      dwell: { steps: 400 },
      watch: "Watch the warm and cool lobes smooth back into the conductive "
        + "profile, and the root mean square velocity (v_rms), bottom left, "
        + "fall toward zero: at Ra = 100 the drive is too weak to convect. "
        + "Ra, on this card's slider, sets how hard the layer is driven, so "
        + "raising it is how to find where convection begins. After each "
        + "change, press seed disturbance: v_rms falling means conduction "
        + "still wins; v_rms rising means convection has started.",
    },
    {
      id: "reading-the-test",
      title: "Reading an onset test",
      body: [
        "These two plots are how to read the test, and \"seed disturbance\" "
        + "clears both so each test starts afresh.",
        "v_rms, above, is how fast the layer is moving. Right after reseeding "
        + "the changes in v_rms can be small. However, what matters is in "
        + "which direction it changes. Falling means the disturbance is dying "
        + "and conduction wins. Rising means the disturbance is feeding "
        + "itself, i.e., convection.",
        "Nu, below, stays at 1 until the flow is strong enough to carry heat "
        + "in earnest, so it responds later than v_rms. A disturbance that has "
        + "grown into a full circulation lifts Nu clearly above 1.",
      ],
      target: "traces",
    },
    {
      id: "approach-threshold",
      title: "Find the threshold yourself",
      body: [
        "Raise Ra with the slider on this card, or type a value into the "
        + "box above it. \"convective vigour\" in the options panel and the Ra "
        + "number in this card are synchronised. After each change, press \"seed disturbance\" and watch v_rms.",
        "Somewhere between 500 and 1,000 v_rms changes from falling to "
        + "rising. Near the threshold, the temperature disturbance changes "
        + "very slowly. Just below it, v_rms falls only gradually; just above "
        + "it, v_rms rises only gradually. The closer Ra is to the threshold, "
        + "the slower the change, so a disturbance that neither clearly grows "
        + "nor clearly fades is a sign you are close.",
      ],
      target: "vigour",
      highlight: "seed",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      patch: ONSET_ANNULUS,
      watch: "Raise Ra a step at a time, pressing seed disturbance after "
        + "each. Below the threshold v_rms falls after the seed; above it, "
        + "v_rms climbs. The lowest Ra at which v_rms climbs is the threshold "
        + "for this seed pattern.",
    },
    {
      id: "onset",
      title: "Crossing into convection",
      body: [
        "The slider now travels from about 320 to about 1,600, starting "
        + "from a fresh disturbance. In this annulus, the mode-4 disturbance "
        + `starts to grow at Ra ≈ ${ANNULUS_RA_C}. Above that, buoyancy `
        + "amplifies a disturbance faster than viscosity and diffusion can "
        + "remove it, and it grows until it is a steady circulation carrying "
        + "heat. Below it, conduction is the stable state. It is a threshold, "
        + "not a gradual increase.",
      ],
      target: "vigour",
      highlight: "seed",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      patch: ONSET_ANNULUS,
      reseed: true,
      ramp: { from: 2.5, to: 3.2, ms: 3000 },
      dwell: { steps: 800 },
      replay: true,
      watch: "Ra now rises on its own. While it is below about "
        + `${ANNULUS_RA_C}, v_rms drifts down; once Ra passes that threshold, `
        + "v_rms turns and climbs. Nu follows later, rising above 1 as the "
        + "disturbance grows into full convection cells.",
    },
    {
      id: "patterns",
      title: "Each pattern has its own threshold",
      body: [
        "\"seed pattern\" chooses the disturbance's harmonic mode, i.e., how "
        + "many times it repeats around the annulus.",
        "Each mode has its own threshold, the Ra it needs before it grows. "
        + "Ra is now 800.",
        ANNULUS_MODES,
        "Mode 4 needs the least. Narrow cells lose their heat sideways to "
        + "their neighbours before it can drive them. Wide cells must push "
        + "rock a long way sideways for every rise and fall, against "
        + "viscosity. The cell width that convects most easily lies in "
        + "between, and in this annulus that is mode 4.",
        "Modes far from 4 do not keep their shape. Above its threshold such a "
        + "pattern grows, but the modes near 4 grow much faster. The flow "
        + "itself stirs a little of them into the layer, and they soon take "
        + "over. Seeded at Ra = 1,500, mode 2 ends as four cells. Seeded at "
        + "Ra = 3,000, mode 10 ends as six. Cells that are too wide split and "
        + "cells that are too narrow merge, until the pattern is close to the "
        + "width that convects most easily.",
        "Pick a mode, press \"seed disturbance\", and watch v_rms.",
      ],
      target: "seedPattern",
      highlight: "seed",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      patch: { ...ONSET_ANNULUS, logRa: Math.log10(800) },
      reseed: true,
      watch: "v_rms doing what the table's last column says. Then raise Ra "
        + "past another mode's threshold and seed it: it grows too, but a "
        + "mode far from 4 soon reorganises into cells near mode 4's width.",
    },
    {
      id: "box-free-slip",
      title: "The textbook threshold",
      body: [
        "For a flat layer of uniform rock with free-slip top and bottom, like "
        + "the rings of the annulus so far, the threshold is known exactly. "
        + "In 1916 Rayleigh found it to be Ra = 27π⁴/4 ≈ 657.5, reached first "
        + "by a wavelength of 2√2 ≈ 2.8 times the layer's depth.",
        "The run is now that layer, in a box exactly one such wavelength "
        + "wide, seeded with mode 1, at Ra ≈ 500. Its rolls are nearly the same width "
        + "as the cells of mode 4 in the annulus, which is why their thresholds "
        + "nearly agree. This app's own test suite checks the solver against 657.5.",
      ],
      target: "vigour",
      highlight: "seed",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      patch: { ...ONSET_BOX, radialWalls: "free-slip", logRa: 2.7 },
      reseed: true,
      watch: "The rolls fading. Bracket 657.5 with the card. At 600 they "
        + "fade and at 720 they grow, both slowly.",
    },
    {
      id: "box-no-slip",
      title: "Boundaries that grip",
      body: [
        "Now the top and bottom are no-slip. Rock touching them cannot move "
        + "at all, neither through them nor along them. Everything else is "
        + "unchanged, and Ra is back to about 800, where the rolls grew a "
        + "moment ago.",
        "Few planetary boundaries grip like this, but it is a useful stand-in "
        + "for a thick lid that does not move. That is roughly the situation "
        + "beneath the stagnant lids of Mars and Mercury, whose cold outer "
        + "shells do not break into moving plates. It is also how convection "
        + "is studied in the laboratory, with a fluid layer heated between "
        + "rigid plates.",
        "Buoyancy drives the flow, and two things work against it. Viscosity "
        + "resists the rock being sheared as it moves, and thermal diffusion "
        + "evens out the temperature differences that make the rock buoyant. "
        + "With free slip, rock slides along the boundaries, so viscosity "
        + "only resists the flow turning inside the layer. With no slip, the "
        + "flow must also come to a stop at each boundary, so the rock beside "
        + "it is sheared hard. That extra viscous drag along both boundaries "
        + "is more resistance for buoyancy to overcome, so a larger Ra is "
        + "needed before the rolls can grow.",
        "For the wavelength that suits these boundaries best, about 2.0 times "
        + "the depth, the threshold rises to about 1,708, the value laboratory "
        + "experiments measure. In this box, sized for free slip, it is about "
        + `${BOX_NO_SLIP_RA_C}.`,
      ],
      target: "vigour",
      highlight: "seed",
      companion: "traces",
      raControl: ONSET_RA_RANGE,
      patch: { ...ONSET_BOX, radialWalls: "no-slip", logRa: 2.9 },
      reseed: true,
      dwell: { steps: 400 },
      watch: `The same rolls fading at the same Ra. Raise Ra past about `
        + `${BOX_NO_SLIP_RA_C}, seed, and they grow again.`,
    },
    {
      id: "why-this-value",
      title: "What this tutorial showed",
      body: [
        "In summary:",
        {
          items: [
            "Convection starts only above a threshold Ra. Below it, a "
            + "disturbance fades and heat crosses the layer by conduction. "
            + "Above it, the disturbance grows into convection cells.",
            "At the threshold, buoyancy just balances viscosity and thermal "
            + "diffusion. Anything that shifts that balance moves the "
            + "threshold.",
            "To test for it, seed a disturbance and watch v_rms. Falling "
            + "means conduction wins and rising means convection. Nu rises "
            + "above 1 once the cells carry heat.",
            "Each seed pattern has its own threshold. In this annulus mode 4 "
            + `needs the least, Ra ≈ ${ANNULUS_RA_C}, and patterns far from `
            + "it reorganise towards its cell width.",
            "A flat free-slip layer's threshold is exactly Ra = 27π⁴/4 ≈ "
            + "657.5. "
            + "The annulus lands close to it but is a different problem, "
            + "because its curvature concentrates the heat entering from "
            + "below.",
            "No-slip boundaries add viscous drag and raise the threshold, to "
            + "Ra ≈ 1,708 for the best-suited cell width.",
          ],
        },
        "The run is left in the box. Pick the annulus again under advanced "
        + "controls → domain, or use \"restore the run I had\" below.",
      ],
      target: "tutorials",
      // Back to the resolution's own cap: the raised one was for watching
      // slow, near-conductive runs, and the run carries on from here.
      patch: { dtMax: PRESETS[DEFAULT_PRESET].dtMax },
    },
  ],
  /**
   * Why only Earth's surface moves: Mars, Venus and Earth under one rock,
   * changing one thing at a time. Mars first with uniform rock (the control),
   * then with rock that stiffens with cold and yields under stress, which
   * gives it a stagnant lid; Venus at the Ra its size gives it, same rock,
   * still a lid; Earth at nearly Venus' Ra with a weaker lid, which breaks.
   *
   * Every planet step states its whole model (`PLANET_SETUP`, `TOSI`), so
   * stepping back lands on the model the card describes, and `tour.ts` drops
   * the fields that already hold, so restating them costs nothing.
   */
  "Three planet tour": [
    {
      id: "mars-uniform",
      title: "Mars, if its rock were uniform",
      body: [
        "This tour visits three rocky planets, Mars, Venus and Earth, to ask "
        + "why only Earth's surface is broken into moving plates.",
        "It starts with Mars as a control experiment. The Rayleigh number is "
        + "10⁵, and the rock is equally stiff everywhere, hot or cold. With "
        + "nothing stiffer at the surface, the whole mantle turns over, and "
        + "the surface moves with it.",
      ],
      target: "planet",
      companion: "traces",
      view: "3d",
      planet: "mars",
      patch: { ...PLANET_SETUP, viscosity: "constant", logRa: PLANET_LOG_RA.mars },
      reseed: true,
      dwell: { steps: 600 },
      watch: "The globe turns to Mars. Then, in the upper plot at the bottom "
        + "left, surface v_rms keeps pace with v_rms: the surface moves as fast "
        + "as the mantle beneath it.",
    },
    {
      id: "mars-lid",
      title: "Mars: a stagnant lid",
      body: [
        "Now the rock follows the Tosi law, named after the 2015 community "
        + "benchmark it comes from. Cold rock is 100,000 times stiffer than hot "
        + "rock, and any rock can yield: where the stress in it exceeds a yield "
        + "stress, σ_Y, it breaks and flows easily.",
        `Here σ_Y is high, ${formatStress(STRONG_LID)}, so the cold rock at the `
        + "surface never breaks. It forms a stagnant lid: a rigid shell the "
        + "mantle convects beneath but cannot move. That is thought to be Mars "
        + "today, and why it has no plate tectonics.",
        "The run restarts from a small disturbance so the lid forms in front "
        + "of you, and the view has zoomed in on the surface at the top of the "
        + "annulus.",
      ],
      target: "rock",
      companion: "traces",
      view: "2d",
      planet: "mars",
      patch: { ...PLANET_SETUP, ...TOSI, sigmaY: STRONG_LID, logRa: PLANET_LOG_RA.mars },
      reseed: true,
      focus: { zoom: 5.5, x: 0, y: radiiFor(MARS).ro - 0.26, ms: 1500 },
      surfaceGuide: true,
      dwell: { steps: 800 },
      watch: "A thick cold band along the surface that never moves, with plumes "
        + "rising beneath it. Surface v_rms stays near zero while v_rms climbs.",
    },
    {
      id: "why-venus",
      title: "Why Venus convects harder",
      body: [
        "Venus is nearly Earth's twin in size, and much larger than Mars. How "
        + "hard a mantle convects is set by its Rayleigh number, and Ra grows "
        + "with gravity and with the cube of the mantle's thickness: Ra ∝ g·d³.",
        PLANET_SCALES,
        "Every other property is held the same, so Venus runs at about "
        + `${formatRatio("venus")} times Mars' Ra, and Earth at about `
        + `${formatRatio("earth")} times.`,
        "Venus keeps Mars' rock and yield stress exactly. Before moving on, "
        + "predict which of these happens:",
        {
          items: [
            "The lid breaks, and the surface starts to move.",
            "The lid holds, and only the mantle beneath it speeds up.",
          ],
        },
      ],
      target: "vigour",
      companion: "traces",
      view: "2d",
      planet: "mars",
      patch: { ...PLANET_SETUP, ...TOSI, sigmaY: STRONG_LID, logRa: PLANET_LOG_RA.mars },
      focus: "reset",
    },
    {
      id: "venus",
      title: "Venus: faster, but still a lid",
      body: [
        `The tour has moved to Venus, with Ra about ${formatRatio("venus")} `
        + "times Mars' and the rock unchanged: the same stiffening with cold, "
        + "the same yield stress.",
        "The mantle convects much harder, with more plumes rising faster. But "
        + "the lid holds. Stress in the lid grows with the vigour of the flow, "
        + "but not enough to reach this yield stress. Venus is thought to have "
        + "a stagnant lid today, like Mars, though its surface may have been "
        + "renewed in episodes of overturn.",
      ],
      target: "planet",
      companion: "traces",
      view: "3d",
      planet: "venus",
      patch: { ...PLANET_SETUP, ...TOSI, sigmaY: STRONG_LID, logRa: PLANET_LOG_RA.venus },
      dwell: { steps: 800 },
      watch: "The globe turns to Venus. v_rms settles several times higher than "
        + "on Mars, while surface v_rms stays far below it.",
    },
    {
      id: "earth",
      title: "Earth: a lid that breaks",
      body: [
        "Earth is about the size of Venus, so its Ra is almost the same, about "
        + `${formatRatio("earth")} times Mars'. The one change that matters is `
        + `the yield stress: ${formatStress(WEAK_LID)} instead of `
        + `${formatStress(STRONG_LID)}.`,
        "Now the lid breaks. Where the stress in it exceeds σ_Y it fails, and "
        + "slabs of cold lid sink back into the mantle. The surface moves with "
        + "the flow beneath it: this model's version of plate tectonics, a "
        + "mobile lid.",
        "Why Earth's lithosphere is weaker than Venus' is still debated. One "
        + "leading idea is water, which weakens rock: Venus' surface is "
        + "extremely dry.",
      ],
      target: "planet",
      companion: "traces",
      view: "3d",
      planet: "earth",
      patch: { ...PLANET_SETUP, ...TOSI, sigmaY: WEAK_LID, logRa: PLANET_LOG_RA.earth },
      dwell: { steps: 800 },
      watch: "The globe turns to Earth. Surface v_rms rises close to v_rms: the "
        + "surface now moves with the mantle.",
    },
    {
      id: "earth-plates",
      title: "Plates, up close",
      body: [
        "Zoomed in on the surface again, the cold band is no longer a fixed "
        + "shell. It thins where the surface pulls apart, and peels away and "
        + "sinks where it converges.",
        "The lit slider is the yield stress σ_Y, in the advanced viscosity "
        + "controls. Dragging it up makes the lid harder to break; dragging it "
        + "down breaks it more easily.",
        "The upper plot's right-hand axis converts surface velocity to cm/yr. "
        + "Earth's surface here moves at up to a few centimetres a year, the "
        + "same order as real plates. Venus' lid crept several times more "
        + "slowly, and Mars' hardly moved at all.",
      ],
      target: "yieldStress",
      companion: "traces",
      advanced: true,
      view: "2d",
      planet: "earth",
      patch: { ...PLANET_SETUP, ...TOSI, sigmaY: WEAK_LID, logRa: PLANET_LOG_RA.earth },
      focus: { zoom: 5.5, x: 0, y: radiiFor(EARTH).ro - 0.26, ms: 1500 },
      surfaceGuide: true,
      dwell: { steps: 600 },
      watch: "Cold rock peeling away from the surface and sinking, and surface "
        + "v_rms close to v_rms.",
    },
    {
      id: "summary",
      title: "What the three planets showed",
      body: [
        "In summary:",
        {
          items: [
            "Where cold rock is much stiffer than hot rock, here 100,000 times, "
            + "a lid forms at the surface.",
            "Mars convects gently beneath a stagnant lid.",
            `Venus convects harder, at about ${formatRatio("venus")} times `
            + "Mars' Ra, because its mantle is thicker and its gravity stronger. "
            + "Its lid still holds.",
            "Earth has almost Venus' Ra. Only its lower yield stress lets the "
            + "lid break into moving plates.",
            "A planet's size sets how hard it convects; the strength of its lid "
            + "decides whether its surface moves.",
          ],
        },
        "Pressing finish leaves Earth running, and puts back the grid, "
        + "Courant number and time-step cap you had before the tour. A change "
        + "of grid restarts the run from a small disturbance. \"Restore the run "
        + "I had\" below puts back everything instead.",
      ],
      target: "tutorials",
      focus: "reset",
      restoreOnFinish: ["resolution", "courant", "dtMax"],
    },
  ],
} as const satisfies Record<string, readonly TourStep[]>;

export type TourName = keyof typeof TOURS;

/** The tour `start()` falls back to when no name is given; the pane's "guided tutorials" folder names one explicitly. */
export const DEFAULT_TOUR: TourName = "First look";

/**
 * A step of a section's help: the subset of `TourStep` that only *says*
 * something. Every field that would drive the model — patch, preset, planet,
 * ramps, reseed, camera, view — is absent from the type rather than merely
 * unused, so a help step that tried to change the run the reader opened it
 * over fails to compile. The help is read over a run, never staged on one.
 *
 * `target` is never `null`: each step is about one control in the section,
 * and that control is the thing lit.
 */
export type SectionHelpStep =
  Pick<TourStep, "id" | "title" | "body" | "highlight" | "watch">
  & { target: TourTargetName };

/**
 * The "?" beside each advanced folder's title (`controls.ts`): a short walk
 * through that folder's controls, one card per control, run by the same
 * overlay as the guided tours. Keyed by the folder's own title.
 *
 * Written for both geometries and every law at once, because the help does
 * not change the run to suit itself: a control that does not apply to the
 * current setup is either greyed out (box width, left / right on the
 * annulus), which its card says, or hidden (the viscosity and tracer knobs a
 * law or mode does not use), which `tour.ts` handles by skipping its card.
 */
export const SECTION_HELP = {
  "domain": [
    {
      id: "geometry",
      title: "Computational geometry",
      body: [
        "The computational geometry is the shape of the region the equations "
        + "are solved in. The model never simulates a whole planet: it solves "
        + "for flow and heat inside this one domain, and its edges are where the "
        + "physics has to be closed off.",
        "The spherical annulus is a slice through a spherical "
        + "shell. It keeps the curvature of a real mantle: the core–mantle "
        + "boundary is much shorter than the surface, so heat entering from "
        + "below is concentrated before it spreads out. The Cartesian box is a "
        + "flat slab, simpler and cheaper, and the setting most published "
        + "benchmark cases are stated in.",
      ],
      target: "geometry",
      watch: "Changing the geometry rebuilds the solver, which takes a second or two.",
    },
    {
      id: "box-width",
      title: "Box width",
      body: [
        "In the Cartesian box, this is the width of the domain in units of its "
        + "depth. The layer is always one unit deep, so a width of 4 is a box "
        + "four times wider than it is tall: its aspect ratio.",
        "Convection cells tend to be roughly as wide as the layer is deep, so "
        + "the width decides how many fit side by side. Benchmark cases fix it "
        + "to match the paper they reproduce.",
      ],
      target: "boxWidth",
      watch: "Box only: greyed out while the annulus is selected.",
    },
    {
      id: "walls",
      title: "Left / right: boundary conditions",
      body: [
        "A boundary condition states what the flow and temperature must do at "
        + "the edge of the domain. Without one the equations have no single "
        + "answer; with a different one, the same interior physics produces a "
        + "different flow.",
        "Periodic joins the left and right edges, so rock leaving one side "
        + "re-enters from the other, as though the box were one repeat of an "
        + "endlessly wide layer. Free-slip walls let nothing cross them and no "
        + "heat escape through them, but rock slides along them without "
        + "friction.",
      ],
      target: "walls",
      watch: "Box only: the annulus closes on itself, so it has no left or right edge.",
    },
    {
      id: "radial-walls",
      title: "Inner / outer: boundary conditions",
      body: [
        "These close the domain at the core–mantle boundary and at the "
        + "surface, labelled inner / outer on the annulus and bottom / top in "
        + "the box. Both are held at a fixed temperature, hot below and cold "
        + "above; this list chooses how the rock may move against them.",
        "Free-slip: no rock crosses the boundary, but it slides along it "
        + "without friction. No-slip: the rock sticks to the boundary and "
        + "cannot move along it either, like a mantle under a rigid lid. The "
        + "extra drag slows the flow beside that boundary and makes convection "
        + "harder to start.",
      ],
      target: "radialWalls",
    },
    {
      id: "resolution",
      title: "Resolution: discretising the domain",
      body: [
        "A computer cannot solve the equations at every point of a continuous "
        + "domain. Instead it discretises it: divides it into a grid of cells "
        + "and solves for values on that grid. The numbers are the size of the "
        + "grid the flow (ψ) is solved on: points across the depth × points "
        + "around the layer.",
        "A finer grid resolves thinner boundary layers and narrower plumes, "
        + "which matters most at high vigour, but the numerical expense grows "
        + "quickly. Doubling the grid in both directions gives four times the "
        + "cells, and smaller cells force a smaller time step too, so the same "
        + "stretch of simulated time costs roughly eight times as much.",
      ],
      target: "resolution",
      watch: "Changing resolution rebuilds the solver and resets the time-step cap.",
    },
  ],
  "numerics": [
    {
      id: "isothermal",
      title: "Isothermal (Ra = 0)",
      body: [
        "Switches thermal buoyancy off: temperature is still carried and "
        + "diffused, but it no longer makes rock rise or sink. The vigour "
        + "slider is greyed out while this is on, since its value is not being "
        + "used.",
        "This is for purely compositional problems, where density differences "
        + "come from chemistry rather than heat, such as the van Keken "
        + "Rayleigh–Taylor benchmark. Those need chemical tracers switched on, "
        + "or nothing drives the flow at all.",
      ],
      target: "isothermal",
    },
    {
      id: "courant",
      title: "Courant number",
      body: [
        "The solver advances in discrete time steps. The Courant number is "
        + "how far the fastest-moving rock may travel in one of them, measured "
        + "in grid cells: 1 means one cell per step. The step is resized "
        + "continually to hold that as the flow speeds up and slows down.",
        "Larger values take bigger steps, so simulated time passes faster, "
        + "but less accurately. The scheme used here stays stable well past 1, "
        + "so the trade is accuracy rather than a crash; the number turns amber "
        + "above 1 and red above 3 as a warning.",
      ],
      target: "courant",
    },
    {
      id: "dt-max",
      title: "Time-step cap",
      body: [
        "An upper limit on the time step, whatever the Courant number would "
        + "allow. When the flow is nearly still the Courant rule alone would "
        + "permit enormous steps, and the cap keeps heat diffusion and the "
        + "first growth of new instabilities accurately resolved.",
        "Each resolution sets its own cap when it is chosen, smaller for finer "
        + "grids. Raise it to cross slow, conductive stretches faster; lower it "
        + "if a run looks too coarse in time.",
      ],
      target: "dtMax",
    },
    {
      id: "dt-initial",
      title: "Initial time step",
      body: [
        "The step size a new solver starts with, before it has measured any "
        + "flow speed to size the step from. After the first few steps the "
        + "Courant rule and the cap take over.",
        "It is only read when the solver is next rebuilt, for example after "
        + "changing the resolution or geometry, so changing it has no effect "
        + "on the run already going.",
      ],
      target: "dtInitial",
    },
  ],
  "viscosity": [
    {
      id: "law",
      title: "Viscosity law",
      body: [
        "Viscosity is the rock's resistance to flow. Mantle rock is solid, "
        + "but over millions of years it creeps like an extremely stiff fluid, "
        + "and how its viscosity varies shapes the convection more than almost "
        + "anything else.",
        "Each law is a formula for that viscosity, written out below the list "
        + "with the value of each symbol and the slider that sets it. Constant "
        + "is the simplest case. The others add stiffening with cold and "
        + "depth, weakening where rock deforms quickly, yielding at high "
        + "stress, or dependence on composition; several are the exact laws "
        + "of published benchmark papers.",
      ],
      target: "law",
      highlight: "equation",
      watch: "Only the sliders the chosen law uses are shown beneath it.",
    },
    {
      id: "contrast",
      title: "Temperature contrast",
      body: [
        "How much stiffer the coldest rock is than the hottest, as a power of "
        + "ten: 3 means cold rock is a thousand times more viscous. It sets γ "
        + "(b in the Blankenbach law) in the equation above.",
        "A large contrast makes the cold top of the mantle too stiff to take "
        + "part in the flow. It forms a stagnant lid, the regime of Mars and "
        + "Venus today, with convection confined beneath it.",
      ],
      target: "contrast",
    },
    {
      id: "depth-contrast",
      title: "Depth contrast",
      body: [
        "How much stiffer rock at the bottom of the mantle is than at the top "
        + "at the same temperature, again as a power of ten. It sets c in the "
        + "equation. Pressure rises with depth and squeezes the rock, so real "
        + "mantle viscosity increases downward.",
        "Stiffening the deep mantle slows the flow there and tends to widen "
        + "the convection cells.",
      ],
      target: "depthContrast",
    },
    {
      id: "power-law-n",
      title: "Power-law exponent n",
      body: [
        "At n = 1 the rock is Newtonian: its viscosity does not depend on how "
        + "fast it is being deformed. Above 1 it is shear-thinning, weaker the "
        + "faster it deforms, like toothpaste. Values near 3 represent "
        + "dislocation creep, thought to dominate the upper mantle.",
        "Shear-thinning concentrates deformation into narrow, fast-moving "
        + "zones, one ingredient of plate-like behaviour.",
      ],
      target: "powerLawN",
    },
    {
      id: "cg-iterations",
      title: "CG iterations",
      body: [
        "With variable viscosity the flow equations can no longer be solved in "
        + "one direct step, so they are solved iteratively: a guess is refined "
        + "repeatedly by the conjugate gradient (CG) method. This sets how many "
        + "refinements each time step gets.",
        "More iterations give a more accurate flow at a direct cost in speed, "
        + "since each one is a full pass over the grid. Larger viscosity "
        + "contrasts converge more slowly and need more.",
      ],
      target: "cgIterations",
    },
    {
      id: "picard-sweeps",
      title: "Picard sweeps",
      body: [
        "When viscosity depends on the strain rate the problem is nonlinear: "
        + "the flow depends on the viscosity, and the viscosity on the flow. A "
        + "Picard sweep solves for the flow, updates the viscosity from it, and "
        + "solves again.",
        "More sweeps bring the two closer to agreement within each step. "
        + "Each costs a full set of CG iterations, so three sweeps are about "
        + "three times the work of one.",
      ],
      target: "picardSweeps",
    },
    {
      id: "yield-stress",
      title: "Yield stress σ_Y",
      body: [
        "The Tackley and Tosi laws add plastic yielding: above a threshold "
        + "stress the rock fails and flows easily, standing in for the faulting "
        + "that breaks Earth's lithosphere into plates. σ_Y is that threshold "
        + "at the surface.",
        "A lower yield stress lets the stiff cold lid break and sink, giving "
        + "mobile, plate-like surface motion. A higher one keeps the lid "
        + "intact and stagnant.",
      ],
      target: "yieldStress",
    },
    {
      id: "yield-gradient",
      title: "Yield-stress gradient σ_b",
      body: [
        "How quickly the yield stress rises with depth, since rock under more "
        + "pressure is harder to break. The yield stress at depth d is "
        + "σ_Y + σ_b·d.",
      ],
      target: "yieldGradient",
    },
    {
      id: "eta-star",
      title: "Minimum plastic viscosity η*",
      body: [
        "A floor on how weak yielding can make the rock. Without it, the "
        + "viscosity in a fast-deforming zone could fall toward zero, which is "
        + "both unphysical and very hard to solve numerically.",
        "Smaller values give sharper, weaker zones where the lid breaks.",
      ],
      target: "etaStar",
    },
    {
      id: "eta-light",
      title: "η light",
      body: [
        "The van Keken law ignores temperature entirely: viscosity depends "
        + "only on the rock's composition, carried by chemical tracers. η light "
        + "is the viscosity of the light, buoyant material.",
      ],
      target: "etaLight",
    },
    {
      id: "eta-dense",
      title: "η dense",
      body: [
        "The viscosity of the dense material. Where the two are mixed, the "
        + "viscosity lies between the two values in proportion to how much of "
        + "each is present.",
        "Equal values give the benchmark's isoviscous case. A contrast "
        + "between them reproduces its other cases, and changes how fast, and "
        + "in what shape, the dense material overturns.",
      ],
      target: "etaDense",
    },
  ],
  "initial condition": [
    {
      id: "seed-mode",
      title: "Seed mode",
      body: [
        "A perfectly undisturbed layer would never start to convect, so every "
        + "run begins with a small temperature perturbation. This sets its "
        + "shape: how many times the pattern repeats around the annulus or "
        + "across the box.",
        "It is used by restart simulation and seed disturbance. The pattern "
        + "that grows is not always the one seeded: the flow settles on the "
        + "cell size it prefers, and how long that takes is part of what a run "
        + "shows.",
        "The simple view offers the same setting as \"seed pattern\".",
      ],
      target: "seedMode",
    },
  ],
  "view": [
    {
      id: "streamline-density",
      title: "Streamline density",
      body: [
        "How many streamlines are drawn. They are contours of the stream "
        + "function ψ, which the solver computes directly: rock moves along "
        + "them, fastest where they crowd together.",
        "Set it to 0 to turn them off. Like everything in this folder, it "
        + "changes only the picture, never the simulation.",
      ],
      target: "streamlineDensity",
    },
    {
      id: "mesh-overlay",
      title: "Mesh overlay",
      body: [
        "Draws the grid the domain is discretised into (see resolution, under "
        + "domain). There are two: the ψ elements the flow is solved on, and "
        + "the finer T grid that temperature is carried on.",
        "Useful for judging resolution: a boundary layer or plume only a cell "
        + "or two across is under-resolved, and a finer grid would change it.",
      ],
      target: "meshOverlay",
    },
    {
      id: "line-width",
      title: "Line width",
      body: ["The thickness, in pixels, of the streamlines and mesh lines."],
      target: "lineWidth",
    },
    {
      id: "plot-window",
      title: "Plot window",
      body: [
        "How much of the run the two corner plots show, counted in solver "
        + "steps: the Nusselt number (heat transport) and the root-mean-square "
        + "velocity.",
        "A short window shows recent detail. A long one shows whether the run "
        + "has settled into a steady state or is still drifting. Changing it "
        + "only rescales the plots; no recorded history is thrown away.",
      ],
      target: "plotWindow",
    },
  ],
  "tracers": [
    {
      id: "tracer-overlay",
      title: "Tracers",
      body: [
        "Tracers are hundreds of thousands of points carried along by the "
        + "flow, like dye dropped into a fluid. They show where rock has been, "
        + "and how the mantle stirs and mixes, which the temperature field "
        + "alone cannot.",
        "Visual tracers only watch. Chemical tracers also carry a composition, "
        + "dense or light material whose weight pushes back on the flow, for "
        + "problems such as a dense layer at the base of the mantle.",
      ],
      target: "tracerOverlay",
    },
    {
      id: "tracer-count",
      title: "Tracer count",
      body: [
        "More tracers give a smoother, more detailed picture and, for chemical "
        + "tracers, a less noisy composition field, at a cost in speed and "
        + "memory. The middle of the list suits the standard resolution; finer "
        + "grids want more.",
      ],
      target: "tracerCount",
    },
    {
      id: "tracer-colour",
      title: "Colour by",
      body: [
        "What each tracer's colour shows. Initial depth and initial φ colour a "
        + "tracer by where it started, which makes stirring visible as the "
        + "colours are drawn out into filaments. Temperature and speed colour "
        + "it by its present surroundings. Species shows the two materials of a "
        + "chemical run.",
      ],
      target: "tracerColour",
    },
    {
      id: "tracer-size",
      title: "Tracer size",
      body: ["The radius each tracer is drawn with, in screen pixels."],
      target: "tracerSize",
    },
    {
      id: "tracer-opacity",
      title: "Tracer opacity",
      body: [
        "How see-through each tracer is. Lower values show the temperature "
        + "field through the cloud and let dense clusters read as brighter "
        + "patches.",
      ],
      target: "tracerOpacity",
    },
    {
      id: "composition",
      title: "Initial composition",
      body: [
        "How the two materials are laid out when chemical tracers are seeded. "
        + "Dense basal layer: a flat layer of dense material on the core–mantle "
        + "boundary, like the chemically distinct piles thought to sit above "
        + "Earth's core.",
        "van Keken interface: a light layer beneath dense material, with a "
        + "gently curved interface between them. This unstable arrangement is "
        + "the starting point of the van Keken et al. (1997) Rayleigh–Taylor "
        + "benchmark.",
      ],
      target: "composition",
      watch: "Changing it reseeds the tracers.",
    },
    {
      id: "log-rb",
      title: "Compositional Rayleigh number",
      body: [
        "How strongly composition drives the flow: the chemical counterpart "
        + "of the convective-vigour slider. Buoyancy is Ra·T − Rb·C, so the "
        + "larger Rb is, the harder the dense material's extra weight pulls it "
        + "down and holds it there.",
        "Against a given thermal vigour, a large Rb keeps a dense layer as a "
        + "stable blanket at the base. A small one lets plumes entrain it and "
        + "stir it into the mantle.",
      ],
      target: "logRb",
    },
    {
      id: "layer-depth",
      title: "Layer depth",
      body: [
        "The thickness of the dense basal layer, or the height of the van "
        + "Keken interface, as a fraction of the mantle's depth.",
      ],
      target: "layerDepth",
      watch: "Only read when tracers are seeded, so changing it reseeds them.",
    },
    {
      id: "reseed-tracers",
      title: "Reseed tracers",
      body: [
        "Scatters a fresh cloud of tracers at the current settings without "
        + "touching the temperature field, to start the mixing picture again "
        + "from a clean state partway through a run.",
      ],
      target: "reseedTracers",
    },
  ],
} as const satisfies Record<string, readonly SectionHelpStep[]>;

export type SectionName = keyof typeof SECTION_HELP;

/** Anything `tour.ts`'s `start` can open: a guided tour, or one section's help. */
export type WalkthroughName = TourName | SectionName;
