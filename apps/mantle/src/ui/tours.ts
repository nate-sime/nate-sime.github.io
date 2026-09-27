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

import type { PlanetId } from "../planets";
import type { BenchmarkName, QuickStartName, State } from "./presets";

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
  "tutorials", "planet", "preset", "vigour", "seed", "rock", "paused", "speed", "flow", "tracers",
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

export interface TourStep {
  /** Stable across edits to the prose — it is what a step is identified by in tests and in a log. */
  id: string;
  title: string;
  /** One paragraph per entry. */
  body: readonly string[];
  /**
   * The control this step is about, lit while everything else is dimmed.
   * `null` dims the whole screen — the opening and closing cards are about
   * the app, not about any one control.
   */
  target: TourTargetName | null;
  /** An additional control to light with `target` when the experiment needs both. */
  highlight?: TourTargetName;
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
  /** Load a documented planetary profile before applying this step's patch. */
  planet?: PlanetId;
  /**
   * Drag the convective-vigour slider for the reader over `ms`, rather than
   * snapping `logRa` to its destination: the point of the step is that the
   * picture responds continuously to a control, and a value that teleports
   * shows the ends without the middle.
   */
  ramp?: { to: number; ms: number };
  /** Smoothly change the Courant-number control over `ms`. */
  courantRamp?: { to: number; ms: number };
  /** Fly the 2-D camera. Ignored in the 3-D view, which has its own camera. */
  focus?: TourFocus;
  /** Draw a guide along the annulus' outer surface for this step. */
  surfaceGuide?: boolean;
  /** Start this step from the standard small thermal perturbation. */
  reseed?: boolean;
  /** Advance on its own once the effect has had time to develop. "Next" still skips ahead. */
  dwell?: TourDwell;
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
        + "laid flat as a ring: the annulus. Its inner edge is the core–mantle "
        + "boundary and its outer edge is the surface.",
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
        "The \"convective vigour\" slider sets the Rayleigh number, Ra: how "
        + "hard buoyancy drives the flow, against the viscosity and thermal "
        + "diffusion resisting it. The slider is logarithmic, so each step "
        + "along it multiplies Ra rather than adding to it.",
        "Watch it rise from about 2,000 to a million. Earth's mantle is "
        + "estimated at 10⁷–10⁸, higher than this ramp goes.",
      ],
      target: "vigour",
      ramp: { to: 6, ms: 5000 },
      dwell: { steps: 300 },
      watch: "The cells break into many thin, fast upwellings and "
        + "downwellings, and the hot and cold layers along the inner and outer "
        + "edges grow thinner.",
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
      watch: "A stiff cold lid forming along the outer edge, and plume stems "
        + "that persist instead of drifting apart.",
    },
    {
      id: "boundary-layer",
      title: "The thermal boundary layer",
      body: [
        "We've zoomed in on the outer edge of the ring. The thin cold band "
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
        + "for plate speed, and the right-hand axis converts it to cm/yr. For "
        + "comparison, the Atlantic Ocean spreads at about 2–5 cm/yr.",
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
        + "planet tour\" compares Mars, Venus and Earth. Published benchmark "
        + "cases from Blankenbach, Tosi and van Keken are in the \"try an "
        + "example\" list, below the ready-made scenes.",
        "The run is left exactly where the tour finished, so you can carry "
        + "on from here, or use \"restore the run I had\" below to put it back.",
      ],
      target: "tutorials",
    },
  ],
  "Convection onset": [
    {
      id: "conduction",
      title: "Conduction before convection",
      body: [
        "This experiment uses uniform-viscosity rock and a deliberately weak thermal drive (Ra = 100). The highlighted seed-disturbance button adds a small, repeatable temperature pattern, so there is something that could grow.",
        "Instead it fades away. Heat travels from the hot inner boundary to the cold outer boundary by diffusion alone, leaving a smooth conductive temperature gradient.",
      ],
      target: "seed",
      view: "2d",
      patch: { viscosity: "constant", logRa: 2, isothermal: false, paused: false, wavenumber: 2 },
      reseed: true,
      dwell: { steps: 500 },
      watch: "Press seed disturbance whenever you want a fresh test: here it smooths out, Nu approaches 1, and the flow dies away.",
    },
    {
      id: "approach-threshold",
      title: "Approach the critical Rayleigh number",
      body: [
        "Raise the convective-vigour slider slowly, then press seed disturbance to test the new state. Buoyancy strengthens with Rayleigh number, while viscosity and thermal diffusion still erase motion.",
        "Near the onset range, the layer is exceptionally sensitive: a fresh disturbance neither clearly grows nor immediately disappears. Repeat this test as you move the slider.",
      ],
      target: "vigour",
      highlight: "seed",
      dwell: { ms: 7000 },
      watch: "Increase vigour a little, press seed disturbance, and watch whether the new pattern fades or grows.",
    },
    {
      id: "onset",
      title: "Crossing into convection",
      body: [
        "Just above the critical value, buoyancy can amplify a temperature disturbance faster than diffusion removes it. The conductive state becomes unstable and organised circulation appears.",
        "This is a threshold, not merely a gradual increase in activity: below it, conduction is the stable outcome; above it, convection can sustain itself.",
      ],
      target: "vigour",
      highlight: "seed",
      ramp: { to: 3.4, ms: 3500 },
      dwell: { steps: 600 },
      watch: "If the field still looks conductive, press seed disturbance: above onset the new pattern grows into a persistent cell and Nu rises above 1.",
    },
    {
      id: "why-this-value",
      title: "Why does onset happen there?",
      body: [
        "Do not treat this slider value as a universal constant. The critical Rayleigh number depends on the shell's shape, its boundary conditions, viscosity law, and which wavelengths fit in the domain.",
        "Ask why this particular circulation pattern is the first one able to grow. What changes if the layer is wider, the boundaries resist slip, or cold rock becomes more viscous?",
      ],
      target: "vigour",
      dwell: { ms: 7000 },
      watch: "Try moving the slider back and forth across onset, then change one physical assumption and test your prediction.",
    },
  ],
  "Three planet tour": [
    {
      id: "mars-constant-viscosity-warmup",
      title: "Mars warm-up: constant viscosity",
      body: [
        "We begin on Mars with Ra = 10⁵, but keep the viscosity constant while the model warms up. Constant-viscosity flow is computationally cheaper than the more complex rheologies that better represent Mars.",
        "Let the circulation settle into a consistent pattern while the Courant number rises gradually to 2.0. This checks that the numerical experiment is developing cleanly before we add the more expensive temperature dependence.",
      ],
      target: "planet",
      planet: "mars",
      view: "3d",
      patch: {
        resolution: "finest · ψ 192×512", courant: 1.0,
        viscosity: "constant", logContrast: 3, logDepthContrast: 0,
        logRa: 5, wavenumber: 5, isothermal: false, paused: false,
      },
      courantRamp: { to: 2.0, ms: 6000 },
      dwell: { steps: 400 },
      watch: "Once the convection looks consistent, click next to add temperature-dependent viscosity.",
    },
    {
      id: "mars-temperature-dependent-viscosity",
      title: "Mars: a sluggish mantle beneath a rigid lid",
      body: [
        "We begin on Mars with Ra = 10⁵. Its viscosity follows η = exp(−bT), with b = ln(10³): cold material is a thousand times stiffer than hot material.",
        "That cold, stiff outer boundary forms a rigid lid. It resists the motion underneath, so the convection is broad, slow, and sluggish despite the hot material's buoyancy.",
      ],
      target: "rock",
      view: "3d",
      patch: {
        courant: 2.0,
        viscosity: "Blankenbach", logContrast: 3, logDepthContrast: 0,
        logRa: 5, isothermal: false, paused: false,
      },
      dwell: { steps: 400 },
      watch: "Look for a cold, stiff lid at the top and slow circulation beneath it.",
    },
    {
      id: "venus-next",
      title: "Next: the same rheology on Venus",
      body: [
        "Next we will move to Venus. We will keep exactly the same temperature-dependent viscosity law, so the comparison is not caused by changing how the rock responds to temperature.",
        "Venus will be much more vigorous because its material coefficients give it a larger Rayleigh number: buoyancy driving wins more strongly over viscous resistance and thermal diffusion.",
      ],
      target: "planet",
      dwell: { ms: 7000 },
      watch: "Predict what changes when the rheology stays fixed but the buoyancy-to-diffusion balance increases.",
    },
    {
      id: "venus-vigour",
      title: "Venus: more vigorous convection",
      body: [
        "Now Venus is loaded with the same η = exp(−bT) viscosity law. Watch the Rayleigh-number slider rise from 10⁵ to 10⁷.",
        "The extra vigour comes from the coefficients gathered in Ra: stronger buoyancy relative to viscosity and thermal diffusion. The viscosity model itself has not changed; the material balance has.",
      ],
      target: "vigour",
      highlight: "rock",
      planet: "venus",
      patch: { viscosity: "Blankenbach", logContrast: 3, logDepthContrast: 0, logRa: 5, isothermal: false, paused: false },
      ramp: { to: 7, ms: 5000 },
      dwell: { steps: 500 },
      watch: "As Ra rises, plumes multiply and sharpen while the cold lid is stirred more energetically.",
    },
    {
      id: "earth-non-newtonian",
      title: "Next: Earth and non-Newtonian flow",
      body: [
        "Next we will move to Earth and add strain-rate dependence to the temperature-dependent viscosity. This is non-Newtonian flow: the resistance is no longer set by temperature alone, but also by how rapidly the material is deforming.",
        "A familiar analogy is toothpaste: it resists a gentle squeeze, but flows readily where you squeeze it hard. In the mantle model, rapidly deforming regions can likewise become easier to deform than slowly moving ones.",
      ],
      target: "rock",
      highlight: "planet",
      planet: "earth",
      patch: { viscosity: "μ(T, d, ε̇)", logContrast: 3, logDepthContrast: 0, isothermal: false, paused: false },
      dwell: { ms: 8000 },
      watch: "The next experiment keeps temperature sensitivity but lets deformation rate change the rock's effective stiffness.",
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
        "The spherical annulus is a ring-shaped slice through a spherical "
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
