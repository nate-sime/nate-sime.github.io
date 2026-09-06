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
  "preset", "vigour", "rock", "paused", "speed", "flow", "tracers",
  "colormap", "restart", "resetView", "view3d", "advanced",
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
  /**
   * Reconciled against the live `Globe3D.viewMode`, and only toggled when the
   * two differ — so stepping backwards through a tour doesn't flip the view
   * on a step that never asked to change it.
   */
  view?: "2d" | "3d";
  /** Applied through `PaneHandle.applyPatch` — literally the path a preset selection takes. */
  patch?: Partial<State>;
  /** Same, by name, for the entries already in `QUICK_STARTS`/`BENCHMARKS`. */
  preset?: QuickStartName | BenchmarkName;
  /**
   * Drag the convective-vigour slider for the reader over `ms`, rather than
   * snapping `logRa` to its destination: the point of the step is that the
   * picture responds continuously to a control, and a value that teleports
   * shows the ends without the middle.
   */
  ramp?: { to: number; ms: number };
  /** Fly the 2-D camera. Ignored in the 3-D view, which has its own camera. */
  focus?: TourFocus;
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
      title: "A planet's mantle, in cross-section",
      body: [
        "This is a two-dimensional slice through the rocky shell between the "
        + "planet's core and its surface. Colour is temperature: hot at the "
        + "core–mantle boundary below, cold at the surface above.",
        "Rock here is solid, but over millions of years it flows. Hot rock "
        + "rises because it is less dense, cold rock sinks, and the whole "
        + "layer turns over. That is mantle convection, and it is what drives "
        + "plate tectonics.",
      ],
      target: "canvas",
      watch: "The slow overturn is happening now. This is a live simulation, not a recording.",
    },
    {
      id: "flat-view",
      title: "Two ways to look at the same field",
      body: [
        "The globe is a cutaway sphere with this field painted on the cut "
        + "face. This gives perspective as to where the 2D slice sits relative to the planet.",
        "The flat annulus behind it is what the solver actually solves, and "
        + "this is where the rest of the tour will take you.",
      ],
      target: "view3d",
      view: "2d",
    },
    {
      id: "sluggish",
      title: "Just past the onset of convection",
      body: [
        "The examples list loads a ready-made scene. This one is barely "
        + "convecting at all: one or two lazy cells, most of the layer doing "
        + "very little.",
        "Below a critical convection vigour nothing moves. Heat crosses the "
        + "layer by conduction alone, and the picture is a simple temperature "
        + "gradient. Sustained convection occurs at a threshold, not a dial that fades in.",
      ],
      target: "preset",
      preset: "Sluggish mantle",
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
        "This slider sets the Rayleigh number: how hard buoyancy drives the "
        + "flow, against the viscosity and thermal diffusion resisting it. "
        + "Watch it increase by three orders of magnitude.",
        "Earth's mantle sits near the top of this range.",
      ],
      target: "vigour",
      ramp: { to: 6, ms: 2600 },
      dwell: { steps: 300 },
      watch: "The cells break into many thin, fast plumes, and the hot and "
        + "cold layers at the two boundaries grow thinner.",
    },
    {
      id: "rock",
      title: "How the rock behaves",
      body: [
        "So far every part of the layer has been equally stiff. Real rock is "
        + "not: it is far stiffer where it is cold. This law gives the cold "
        + "rock a thousand times the viscosity of the hot rock.",
        "That one change is what turns a field of interchangeable cells into "
        + "a few broad, long-lived upwellings. These upwellings compose "
        + "\"plume\"s.",
      ],
      target: "rock",
      patch: { viscosity: "μ(T, d)", logRa: 5, logContrast: 3, logDepthContrast: 0 },
      dwell: { steps: 500 },
      watch: "A stiff cold lid forming across the top, and plume stems that "
        + "persist instead of drifting apart.",
    },
    {
      id: "boundary-layer",
      title: "The thermal boundary layer",
      body: [
        "Zooming in on the top of the annulus. The thin cold band against the "
        + "outer edge is the thermal boundary layer which composes this model's version of "
        + "a planet's lithosphere.",
        "Nearly the entire temperature drop across the mantle happens inside "
        + "that band. The interior below it is close to uniform in "
        + "temperature, because convection stirs it faster than conduction "
        + "can build a gradient. The layer thickens until it is too heavy, and then falls.",
      ],
      target: "canvas",
      view: "2d",
      // Just inside the outer boundary, at the top of the annulus. The zoom
      // is limited by `clampPan` in main.ts, which keeps at least 1/zoom of
      // the domain on screen — at 5.5 this frames the cold lid together with
      // the upper half of the interior it is falling into, which is the
      // comparison the step is making.
      focus: { zoom: 5.5, x: 0, y: OUTER_RADIUS - 0.26, ms: 1500 },
      dwell: { ms: 6000 },
      watch: "A cold finger thickening, detaching, and sinking away from the surface.",
    },
    {
      id: "flow-lines",
      title: "Flow lines",
      body: [
        "These are streamlines: contours of the stream function, the "
        + "quantity the solver actually solves for. The mantle's velocity "
        + "is recovered from the stream function. A parcel of rock moves along the streamlines.",
        "Closed loops are convection cells. Where the lines crowd together "
        + "the flow is fast. Where they are far apart it is nearly still.",
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
        "The lower panel records the Nusselt number: total heat crossing the "
        + "inner (core-mantle bounary) or outher (lithosphere) boundary, divided by what conduction alone would have carried. "
        + "Nu = 1 means no convection at all. The value here says how much "
        + "the mantle overturn is adding.",
        "Above it, the root-mean-square velocity indicates how fast the layer is "
        + "moving overall. The surface root-mean-square velocity indicates "
        + "the plate-motion velocity. The modern day Earth's Atlantic Ocean spreads at an average rate of 2-5 cm/yr.",
      ],
      target: "traces",
      dwell: { ms: 7000 },
      watch: "A flat trace is a steady state result. A wobbling one is "
        + "time-dependent flow, not noise in the measurement.",
    },
    {
      id: "done",
      title: "That is the tour",
      body: [
        "Everything else is behind \"advanced controls\": viscosity "
        + "laws including yielding and power-law creep, the domain and its "
        + "boundary conditions, the numerics, and published benchmark cases "
        + "from Blankenbach, Tosi and van Keken for solver validation.",
        "Scroll to zoom and drag to pan the canvas at any time, and press H "
        + "to hide the interface.",
        "The run is left exactly where the tour finished, so you can carry "
        + "on from here — or put it back the way you found it.",
      ],
      target: "advanced",
    },
  ],
} as const satisfies Record<string, readonly TourStep[]>;

export type TourName = keyof typeof TOURS;

/** The tour the "Tutorial" chip opens. A second entry above needs a picker; one does not. */
export const DEFAULT_TOUR: TourName = "First look";
