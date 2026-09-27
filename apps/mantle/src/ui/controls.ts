/**
 * Tweakpane controls.
 *
 * Two layers of grouping sit on top of each other here. The one a reader
 * meets first is *audience*: guided tutorials and the planet library open
 * the pane, then the "simulation" folder holds the plain-language controls —
 * try an example, convective vigour, seed disturbance, how the rock behaves,
 * playback, show flow lines, show tracers (and, once that is checked, colour
 * tracers by), temperature colour map, restart simulation, reset view, 3D
 * view, hide UI —
 * and everything else (still every field the solver reads; nothing below is
 * removed) sits in folders hidden behind the "advanced controls" toggle,
 * named for what they let a reader who already knows the physics reach.
 * Checking that toggle is a switch between two complete views, not an
 * addition: the four simple controls that are proxies for an advanced one
 * (how the rock behaves, show flow lines, show tracers, colour tracers by)
 * hide while the advanced folders show, so no setting is on screen twice. The
 * friendly names layered over the technical ones live in `presets.ts`
 * (`QUICK_STARTS`, `SIMPLE_VISCOSITY`) rather than here, for the same reason
 * `BENCHMARKS` does: they are data this file renders, not logic of their
 * own, and worth regression-testing without a DOM.
 *
 * The layer this file has always used, and still uses beneath that, is
 * *cost* — because it is the honest distinction here and it is invisible
 * from the labels:
 *
 *   Ra, contours, line width,   — a 160-byte uniform write; next frame.
 *   mesh, n
 *   σ_Y, σ_b, η* (Tackley,
 *   Tosi)
 *   speed, pause, iterations,   — free; they change how often, or how hard, the
 *   Picard sweeps                 frame loop works, and nothing is precomputed.
 *   Courant number, dt cap      — free here too: dt itself is sized every poll
 *                                 from the GPU's CFL reduction (`adaptiveDt`),
 *                                 so these two only bound that computation
 *                                 rather than triggering it. The f64
 *                                 refactorisation they eventually cause runs
 *                                 in `main.ts`'s frame loop, hysteresis-gated.
 *   reseed                      — rewrites T and re-solves Stokes; one frame.
 *   contrast, depth contrast    — re-invert the μ̄(r) radial blocks in f64, the
 *                                 same job as start-up; announced. Both act on
 *                                 the one profile, for every variable law
 *                                 including Tosi, which reads them as its own
 *                                 γ_T, γ_z.
 *   viscosity tier, resolution, — rebuilds every table and pipeline, 1.3–2.7 s,
 *   geometry, box length          and the page says so rather than appearing to
 *                                 hang. Only *entering or leaving* the Krylov
 *                                 tier does this; μ(T, d) ↔ μ(T, d, ε̇) is a
 *                                 uniform. Tackley and Tosi each use a
 *                                 different pointwise kernel, so entering or
 *                                 leaving *either* is also a rebuild, even
 *                                 though both stay in the tier.
 *                                 Geometry is a rebuild because the metric is
 *                                 compiled into the shaders and the box length
 *                                 reaches the knot vector — see `presets.ts`.
 *
 * This module owns no simulation state: it mutates `State` and calls back. That
 * keeps the pane replaceable (and absent, in tests) without the solver noticing.
 */

import {
  Pane, type ButtonApi, type FolderApi, type ListInputBindingApi,
} from "tweakpane";
import type { BindingApi } from "@tweakpane/core";
import { COLORMAPS, type ColormapName } from "../colormaps";
import { boundaryNames } from "../geometry";
import { isPlanetProfileModified, PLANETS, planetFor, type PlanetId } from "../planets";
import {
  PARTICLE_TINT, SIMPLE_PARTICLE_TINT, SPECIES_CONDITIONS, type TintMode,
} from "../particles";
import { colorbarBlock } from "./colorbar";
import { EQUATION, parseFormula } from "./equation";
import { applyOptgroups, deriveGroups } from "./preset-optgroups";
import type { SectionName, TourName, TourTargetName } from "./tours";
import {
  BENCHMARKS, BOX_LENGTH, CONTRAST, DEPTH_CONTRAST, ETA_VAN_KEKEN, GEOMETRY,
  LABELS, LAYER_DEPTH, LOG_RA, LOG_RB, MESH, NU_WINDOWS, PARTICLE_COUNTS,
  PARTICLE_OPACITY, PARTICLE_SIZE, PARTICLES, PRESETS, QUICK_STARTS, MIN_DT_INITIAL,
  RADIAL_WALLS, SIMPLE_VISCOSITY, SPEEDS, VISCOSITY, type BenchmarkName,
  type CustomSurfaceSource, type GeometryName, type MeshName, type ParticlesName, type PresetName,
  type QuickStartName, type RadialWallsName, type State, type ViscosityName,
  type WallsName,
} from "./presets";

export type {
  BenchmarkName, GeometryName, MeshName, ParticlesName, PresetName,
  RadialWallsName, State, ViscosityName, WallsName,
};
export {
  BENCHMARKS, GEOMETRY, MESH, NU_WINDOWS, PARTICLES, PRESETS, RADIAL_WALLS,
  SPEEDS, VISCOSITY, WALLS, defaultState, geometryFor,
} from "./presets";

export interface Hooks {
  /** Open the guided tour — see `ui/tour.ts`, which `main.ts` builds after this pane and hands back through here. */
  onTutorial(name: TourName): void;
  /** Open one advanced folder's help — the same overlay, over the run as it stands; see `SECTION_HELP` in `tours.ts`. */
  onSectionHelp(name: SectionName): void;
  /** A benchmark case has just written its fields onto `state`; rebuild from it. */
  onBenchmark(): void;
  /** A complete planetary profile has replaced the planet-owned solver fields. */
  onPlanet(id: PlanetId, resumeAfterBuild: boolean): void;
  /** A reader-created annulus changed geometry; existing solver controls remain intact. */
  onCustomPlanet(resumeAfterBuild: boolean): void;
  onRa(v: number): void;
  /** `Ra` forced to 0 regardless of the slider, or released back to it — a pure uniform write either way. */
  onIsothermal(v: boolean): void;
  onStreamlines(levels: number, lineW: number): void;
  onMesh(m: MeshName): void;
  onColormap(v: ColormapName): void;
  onNuWindow(steps: number): void;
  onReseed(): void;
  /** Place a small deterministic temperature disturbance without resetting tracers. */
  onSeedDisturbance(): void;
  onResolution(p: PresetName): void;
  /** Either half of the domain — the list, or the box's length. Both rebuild. */
  onGeometry(): void;
  onViscosity(v: ViscosityName): void;
  /** Either contrast: both re-invert the μ̄(r) blocks, so both take this path. */
  onContrast(log10: number, log10Depth: number): void;
  onIters(n: number): void;
  onExponent(n: number): void;
  onPicard(n: number): void;
  /** Yield-stress parameters, shared by the Tackley and Tosi laws. Pure uniform writes — see `presets.ts`. */
  onSigmaY(v: number): void;
  onSigmaB(v: number): void;
  onEtaStar(v: number): void;
  /** η_light and η_dense together — re-inverts the preconditioner; see `GpuSimulation.setViscosity`. */
  onEtaVanKeken(etaLight: number, etaDense: number): void;
  onResetView(): void;
  /** Start (or reverse) the animated transition to the 3D cutaway-globe view — annulus only, see `Globe3D`. */
  onToggle3D(): void;
  /** Toggle the chrome (HUD, pane, corner plots, scale slider) off for a clean read of the canvas — see `.chrome-hidden` in index.html and `toggleChrome` in main.ts. Also bound to the `H` key there, since this button is one of the things it hides. */
  onToggleChrome(): void;
  /** Toggles the text readout — see `debug` on `State`. */
  onDebug(v: boolean): void;
  /**
   * Tracer overlay mode — see `PARTICLES`. `off ↔` either other value
   * constructs or tears down a `GpuParticles`; `visual ↔ chemical` is one
   * uniform write (the buoyancy coupling), which is why this hook is the one
   * place both costs are decided together rather than split across two.
   */
  onParticles(mode: ParticlesName): void;
  /** Tracer count changed — reallocates the tracer buffers and redraws the cloud at the new count. */
  onParticleCount(): void;
  /** Colour mode changed — rebuilds the push/render pipelines for the new mode's WGSL expression and colour map. */
  onParticleTint(): void;
  /** Dot radius (px) and opacity together — a single uniform write, like `onSigmaY` et al. */
  onParticleStyle(radius: number, opacity: number): void;
  /** The compositional Rayleigh number. A uniform write, read by the buoyancy load only in chemical mode. */
  onRb(v: number): void;
  /** Initial composition profile changed — only ever read at seeding, so this reseeds the cloud with the new profile. */
  onParticleSpecies(): void;
  /** Dense-layer thickness/interface height changed — only ever read at seeding, so this reseeds the cloud with the new profile. */
  onLayerDepth(): void;
  /** Draw a fresh cloud at the current settings, without touching T — `onReseed` ("restart simulation") already redraws the cloud along with T. */
  onReseedParticles(): void;
}

/** Tweakpane list options want `{ label: value }`. */
const nameOptions = <T extends object>(o: T) =>
  Object.fromEntries(Object.keys(o).map((k) => [k, k]));

/**
 * Swap a numeric binding's linear slider for one dragged in log₁₀ space, for
 * the bindings whose range spans several decades (Courant number, dt cap, dt
 * initial) — a linear strip would give the top decade nearly all the travel.
 * There is no `log` option on a plain binding, so Tweakpane's own strip
 * (`.tp-sldv`, inside the `.tp-sldtxtv_s` half of the slider+text view) is
 * hidden and a native `<input type="range">` stands in for it, styled by
 * `.co-log-slider` in index.html to read as one of Tweakpane's own. The number
 * field stays Tweakpane's and keeps reading/accepting real values.
 *
 * Those class names are internal, not public API, so this is only as stable
 * as the Tweakpane version pinned in package.json — and if either is missing
 * it does nothing, leaving the linear slider in place rather than losing it.
 *
 * The binding's `"change"` fires on every drag tick, on committing a typed
 * value and on any refresh that moves it, so the handle follows all three.
 */
function logSlider(
  binding: BindingApi<unknown, number>,
  min: number, max: number, initial: number,
  write: (v: number) => void,
): void {
  const wrap = binding.element.querySelector<HTMLElement>(".tp-sldtxtv_s");
  const linear = wrap?.querySelector<HTMLElement>(".tp-sldv");
  if (!wrap || !linear) return;
  linear.style.display = "none";
  const input = document.createElement("input");
  input.type = "range";
  input.className = "co-log-slider";
  input.min = String(Math.log10(min));
  input.max = String(Math.log10(max));
  input.step = "0.001";
  input.value = String(Math.log10(initial));
  wrap.appendChild(input);
  input.addEventListener("input", () => {
    write(Math.min(max, Math.max(min, 10 ** input.valueAsNumber)));
    binding.refresh();
  });
  binding.on("change", (e) => { input.value = String(Math.log10(e.value)); });
}

const para = (cls: string): HTMLParagraphElement => {
  const e = document.createElement("p");
  e.className = cls;
  return e;
};

/**
 * The equation for the selected law, and a legend naming the slider behind each
 * symbol in it (see `equation.ts` for why that is not decoration). Plain
 * elements rather than a Tweakpane blade: none of the built-in views renders a
 * superscript, and this is read, never edited.
 *
 * `redraw` is called for the law *and* for the sliders whose symbols appear in
 * it, so the legend reads as the handle moves. It follows the *slider*, not the
 * solver: the contrast is applied on release, so mid-drag the γ shown is the one
 * about to be solved with. That is the useful reading — it is what tells you
 * where you are dragging to.
 */
function equationBlock(state: State): { el: HTMLElement; redraw: () => void } {
  const el = document.createElement("div");
  el.className = "eq";

  const redraw = (): void => {
    const eq = EQUATION[state.viscosity];
    el.replaceChildren();

    for (const line of eq.lines) {
      const row = para("eq-f");
      for (const { text, sup } of parseFormula(line)) {
        if (!sup) { row.append(text); continue; }
        const s = document.createElement("sup");
        s.textContent = text;
        row.append(s);
      }
      el.append(row);
    }
    for (const prm of eq.params) {
      const row = para("eq-p");
      const sym = document.createElement("b");
      sym.textContent = `${prm.sym} = ${prm.value(state)}`;
      const from = document.createElement("span");
      from.textContent = prm.control;   // verbatim: it must be findable in the pane
      row.append(sym, from);
      el.append(row);
    }
    if (eq.note) el.append(Object.assign(para("eq-n"), { textContent: eq.note }));
  };

  redraw();
  return { el, redraw };
}

/**
 * The subset of `TOUR_TARGETS` (`tours.ts`) that lives in this pane. The rest
 * — the canvas, the corner traces, the caption — are static elements in
 * index.html that only `main.ts` can resolve, so it is that file which merges
 * the two halves into one resolver.
 */
export type PaneTargetName = Exclude<TourTargetName, "canvas" | "traces" | "caption">;

/**
 * The "?" in an advanced folder's title bar, which opens that folder's help
 * (`SECTION_HELP`, `tours.ts`).
 *
 * A sibling of Tweakpane's title button (`.tp-fldv_b`) rather than a child of
 * it: a button inside a button is invalid HTML, and a click on it would fold
 * the folder too. `.section-help` in index.html positions it over the title
 * bar, just left of the fold mark. `.tp-fldv_b` and `.tp-fldv_c` are internal
 * class names, not public API, so this is only as stable as the Tweakpane
 * version pinned in package.json — the same caveat `logSlider` carries.
 *
 * A collapsed folder's contents are `display: none`, so every control in it
 * would measure as nothing and its help would have nothing to point at. The
 * folder is opened first, and the help waits for the opening to finish: the
 * tour scrolls each control into view as its card opens, and mid-animation
 * the folder is still too short to scroll to. The timeout is the fallback for
 * a height transition that never reports its end.
 */
function sectionHelpButton(folder: FolderApi, open: () => void): void {
  const title = folder.element.querySelector(".tp-fldv_b");
  if (!title) return;
  const b = document.createElement("button");
  b.type = "button";
  b.className = "section-help";
  b.textContent = "?";
  b.title = `What do the ${folder.title} controls do?`;
  b.setAttribute("aria-label", b.title);
  b.addEventListener("click", () => {
    if (folder.expanded) { open(); return; }
    const content = folder.element.querySelector(".tp-fldv_c");
    let fired = false;
    const once = (): void => {
      if (fired) return;
      fired = true;
      content?.removeEventListener("transitionend", onEnd);
      open();
    };
    const onEnd = (e: Event): void => {
      if ((e as TransitionEvent).propertyName === "height") once();
    };
    content?.addEventListener("transitionend", onEnd);
    folder.expanded = true;
    window.setTimeout(once, 400);
  });
  title.after(b);
}

/**
 * Setting a control the way a click on it would, rather than by writing
 * `state` and hoping. Each of these lands in the same function the binding's
 * own `"change"` handler lands in, so a tour driving the pane and a reader
 * driving it are running identical code — the property that makes it safe for
 * `tour.ts` to know nothing about what any given control costs.
 *
 * Only the controls a tour actually *animates* are here; everything else a
 * step wants set goes through `applyPatch` below, which covers every field of
 * `State` at once.
 */
export interface PaneSetters {
  /**
   * Refreshes the one blade rather than the whole pane: `tour.ts`'s vigour
   * ramp calls this every frame for a couple of seconds, and `pane.refresh()`
   * walks every binding in the rack.
   */
  logRa(v: number): void;
  /** Move the Courant control without refreshing the whole pane. */
  courant(v: number): void;
}

/**
 * The pane, plus what has to be reachable from outside it.
 *
 * `view3d` is here because that button's label names the click's destination
 * rather than today's mode, and `Globe3D` — the source of truth for which
 * mode is live — lives outside this module.
 *
 * The rest is what a guided tour needs (`ui/tour.ts`): something to point at,
 * and a way to set a control that goes through this file's own handlers
 * instead of around them.
 */
export interface PaneHandle {
  pane: Pane;
  view3d: ButtonApi;
  /**
   * The rendered element behind each control a tour can point at. Tweakpane
   * renders no IDs and this file's labels change under the "advanced
   * controls" toggle, so a handle captured at build time is the only stable
   * way to find one. `BladeApi.element` is public API, unlike the
   * `.tp-sldtxtv_t`-style reaches elsewhere in this file.
   */
  targets: Record<PaneTargetName, HTMLElement>;
  /** See `applyPatch` below. */
  applyPatch(patch: Partial<State>): void;
  /** Lock solver-mutating blades while retaining the planetary destination selector. */
  setPlanetTraveling(traveling: boolean): void;
  /** Select a supported planetary profile through the same path as the picker. */
  selectPlanet(id: PlanetId): void;
  set: PaneSetters;
}

export function buildPane(state: State, hooks: Hooks): PaneHandle {
  // Mounted into a container the page sizes (see index.html): the app is
  // embedded in an iframe of unknown width, where Tweakpane's default fixed
  // 256px would overlap the readout.
  const pane = new Pane({
    title: "mantle convection",
    container: document.getElementById("pane") ?? undefined,
  });

  // -------------------------------------------------------------------
  // Mode. Every folder below that only a reader who already knows the
  // physics needs is built the same as it always was, then pushed here and
  // hidden until the "advanced controls" checkbox (built right after "reset
  // view", below — see that binding's own note on why *there* and not after
  // the folders) is checked — same ~27 controls as before, just not
  // competing for attention with the seven that actually orient a first
  // visit. Nothing here is deleted; `advancedFolders` is populated as those
  // folders are created further down, and hidden as a batch once the whole
  // pane exists.
  // -------------------------------------------------------------------
  const advancedFolders: FolderApi[] = [];

  // ---- guided tutorials ----
  //
  // First in the rack, above even "try an example": it is the one control
  // here aimed at somebody who does not yet know what any of the others do,
  // and every one of them is easier to find afterwards than before. It walks
  // the pane itself (`ui/tour.ts`) — dimming everything but one control at a
  // time, driving the model, and saying what each does physically — so it
  // belongs with the controls it is about rather than off in the HUD's
  // corner among the links.
  //
  // The one thing that costs: `#pane` is one of the four containers
  // `.chrome-hidden` hides, so this button goes with them under "hide UI".
  // That is the same trade "hide UI" itself already makes by hiding the
  // button that turns it back on, and the tour restores the chrome on its
  // way in regardless.
  //
  // Every tour switches the pane back to its simple view first: the steps
  // point at simple controls ("how the rock behaves", "show flow lines",
  // "show tracers") that the advanced view hides — see `setAdvanced` below.
  const startTutorial = (name: TourName): void => {
    setAdvanced(false);
    hooks.onTutorial(name);
  };
  const tutorials = pane.addFolder({ title: "guided tutorials" });
  tutorials.addButton({ title: "first look" }).on("click", () => startTutorial("First look"));
  tutorials.addButton({ title: "convection onset" }).on("click", () => startTutorial("Convection onset"));
  const tourWarningDialog = document.createElement("dialog");
  tourWarningDialog.className = "pane-dialog tour-warning-dialog";
  tourWarningDialog.setAttribute("aria-labelledby", "tour-warning-title");
  const tourWarningForm = document.createElement("form");
  tourWarningForm.method = "dialog";
  const tourWarningTitle = document.createElement("h2");
  tourWarningTitle.id = "tour-warning-title";
  tourWarningTitle.textContent = "Numerical-resolution warning";
  const tourWarningText = document.createElement("p");
  tourWarningText.textContent = "Accurately representing the models in this tour requires increased numerical resolution.";
  const tourWarningPerformance = document.createElement("p");
  tourWarningPerformance.textContent = "This is computationally expensive and may not perform well on smaller devices, such as phones and laptops without dedicated GPUs.";
  const tourWarningActions = document.createElement("div");
  tourWarningActions.className = "dialog-actions";
  const tourWarningCancel = document.createElement("button");
  tourWarningCancel.type = "button";
  tourWarningCancel.textContent = "cancel";
  tourWarningCancel.addEventListener("click", () => tourWarningDialog.close());
  const tourWarningContinue = document.createElement("button");
  tourWarningContinue.type = "submit";
  tourWarningContinue.textContent = "start tour";
  tourWarningActions.append(tourWarningCancel, tourWarningContinue);
  tourWarningForm.append(tourWarningTitle, tourWarningText, tourWarningPerformance, tourWarningActions);
  tourWarningDialog.append(tourWarningForm);
  document.body.append(tourWarningDialog);
  tourWarningForm.addEventListener("submit", (event) => {
    event.preventDefault();
    tourWarningDialog.close();
    startTutorial("Three planet tour");
  });
  tutorials.addButton({ title: "three planet tour" }).on("click", () => tourWarningDialog.showModal());

  // ---- planet library -----------------------------------------------------
  //
  // A planet is a complete, sourced profile, not another solver. Its selector
  // has a section separate from the one-off examples and live-run controls.
  // The list holds planets only — creating one is an action, so it is a
  // button below the list rather than an entry in it. The numerical-benchmark
  // entry is a placeholder, listed (first, like "— custom —" in "try an
  // example") only while it is the current state: it honestly names what a
  // Cartesian benchmark leaves behind, and is never something to pick.
  const planetLibrary = pane.addFolder({ title: "planet library" });
  const NO_PLANET = "— numerical benchmark —";
  type CustomPlanetChoice = `custom-${number}`;
  type PlanetChoice = PlanetId | CustomPlanetChoice | typeof NO_PLANET;
  const planetState: { planet: PlanetChoice } = {
    planet: state.activePlanet ?? NO_PLANET,
  };
  // Custom bodies live for this page session.  Their full definitions stay in
  // the picker so selecting one later rebuilds exactly the planet that was made.
  const customPlanets = new Map<CustomPlanetChoice, NonNullable<State["customPlanet"]>>();
  let nextCustomPlanet = 1;
  const planetOptions = (): { text: string; value: PlanetChoice }[] => [
    ...(planetState.planet === NO_PLANET ? [{ text: NO_PLANET, value: NO_PLANET } as const] : []),
    ...Object.values(PLANETS).map((planet) => ({ text: planet.label, value: planet.id })),
    ...[...customPlanets].map(([id, planet]) => ({ text: `Custom — ${planet.name}`, value: id })),
  ];
  const planetSelect = planetLibrary.addBinding(planetState, "planet", {
    options: planetOptions(), label: "planet",
  }) as unknown as ListInputBindingApi<PlanetChoice>;
  planetSelect.element.classList.add("planet-selector");
  /** Re-derive the list (the placeholder comes and goes with it) and repaint the selection. */
  const refreshPlanetOptions = (): void => {
    planetSelect.options = planetOptions();
    planetSelect.refresh();
  };
  // `openCreatePlanetDialog` is defined with the dialog below; the click
  // cannot arrive before the whole pane exists.
  planetLibrary.addButton({ title: "create a planet…" }).on("click", () => openCreatePlanetDialog());
  const currentPlanetChoice = (): PlanetChoice =>
    state.activePlanet ?? [...customPlanets].find(([, planet]) => planet === state.customPlanet)?.[0] ?? NO_PLANET;

  const planetDialog = document.createElement("dialog");
  planetDialog.className = "pane-dialog create-planet-dialog";
  planetDialog.setAttribute("aria-labelledby", "create-planet-title");
  const form = document.createElement("form");
  form.method = "dialog";
  const title = document.createElement("h2");
  title.id = "create-planet-title";
  title.textContent = "Create a planet";
  const intro = document.createElement("p");
  intro.textContent = "Set the concentric inner and outer boundaries of the mantle annulus.";
  const nameLabel = document.createElement("label");
  nameLabel.textContent = "name";
  const nameInput = document.createElement("input");
  nameInput.type = "text";
  nameInput.required = true;
  nameInput.maxLength = 64;
  nameInput.value = "Custom planet";
  nameLabel.append(nameInput);
  const field = (label: string, value: string): HTMLInputElement => {
    const row = document.createElement("label");
    row.textContent = label;
    const input = document.createElement("input");
    input.type = "number";
    input.required = true;
    input.min = "1";
    input.step = "1";
    input.value = value;
    row.append(input);
    form.append(row);
    return input;
  };
  form.append(title, intro, nameLabel);
  const innerInput = field("inner radius (km)", "3486");
  const outerInput = field("outer radius (km)", "6371");

  const appearanceTitle = document.createElement("h3");
  appearanceTitle.textContent = "Surface appearance";
  const surfaceChoiceLabel = document.createElement("label");
  surfaceChoiceLabel.textContent = "surface source";
  const surfaceChoice = document.createElement("select");
  const proceduralOption = document.createElement("option");
  proceduralOption.value = "procedural";
  proceduralOption.textContent = "procedural generator";
  const earthOption = new Option("Earth image", "earth-daymap");
  const venusOption = new Option("Venus image (Magellan)", "venus-magellan");
  const marsOption = new Option("Mars image (Viking)", "mars-viking");
  surfaceChoice.append(earthOption, venusOption, marsOption, proceduralOption);
  surfaceChoiceLabel.append(surfaceChoice);
  const generator = document.createElement("fieldset");
  generator.className = "procedural-generator";
  const generatorLegend = document.createElement("legend");
  generatorLegend.textContent = "Procedural generator";
  generator.append(generatorLegend);
  const seedLabel = document.createElement("label");
  seedLabel.textContent = "world seed";
  const seedInput = document.createElement("input");
  seedInput.type = "number";
  seedInput.min = "0";
  seedInput.max = "999999999";
  seedInput.step = "1";
  seedInput.value = "424242";
  seedLabel.append(seedInput);
  generator.append(seedLabel);
  const range = (label: string, value: number): HTMLInputElement => {
    const row = document.createElement("label");
    row.className = "surface-range";
    const text = document.createElement("span");
    text.textContent = label;
    const output = document.createElement("output");
    const input = document.createElement("input");
    input.type = "range";
    input.min = "0";
    input.max = "1";
    input.step = "0.01";
    input.value = String(value);
    const update = () => { output.value = `${Math.round(Number(input.value) * 100)}%`; };
    input.addEventListener("input", update);
    update();
    row.append(text, output, input);
    generator.append(row);
    return input;
  };
  const rockinessInput = range("rockiness", 0.62);
  const terrainScaleInput = range("terrain scale", 0.55);
  const oceanInput = range("ocean coverage", 0.58);
  const plantInput = range("plant life", 0.45);
  const iceInput = range("ice caps", 0.18);
  const cloudInput = range("cloud cover", 0.36);
  const atmosphereHueLabel = document.createElement("label");
  atmosphereHueLabel.textContent = "atmosphere hue";
  const atmosphereHueInput = document.createElement("input");
  atmosphereHueInput.type = "color";
  atmosphereHueInput.value = "#73b8ff";
  atmosphereHueLabel.append(atmosphereHueInput);
  generator.append(atmosphereHueLabel);
  const atmosphereDensityInput = range("atmosphere density", 0.3);
  const appearanceNote = document.createElement("p");
  appearanceNote.className = "appearance-note";
  appearanceNote.textContent = "Earth, Venus, and Mars use the existing credited maps. Procedural settings generate this planet's exterior.";
  const syncSurfaceControls = (): void => {
    const procedural = surfaceChoice.value === "procedural";
    generator.disabled = !procedural;
    generator.hidden = !procedural;
    generator.classList.toggle("is-disabled", !procedural);
  };
  surfaceChoice.addEventListener("change", syncSurfaceControls);
  syncSurfaceControls();
  const actions = document.createElement("div");
  actions.className = "dialog-actions";
  const cancel = document.createElement("button");
  cancel.type = "button";
  cancel.textContent = "cancel";
  cancel.addEventListener("click", () => planetDialog.close());
  const submit = document.createElement("button");
  submit.type = "submit";
  submit.textContent = "create";
  actions.append(cancel, submit);
  form.append(appearanceTitle, surfaceChoiceLabel, generator, appearanceNote, actions);
  planetDialog.append(form);
  document.body.append(planetDialog);
  const openCreatePlanetDialog = (): void => {
    form.reset();
    nameInput.value = "Custom planet";
    innerInput.value = "3486";
    outerInput.value = "6371";
    surfaceChoice.value = "procedural";
    seedInput.value = "424242";
    rockinessInput.value = "0.62";
    terrainScaleInput.value = "0.55";
    oceanInput.value = "0.58";
    plantInput.value = "0.45";
    iceInput.value = "0.18";
    cloudInput.value = "0.36";
    atmosphereHueInput.value = "#73b8ff";
    atmosphereDensityInput.value = "0.3";
    // Repaint the range outputs, which are not form fields themselves.
    for (const input of [rockinessInput, terrainScaleInput, oceanInput, plantInput, iceInput, cloudInput, atmosphereDensityInput])
      input.dispatchEvent(new Event("input"));
    syncSurfaceControls();
    // This intentionally stays non-modal: the UI-scale pane is a global
    // accessibility control and must remain reachable while editing a planet.
    planetDialog.show();
  };
  form.addEventListener("submit", (event) => {
    event.preventDefault();
    const innerRadiusKm = Number(innerInput.value);
    const outerRadiusKm = Number(outerInput.value);
    if (!(innerRadiusKm > 0) || !(outerRadiusKm > innerRadiusKm)) {
      outerInput.setCustomValidity("Outer radius must be greater than the inner radius.");
      outerInput.reportValidity();
      return;
    }
    outerInput.setCustomValidity("");
    const resumeAfterBuild = !state.paused;
    const customPlanet = {
      name: nameInput.value.trim() || "Custom planet", innerRadiusKm, outerRadiusKm,
      surfaceSource: surfaceChoice.value as CustomSurfaceSource,
      surface: {
        kind: surfaceChoice.value as "procedural",
        seed: Math.max(0, Math.floor(Number(seedInput.value) || 0)),
        rockiness: Number(rockinessInput.value),
        terrainScale: Number(terrainScaleInput.value),
        oceanCoverage: Number(oceanInput.value),
        plantLife: Number(plantInput.value),
        iceCaps: Number(iceInput.value),
        cloudCover: Number(cloudInput.value),
        atmosphereHue: atmosphereHueInput.value,
        atmosphereDensity: Number(atmosphereDensityInput.value),
      },
    };
    const customChoice = `custom-${nextCustomPlanet++}` as CustomPlanetChoice;
    customPlanets.set(customChoice, customPlanet);
    state.customPlanet = customPlanet;
    state.activePlanet = null;
    state.geometry = "spherical annulus";
    state.paused = true;
    planetState.planet = customChoice;
    refreshPlanetOptions();
    enableBox(state.geometry);
    refreshBulk();
    planetDialog.close();
    hooks.onCustomPlanet(resumeAfterBuild);
  });

  // The everyday controls are one section of their own: tutorials and the
  // planet library are optional ways into the app, while these are the
  // controls for the live run.
  const simulation = pane.addFolder({ title: "simulation" });

  // ---- try an example: three plain pictures, then the literature ----
  //
  // One dropdown over two tables (`QUICK_STARTS`, `BENCHMARKS`): both are
  // partial `State`s applied the same way, and a reader picking "vigorous
  // convection" shouldn't have to know it lives in a different list than
  // "Blankenbach 1a". Snaps back to "— custom —" immediately after applying,
  // like the benchmark list this replaces always did — a preset is a
  // one-shot load, not a mode the pane keeps asserting once a slider moves.
  const CUSTOM = "— custom —";
  type PresetChoice = QuickStartName | BenchmarkName | typeof CUSTOM;
  const presetTable = { ...QUICK_STARTS, ...BENCHMARKS } as
    Record<QuickStartName | BenchmarkName, Partial<State>>;
  const presetState: { preset: PresetChoice } = { preset: CUSTOM };
  const preset = simulation.addBinding(presetState, "preset", {
    options: {
      [CUSTOM]: CUSTOM, ...nameOptions(QUICK_STARTS), ...nameOptions(BENCHMARKS),
    } as Record<string, PresetChoice>,
    label: "try an example",
  });
  // Purely cosmetic — see `preset-optgroups.ts`'s own header for why this is
  // a separate, isolated reach into Tweakpane's rendered DOM rather than
  // something built into the binding above.
  const presetSelect = preset.element.querySelector("select");
  if (presetSelect) applyOptgroups(presetSelect, deriveGroups(Object.keys(BENCHMARKS)));
  /**
   * Fields of `State` a change to which can only be a rebuild — the metric is
   * compiled into the shaders and the box length reaches the knot vector (see
   * this file's own cost table above). `viscosity` is deliberately *not* here:
   * `hooks.onViscosity` already decides rebuild-versus-uniform for itself, and
   * the rebuild it chooses carries the outgoing law's settled temperature
   * field across rather than reseeding, which is the better answer.
   */
  const REBUILD_KEYS: ReadonlySet<keyof State> = new Set<keyof State>([
    "geometry", "boxLength", "walls", "radialWalls", "resolution",
  ]);

  /**
   * `pane.refresh()` fires the `"change"` handler of every binding whose value
   * it moves — including the five domain bindings behind `REBUILD_KEYS`, each
   * of which rebuilds the solver. A bulk write (a preset, a tour step, a
   * planet) already makes exactly one rebuild of its own, so those five
   * handlers check this flag and stand down while it is set; without it, a
   * single preset that changed geometry, width and walls queued four
   * second-long rebuilds back to back.
   */
  let bulkWrite = false;
  const refreshBulk = (): void => {
    bulkWrite = true;
    try { pane.refresh(); } finally { bulkWrite = false; }
  };

  /**
   * Write a partial `State` onto the live one and make the pane and the
   * solver agree with it — the path a preset selection takes, and the only
   * one `ui/tour.ts` uses to change anything.
   *
   * Geometry, viscosity and their dependent visibility can all move in a
   * single patch, so the housekeeping each of their own change handlers does
   * below has to run here too: those handlers fire from pointer and list
   * events on the pane, none of which a bulk write goes through.
   * `enable`/`eq`/`enableBox`/`enableRa`/`enableParticles`/`pcbar` are all
   * defined further down this function; this only ever runs after the whole
   * pane (and so every one of those consts) exists.
   *
   * `rebuild` is what a preset always wants — it has just written a whole
   * problem statement and `build` reads the lot fresh regardless of which
   * fields moved. A tour patching one or two things wants the opposite: a
   * 1.3–2.7 second rebuild notice in the middle of a step about streamlines
   * would put the sentence explaining them over a blank screen. So the
   * default is to dispatch only the hooks for the keys that actually moved,
   * and to reach for the rebuild only when a key in `REBUILD_KEYS` did.
   */
  const applyPatch = (patch: Partial<State>, rebuild = false): void => {
    const resolutionBefore = state.resolution;
    Object.assign(state, patch);
    adoptState();
    // The same dt-cap adoption a resolution picked by hand gets (see
    // `onResolution` in main.ts), unless the patch states its own cap.
    if (state.resolution !== resolutionBefore && !("dtMax" in patch))
      state.dtMax = PRESETS[state.resolution].dtMax;
    // Literature benchmarks describe numerical domains, not a planet. A
    // quick start leaves the active annulus profile in place and therefore
    // correctly reads as "modified" instead.
    if (patch.geometry === "Cartesian box") {
      state.activePlanet = null;
      state.customPlanet = null;
      planetState.planet = NO_PLANET;
      refreshPlanetOptions();
    }
    enableBox(state.geometry);
    enableRa(state.isothermal);
    enable(state.viscosity);
    enableParticles(state.particles);
    pcbar.setColormap(state.particleColormap);
    eq.redraw();
    syncSimpleControls();
    // Before the hooks below, not after: this fires `"change"` on every simple
    // proxy `syncSimpleControls` just moved (see the note on "show flow
    // lines"), and those handlers call their own hooks. Everything dispatched
    // after it is either a key no proxy covers, or a second, identical call —
    // each of these hooks is a uniform write or an already-guarded no-op.
    // The domain lists' rebuilds are muted for it — see `refreshBulk`.
    refreshBulk();
    if (rebuild || Object.keys(patch).some((k) => REBUILD_KEYS.has(k as keyof State))) {
      hooks.onBenchmark();
      return;
    }
    const has = (k: keyof State): boolean => k in patch;
    // Ra before the isothermal override, so that when a patch moves both, the
    // override is what lands — `onIsothermal` forces Ra to 0 regardless of
    // the slider (see that flag's own header in presets.ts).
    if (has("logRa") && !state.isothermal) hooks.onRa(10 ** state.logRa);
    if (has("isothermal")) hooks.onIsothermal(state.isothermal);
    if (has("viscosity")) hooks.onViscosity(state.viscosity);
    if (has("logContrast") || has("logDepthContrast"))
      hooks.onContrast(state.logContrast, state.logDepthContrast);
    if (has("n")) hooks.onExponent(state.n);
    if (has("picard")) hooks.onPicard(state.picard);
    if (has("iters")) hooks.onIters(state.iters);
    if (has("sigmaY")) hooks.onSigmaY(state.sigmaY);
    if (has("sigmaB")) hooks.onSigmaB(state.sigmaB);
    if (has("etaStar")) hooks.onEtaStar(state.etaStar);
    if (has("etaLight") || has("etaDense"))
      hooks.onEtaVanKeken(state.etaLight, state.etaDense);
    if (has("contours") || has("lineWidth"))
      hooks.onStreamlines(state.contours, state.lineWidth);
    if (has("mesh")) hooks.onMesh(state.mesh);
    if (has("colormap")) hooks.onColormap(state.colormap);
    if (has("nuWindow")) hooks.onNuWindow(state.nuWindow);
    if (has("debug")) hooks.onDebug(state.debug);
    if (has("particles")) hooks.onParticles(state.particles);
    if (has("particleCount")) hooks.onParticleCount();
    if (has("particleTint")) hooks.onParticleTint();
    if (has("particleSpecies")) hooks.onParticleSpecies();
    if (has("layerDepth")) hooks.onLayerDepth();
    if (has("particleSize") || has("particleOpacity"))
      hooks.onParticleStyle(state.particleSize, state.particleOpacity);
    if (has("logRb")) hooks.onRb(10 ** state.logRb);
    // `paused` and `speed` have no hook by design — `main.ts`'s frame loop
    // reads both off `state` directly every frame.
  };

  const selectPlanet = (choice: PlanetChoice): void => {
    if (choice === NO_PLANET) {
      // Only listed while it is already the current state, so there is
      // nothing to switch to — just keep the display honest.
      planetState.planet = currentPlanetChoice();
      refreshPlanetOptions();
      return;
    }
    const customPlanet = customPlanets.get(choice as CustomPlanetChoice);
    if (customPlanet) {
      if (state.customPlanet === customPlanet) {
        planetState.planet = choice;
        refreshPlanetOptions();
        return;
      }
      const resumeAfterBuild = !state.paused;
      state.paused = true;
      state.activePlanet = null;
      state.customPlanet = customPlanet;
      state.geometry = "spherical annulus";
      planetState.planet = choice;
      refreshPlanetOptions();
      enableBox(state.geometry);
      refreshBulk();
      hooks.onCustomPlanet(resumeAfterBuild);
      return;
    }
    const id = choice as PlanetId;
    if (id === state.activePlanet && !isPlanetProfileModified(state)) {
      planetState.planet = id;
      refreshPlanetOptions();
      return;
    }
    const profile = planetFor(id);
    const resumeAfterBuild = !state.paused;
    // Pause synchronously, before the async asset fetch and rebuild hand off
    // to another frame. The main hook restores the reader's prior preference
    // after the new profile is fully live.
    state.paused = true;
    state.activePlanet = id;
    state.customPlanet = null;
    Object.assign(state, profile.solver.state, {
      isothermal: false,
      wavenumber: profile.solver.initialWavenumber,
    });
    // `onPlanet` rebuilds from the whole of `state`; a law the profile
    // changed must not also reach `onViscosity` through the refresh below.
    adoptState();
    planetState.planet = id;
    refreshPlanetOptions();
    enableBox(state.geometry);
    enableRa(state.isothermal);
    enable(state.viscosity);
    eq.redraw();
    syncSimpleControls();
    refreshBulk();
    hooks.onPlanet(id, resumeAfterBuild);
  };
  planetSelect.on("change", (e) => selectPlanet(e.value as PlanetChoice));

  preset.on("change", (e) => {
    const name = e.value;
    if (name === CUSTOM) return;
    // Snapped back before `applyPatch`, so that function's own
    // `pane.refresh()` is the one that repaints the list — a preset is a
    // one-shot load, not a mode the pane keeps asserting once a slider moves,
    // and refreshing twice to say so would be one refresh too many.
    presetState.preset = CUSTOM;
    applyPatch(presetTable[name], true);
  });

  // ---- convective vigour / log₁₀ Ra ----
  //
  // One control, two faces. Simple: "convective vigour", slider only — named
  // for what dragging it does to the picture (more plumes, faster overturn)
  // rather than for what it literally is, and with no number, since a first
  // reader has no unit to read it against. Advanced: "log₁₀ Ra", with the
  // number restored — the label and the format the control had before this
  // pane grew a simple mode at all, for the reader who came to enter an
  // exact value. Both faces share the one binding and the one `logRa`, so
  // there is nothing to keep in sync between them; the "advanced controls"
  // binding below just flips which face is showing, the same switch it
  // throws for every folder it un-hides. Bound to the same `logRa` a linear
  // Ra slider would waste most of its travel on (three decades of
  // interesting behaviour — onset, then plume count) — dragging is
  // log-scale in both faces.
  const vigour = simulation.addBinding(state, "logRa", {
    min: LOG_RA.min, max: LOG_RA.max, step: LOG_RA.step, label: "convective vigour",
  });
  vigour.on("change", (e) => hooks.onRa(10 ** e.value));
  // `.tp-sldtxtv_t` is the number half of the slider+text composite view
  // (`.tp-sldtxtv_s`, alongside it, is the slider half the Courant/dt-cap
  // log-sliders elsewhere in this file already reach into) — an internal
  // class name, not a public API, so it is only as stable as the Tweakpane
  // version pinned in package.json. Hidden by simple default; the "advanced
  // controls" binding below restores it alongside the label.
  const vigourNumber = vigour.element.querySelector<HTMLElement>(".tp-sldtxtv_t");
  if (vigourNumber) vigourNumber.style.display = "none";

  // An exactly conductive numerical field has no non-conductive mode for an
  // instability to amplify. This deliberately does less than "restart
  // simulation": it replaces only T with the standard, reproducible seed,
  // leaving tracers and every physical setting alone for a fair onset test.
  const seed = simulation.addButton({ title: "seed disturbance" });
  seed.on("click", () => hooks.onSeedDisturbance());

  // ---- how the rock behaves ----
  //
  // Three of `VISCOSITY`'s seven laws, under `SIMPLE_VISCOSITY`'s plain
  // names. Bound to its own object rather than `state.viscosity` directly, so
  // its option list can be its own: the four laws this list does not offer
  // (Tackley, Tosi, Blankenbach, van Keken) are reachable under "law" in the
  // advanced viscosity folder, or set by a benchmark or a tour. When one of
  // those is live it is appended to this list under its own name, so the
  // control always shows the law actually being solved rather than the last
  // plain one picked. `applyViscosity` is the one place either list's change
  // lands, so the two can never disagree about what selecting a law costs.
  const isSimpleLaw = (v: ViscosityName): boolean =>
    (Object.values(SIMPLE_VISCOSITY) as ViscosityName[]).includes(v);
  const rockOptions = (v: ViscosityName): { text: string; value: ViscosityName }[] => [
    ...Object.entries(SIMPLE_VISCOSITY).map(([text, value]) => ({ text, value })),
    ...(isSimpleLaw(v) ? [] : [{ text: v, value: v }]),
  ];
  const simpleLaw: { law: ViscosityName } = { law: state.viscosity };
  const rock = simulation.addBinding(simpleLaw, "law", {
    options: rockOptions(state.viscosity), label: "how the rock behaves",
  }) as unknown as ListInputBindingApi<ViscosityName>;
  /** Point the proxy at `state.viscosity`, listing it by name if it is not one of the plain three. Repainted by the caller's `pane.refresh()`. */
  const syncRockOptions = (): void => {
    simpleLaw.law = state.viscosity;
    rock.options = rockOptions(state.viscosity);
  };
  // `enable` and `eq` are defined in the advanced viscosity folder below;
  // referenced here only inside a callback, which never runs before the
  // whole pane (and so both consts) exists.
  //
  // `appliedLaw` drops echoes. `pane.refresh()` fires the "change" handler of
  // whichever list it has just moved to match `state` — the proxy here after
  // a pick under "law", or "law" after a pick here — and re-dispatching
  // `onViscosity` for a law already applied would queue a second rebuild
  // behind the first. `adoptState` (below) moves it for bulk writes, which
  // dispatch their own hooks.
  let appliedLaw = state.viscosity;
  const applyViscosity = (v: ViscosityName): void => {
    if (v === appliedLaw) return;
    appliedLaw = v;
    state.viscosity = v;
    enable(v);
    eq.redraw();
    syncRockOptions();
    hooks.onViscosity(v);
    pane.refresh();
  };
  rock.on("change", (e) => applyViscosity(e.value));

  // ---- playback ----
  //
  // Both held in a const rather than added and forgotten: `PaneHandle.targets`
  // below hands their rendered elements to `ui/tour.ts`, which has no other
  // way to find a blade (Tweakpane renders no IDs).
  const paused = simulation.addBinding(state, "paused");
  // A list rather than a slider: the useful settings span 1/16 to 16 steps per
  // frame, and the labels say what happens far better than a number would.
  const speed = simulation.addBinding(state, "speed", { options: SPEEDS });

  // ---- show flow lines ----
  //
  // Stands in for the density slider (`contours`, 0–60) and the mesh/line-
  // width pair beneath it in the advanced "view" folder: a first visit needs
  // to know streamlines exist, not how many. `24` is an arbitrary but
  // reasonable mid-ladder density — the exact count is exactly what the
  // advanced slider is for.
  const SIMPLE_STREAMLINE_DENSITY = 24;
  const simpleFlow = { on: state.contours > 0 };
  // Preserves a nonzero `contours` a benchmark set rather than snapping it to
  // `SIMPLE_STREAMLINE_DENSITY` — same reason `simpleParticles`'s own "on"
  // handler now guards `state.particles`: `pane.refresh()` fires this "change"
  // event on any programmatic flip of `on`, not only a real click.
  const flow = simulation.addBinding(simpleFlow, "on", { label: "show flow lines" });
  flow.on("change", (e) => {
    state.contours = e.value ? (state.contours > 0 ? state.contours : SIMPLE_STREAMLINE_DENSITY) : 0;
    hooks.onStreamlines(state.contours, state.lineWidth);
    pane.refresh();
  });

  // ---- show tracers ----
  //
  // Off ↔ "visual" — the picture worth a first look. "chemical" (the
  // buoyancy-coupled mode) stays reachable only from the full three-way list
  // in the advanced "tracers" folder: turning tracers off here always lands
  // on "off" outright, the same one-click reset the mockup this was built
  // from settled on, rather than trying to remember which coupled mode to
  // return to.
  //
  // Turning "on" preserves "chemical" rather than collapsing it to "visual".
  // `pane.refresh()` fires this exact "change" event whenever `syncSimpleControls`
  // has just flipped `simpleParticles.on` from false to true to reflect a
  // benchmark's own `Object.assign(state, ...)` — not only on an actual click
  // here — so a two-way `e.value ? "visual" : "off"` would silently zero a
  // benchmark's own `Rb` the moment it finished loading (found reproducing "van
  // Keken 1a": the tracer cloud attached, but the buoyancy load it was meant to
  // drive never coupled, so nothing in the flow ever moved). Turning tracers
  // back *off* is still a one-click reset to "off" outright, same as before.
  const simpleParticles = { on: state.particles !== "off" };
  const tracers = simulation.addBinding(simpleParticles, "on", { label: "show tracers" });
  tracers.on("change", (e) => {
    const mode: ParticlesName = e.value ? (state.particles === "chemical" ? "chemical" : "visual") : "off";
    state.particles = mode;
    enableParticles(mode);
    hooks.onParticles(mode);
    pane.refresh();
  });

  // ---- colour tracers by ----
  //
  // Only worth showing once there is a cloud to colour — hidden until "show
  // tracers" is checked, the same way the advanced folder's own copy of this
  // (`tint`, below) is hidden until `particles` is attached; both read the
  // same condition; see `syncVisibility`. Bound to its own proxy rather than
  // `state.particleTint` directly, the same reason "how the rock behaves" is:
  // `SIMPLE_PARTICLE_TINT` offers two of the full list's seven rows, and the
  // live mode is appended under its own name whenever it is one of the other
  // five, so the control never shows a mode that is not the one drawn.
  // `applyTint` is the one place either list's change lands, so the two can
  // never disagree about what colouring a tracer by X means.
  const isSimpleTint = (t: TintMode): boolean =>
    (Object.values(SIMPLE_PARTICLE_TINT) as TintMode[]).includes(t);
  const tintOptions = (t: TintMode): { text: string; value: TintMode }[] => [
    ...Object.entries(SIMPLE_PARTICLE_TINT).map(([text, value]) => ({ text, value })),
    ...(isSimpleTint(t) ? [] : [{ text: PARTICLE_TINT[t].label, value: t }]),
  ];
  const simpleTintState: { tint: TintMode } = { tint: state.particleTint };
  const simpleTint = simulation.addBinding(simpleTintState, "tint", {
    options: tintOptions(state.particleTint), label: "colour tracers by",
  }) as unknown as ListInputBindingApi<TintMode>;
  /** Point the proxy at `state.particleTint`, listing it by name if it is not one of the plain two. */
  const syncTintOptions = (): void => {
    simpleTintState.tint = state.particleTint;
    simpleTint.options = tintOptions(state.particleTint);
  };
  // Echo guard — see `appliedLaw` above; a repeated tint rebuilds the cloud.
  let appliedTint = state.particleTint;
  const applyTint = (t: TintMode): void => {
    if (t === appliedTint) return;
    appliedTint = t;
    state.particleTint = t;
    state.particleColormap = PARTICLE_TINT[t].colormap;
    pcbar.setColormap(state.particleColormap);
    syncTintOptions();
    hooks.onParticleTint();
    pane.refresh();
  };
  simpleTint.on("change", (e) => applyTint(e.value));

  // ---- colour map ----
  //
  // Already plain — a swatch, not a formula — so it stays at the root next
  // to the legend it always sat beside.
  const cbar = colorbarBlock(state.colormap, ["0 (cold)", "1 (hot)"]);
  const cmap = simulation.addBinding(state, "colormap",
    { options: nameOptions(COLORMAPS), label: "temperature colour map" });
  cmap.on("change", (e) => {
    cbar.setColormap(e.value as ColormapName);
    hooks.onColormap(e.value as ColormapName);
  });

  // ---- restart simulation ----
  //
  // Re-seeds T, resets the clock and re-solves Stokes from it, and — if a
  // tracer cloud is attached — redraws that too (`GpuSimulation.reseed` does
  // both), so every field the picture shows starts over together. What it
  // restarts *into* (seed mode, composition, species) is set from the
  // advanced "initial condition" and "tracers" folders; this is the one
  // place in the pane that restarts the whole run. Compare "seed
  // disturbance", which replaces T but leaves the tracers where they are.
  const restart = simulation.addButton({ title: "restart simulation" });
  restart.on("click", () => hooks.onReseed());

  // Scroll to zoom, drag to pan (see main.ts) — this is the way back from
  // either with no pointer precision required.
  const resetView = simulation.addButton({ title: "reset view" });
  resetView.on("click", () => hooks.onResetView());

  // The wedge cutaway is a cross-section of a sphere, which a Cartesian box
  // has no version of — disabled there rather than hidden, the same policy
  // `len`/`walls` follow below, so picking a geometry never moves this
  // button out from under the pointer. The label names the *destination*,
  // not the current state ("3D view" while flat, "scientific view" once in
  // the globe) — `main.ts` owns `Globe3D`, the source of truth for which
  // mode is live, so it writes this title directly through the returned
  // `view3d` handle rather than this module tracking a copy of state it
  // cannot itself observe.
  const view3d = simulation.addButton({ title: "3D view" })
    .on("click", () => hooks.onToggle3D());

  // ---- hide UI ----
  //
  // A clean read of the canvas alone — the flow field, or the globe, with
  // nothing drawn over it — for a screenshot or a recording. The one place
  // this button ends up once clicked is nowhere: `#pane` is one of the four
  // containers `.chrome-hidden` (index.html) hides, so the button that turns
  // the chrome off disappears along with it. `main.ts` binds the `H` key to
  // the same toggle as the way back in, and shows a brief on-canvas reminder
  // of it the moment the chrome disappears — the only affordance left once
  // this fires.
  simulation.addButton({ title: "hide UI (H)" }).on("click", () => hooks.onToggleChrome());

  // ---- UI scale ----
  //
  // Resizes the pane, HUD and corner plots together, via the `--ui-scale`
  // custom property `index.html` defines each of those three around — the
  // app's chrome, not the simulation canvas, which already has its own
  // scroll-to-zoom (see main.ts).
  //
  // Its own small, title-less Tweakpane instance in `#scale` (index.html) —
  // deliberately *not* a binding on `pane` itself. `pane`'s own element is
  // one of the three `--ui-scale` resizes, and this control sets that
  // variable on every drag tick, so mounting it there would have it resize
  // its own drag surface out from under the pointer mid-gesture — a small,
  // deliberate movement turning into a runaway jump toward whichever end
  // the resize was already leaning (see `#scale`'s own note in index.html).
  // A second `Pane` gets Tweakpane's own slider/number/keyboard handling for
  // free, styled identically to the first, rather than a bespoke widget.
  const scalePane = new Pane({ container: document.getElementById("scale") ?? undefined });
  scalePane.addBinding(state, "uiScale", {
    min: 0.75, max: 1.75, step: 0.05, label: "UI scale",
  }).on("change", (e) => {
    document.documentElement.style.setProperty("--ui-scale", String(e.value));
  });

  // ---- advanced controls ----
  //
  // Built here, right after the last plain-language control and before any
  // advanced folder — not after them — so this stays put in the rack
  // regardless of whether it is checked. Placed *after* the folders (as
  // "built last" once was), checking it un-hides several screens of content
  // that were sitting between this checkbox and the controls above it, which
  // shoves the checkbox itself far down the pane the instant it is clicked —
  // the opposite of what a fixed anchor is for. Here, every advanced folder
  // renders *below* this line whether hidden or not, so this is always the
  // last thing in the "simulation" folder.
  const ui = { advanced: false };
  const advanced = simulation.addBinding(ui, "advanced", { label: "advanced controls" });
  advanced.on("change", () => syncVisibility());
  /** Switch views from outside the checkbox — the tutorials use it to start from the simple view. */
  const setAdvanced = (on: boolean): void => {
    if (ui.advanced === on) return;
    ui.advanced = on;
    advanced.refresh();
    syncVisibility();
  };
  /**
   * The one place that decides what the pane shows in each view. Advanced:
   * every advanced folder, and none of the four simple proxies — each has
   * its full-range counterpart in one of those folders ("law", "streamline
   * density", "tracer overlay", "colour by"), so showing both would put one
   * setting on screen twice. Simple: the reverse. The vigour slider is the
   * one control in both views; it swaps faces instead (see its own note).
   * "colour tracers by" additionally needs a cloud to colour.
   */
  const syncVisibility = (): void => {
    const adv = ui.advanced;
    for (const f of advancedFolders) f.hidden = !adv;
    rock.hidden = flow.hidden = tracers.hidden = adv;
    simpleTint.hidden = adv || !PARTICLES[state.particles].attached;
    vigour.label = adv ? "log₁₀ Ra" : "convective vigour";
    if (vigourNumber) vigourNumber.style.display = adv ? "" : "none";
  };

  /**
   * Re-reads `state` into the four simple proxies above, for whichever of
   * them a change made elsewhere (a preset, or an advanced control) may have
   * moved without going through the simple control itself. Cheap and always
   * safe to call — each line is a no-op unless the two actually disagree.
   */
  const syncSimpleControls = (): void => {
    syncRockOptions();
    simpleFlow.on = state.contours > 0;
    simpleParticles.on = state.particles !== "off";
    syncTintOptions();
    syncVisibility();
  };

  /**
   * Bulk writes (`applyPatch`, a planet selection) dispatch their own hooks,
   * so the echo guards on the two list proxies are moved up to the new
   * values first — otherwise the `pane.refresh()` that follows would have
   * those lists' handlers dispatch the same hook a second time.
   */
  const adoptState = (): void => {
    appliedLaw = state.viscosity;
    appliedTint = state.particleTint;
  };

  // =====================================================================
  // Advanced. Folders below are built exactly as the pane has always built
  // them and hidden as a batch at the end of this function — see `ui`
  // above.
  // =====================================================================

  // *What* is being solved, above everything about how. Both controls in here
  // rebuild every table and pipeline (see `presets.ts`), so both are announced;
  // the length is disabled rather than hidden on the annulus, so selecting a
  // geometry does not move the rest of the pane out from under the pointer.
  const dom = pane.addFolder({ title: "domain" });
  advancedFolders.push(dom);
  const geom = dom.addBinding(state, "geometry",
    { options: nameOptions(GEOMETRY), label: "geometry" });
  const len = dom.addBinding(state, "boxLength", {
    min: BOX_LENGTH.min, max: BOX_LENGTH.max, step: BOX_LENGTH.step,
    format: BOX_LENGTH.format, label: "box width",
  });
  // Below the width, because it is a statement about the domain's *edges* and
  // reads as one only once there is a width for them to be the edges of.
  // Displayed without the `WALLS` key's "walls" suffix, so the two boundary
  // lists read in parallel ("free-slip" here, "free-slip"/"no-slip" below);
  // the key itself stays, since benchmarks and saved state name it.
  const walls = dom.addBinding(state, "walls", {
    options: { "periodic": "periodic", "free-slip": "free-slip walls" } satisfies Record<string, WallsName>,
    label: "left / right",
  });
  // Unlike `walls`, legal on *both* geometries — a no-slip radial condition
  // means the same thing on an annulus (inner/outer) as on a box (top/
  // bottom), so this is never disabled the way `len`/`walls` are, only
  // relabelled to say which pair of boundaries it closes. See
  // `boundaryNames` in `geometry.ts`, the same names the Nusselt readout uses.
  const radialWalls = dom.addBinding(state, "radialWalls",
    { options: nameOptions(RADIAL_WALLS), label: "top / bottom" });
  const enableBox = (g: GeometryName) => {
    len.disabled = walls.disabled = GEOMETRY[g] !== "box";
    view3d.disabled = GEOMETRY[g] !== "annulus";
    const bn = boundaryNames(GEOMETRY[g]);
    radialWalls.label = `${bn.inner} / ${bn.outer}`;
  };
  // Every rebuild below stands down during a bulk write, which makes its own
  // — see `refreshBulk`.
  geom.on("change", (e) => {
    enableBox(e.value as GeometryName);
    if (!bulkWrite) hooks.onGeometry();
  });
  // On release only. The width changes the azimuthal knot vector, so it is the
  // same second-or-two rebuild the resolution list is; firing it per pointer
  // move would queue one for every pixel dragged. The list has no drag to wait
  // for, so it fires on change like every other list in the pane.
  len.on("change", (e) => { if (e.last && !bulkWrite) hooks.onGeometry(); });
  walls.on("change", () => { if (!bulkWrite) hooks.onGeometry(); });
  radialWalls.on("change", () => { if (!bulkWrite) hooks.onGeometry(); });
  enableBox(state.geometry);
  const resolution = dom.addBinding(state, "resolution",
    { options: nameOptions(PRESETS), label: "resolution" });
  resolution.on("change", (e) => {
    if (bulkWrite) return;
    hooks.onResolution(e.value as PresetName);
    // `onResolution` adopts the preset's own dt cap onto `state`; show it.
    // `dtMax` is built in the numerics folder below — a click cannot
    // arrive before it exists.
    dtMax.refresh();
  });

  // What actually drives the step, plus the ceiling it is held under, and the
  // isothermal override — everything about the solve that isn't "how vigorous"
  // or "which law", both of which live in the simulation folder above.
  const numerics = pane.addFolder({ title: "numerics" });
  advancedFolders.push(numerics);
  // Forces Ra = 0 regardless of the convection-vigour slider above — the
  // purely compositional (isothermal) buoyancy the van Keken Rayleigh–Taylor
  // benchmark needs (see `isothermal`'s own header in presets.ts on why this
  // is a checkbox rather than a widened `logRa` floor). Disables the vigour
  // slider while checked: `logRa`'s value is not what is being solved with,
  // and a slider that is still draggable but silently ignored would be worse.
  // Disabled rather than hidden, the pane's policy for a control that does
  // not apply to the current setup (box width, side walls, 3D view), so the
  // rack does not shift under the pointer and a tour can still point at it.
  const iso = numerics.addBinding(state, "isothermal", { label: "isothermal (Ra = 0)" });
  const enableRa = (isothermal: boolean): void => { vigour.disabled = isothermal; };
  iso.on("change", (e) => { enableRa(e.value); hooks.onIsothermal(e.value); });
  enableRa(state.isothermal);
  // 0.1–100: three decades, so the number field alone would need three
  // regimes of care from the reader — fine near 0.1, coarse near 100 — while
  // showing the same digit count throughout. `format` gives each decade one
  // more decimal than the one above it instead. `step` is a granularity floor
  // (Tweakpane snaps the bound value to its nearest multiple, including on
  // programmatic writes — see the slider below), not an editing increment: at
  // 0.001 it is finer than the display ever shows, so it never visibly bites.
  const courant = numerics.addBinding(state, "courant", {
    min: 0.1, max: 100, step: 0.001, label: "Courant number",
    format: (v) => v.toFixed(v < 1 ? 3 : v < 10 ? 2 : 1),
  });
  // Co ≤ 1 is the conventional, dt-limited-by-nothing-but-accuracy regime;
  // above 1 the step is coarser than one cell crossing per step, which is
  // still fine here (SL advection + implicit diffusion are unconditionally
  // stable — see `gpu/sim.ts`) but trades accuracy for it, more so past 3.
  // Tweakpane has no per-value text colour on a binding, so this reaches into
  // its own DOM: `.tp-txtv_i` is the number field inside the combined
  // slider+text view a `min`/`max`/`step` binding renders as — an internal
  // class name, not a public API, so it is only as stable as the Tweakpane
  // version pinned in package.json.
  const courantInput = courant.element.querySelector<HTMLInputElement>(".tp-txtv_i");
  const courantColour = (v: number): string =>
    v > 3 ? "#ff5c5c" : v > 1 ? "#e8a33d" : "#ffffff";
  const paintCourant = (v: number): void => {
    if (courantInput) courantInput.style.color = courantColour(v);
  };
  // Tweakpane's own slider is linear in the bound value, which across three
  // decades would put all the usable travel in the top decade and leave 0.1–1
  // a couple of pixels wide — see `logSlider` (top of this file), which the
  // two dt bindings below share.
  logSlider(courant, 0.1, 100, state.courant, (v) => { state.courant = v; });
  courant.on("change", (e) => paintCourant(e.value));
  paintCourant(state.courant);
  // 1e-4 to 1e3: seven decades, wider even than Courant's three above — a
  // run's own accuracy ceiling can sit orders of magnitude above the
  // resolution ladder's default, so the slider has to reach past it — so a
  // linear slider is even less usable here than Courant's was. `format`'s
  // fixed-decimal digit count would be unreadable across that range too, so
  // this reads in the same scientific notation the codebase's own comments
  // state these ceilings in.
  const dtMax = numerics.addBinding(state, "dtMax", {
    min: 1e-4, max: 1e3, step: 1e-6, label: "dt cap",
    format: (v) => v.toExponential(1),
  });
  logSlider(dtMax, 1e-4, 1e3, state.dtMax, (v) => { state.dtMax = v; });
  // 1e-6 to 1e3: the step `GpuSimulation.create` is seeded with, before the
  // first `pollStats` readback gives `adaptiveDt` a CFL-implied value to work
  // from (see `dtInitial` on `State`). Only read at build time — unlike the
  // cap immediately above, changing it has no effect on a solver already
  // running, only the next one built — which the label says, since nothing
  // else in the pane would.
  const dtInitial = numerics.addBinding(state, "dtInitial", {
    min: MIN_DT_INITIAL, max: 1e3, step: MIN_DT_INITIAL, label: "dt initial (on rebuild)",
    format: (v) => v.toExponential(1),
  });
  logSlider(dtInitial, MIN_DT_INITIAL, 1e3, state.dtInitial, (v) => { state.dtInitial = v; });

  // Viscosity: the full law list (all seven — the three the simple "how the
  // rock behaves" control offers, plus Tackley, Tosi, Blankenbach and van
  // Keken under the names their own papers use) picks the rheology, and the
  // knobs below it only mean anything for some of them — so they are hidden
  // rather than disabled. With six laws' worth of knobs in play, greying out
  // the ones that don't apply still leaves them taking up space and
  // competing for attention; hiding them is what actually keeps the pane
  // readable. Two levels of that: contrast and the CG budget need the Krylov
  // tier, n and the Picard sweeps need the power law on top of it.
  const rheo = pane.addFolder({ title: "viscosity" });
  advancedFolders.push(rheo);
  const law = rheo.addBinding(state, "viscosity",
    { options: nameOptions(VISCOSITY), label: "law" });
  const eq = equationBlock(state);
  // Bounds and step are `CONTRAST`/`DEPTH_CONTRAST` in presets.ts — see there
  // for why the step is as fine as it is.
  const contrast = rheo.addBinding(state, "logContrast",
    { ...CONTRAST, label: LABELS.contrast });
  // Directly below the thermal contrast, because they are the same kind of
  // number — a log₁₀ ratio across the layer — and reading them as a pair is what
  // says the total contrast is their product. Its floor is 0 (no depth
  // dependence, the law the app opens with) and its ceiling is lower than the
  // thermal one's: the two multiply inside one clamp, and 10⁵ of each is a
  // contrast no fixed Krylov budget is going to hold.
  const depth = rheo.addBinding(state, "logDepthContrast",
    { ...DEPTH_CONTRAST, label: LABELS.depth });
  const nExp = rheo.addBinding(state, "n",
    { min: 1, max: 5, step: 0.25, label: LABELS.n });
  const iters = rheo.addBinding(state, "iters",
    { min: 1, max: 40, step: 1, label: "CG iterations" });
  const picard = rheo.addBinding(state, "picard",
    { min: 1, max: 3, step: 1, label: "Picard sweeps" });
  // Tackley's own parameters — γ, c and n mean nothing to Tackley, so it gets
  // its own three rather than reusing contrast/depth/n under a different name.
  // Tosi states the identical yielding branch, so it reuses these three
  // rather than getting a second copy under different names.
  const sigmaY = rheo.addBinding(state, "sigmaY",
    { min: 0, max: 5, step: 0.1, label: LABELS.sigmaY });
  const sigmaB = rheo.addBinding(state, "sigmaB",
    { min: 0, max: 5, step: 0.1, label: LABELS.sigmaB });
  const etaStar = rheo.addBinding(state, "etaStar",
    { min: 1e-4, max: 1e-2, step: 1e-4, label: LABELS.etaStar });
  // van Keken's own two parameters — no T dependence at all, so γ/c mean
  // nothing to it either, the same reason Tackley gets its own set instead
  // of reusing contrast/depth.
  const etaLight = rheo.addBinding(state, "etaLight",
    { ...ETA_VAN_KEKEN, label: LABELS.etaLight });
  const etaDense = rheo.addBinding(state, "etaDense",
    { ...ETA_VAN_KEKEN, label: LABELS.etaDense });

  const enable = (v: ViscosityName) => {
    const { variable, strainRate, tackley, tosi, vanKeken } = VISCOSITY[v];
    contrast.hidden = depth.hidden = !variable || tackley || vanKeken;
    iters.hidden = !variable;
    nExp.hidden = !strainRate || tackley || tosi;
    picard.hidden = !strainRate;
    sigmaY.hidden = sigmaB.hidden = etaStar.hidden = !(tackley || tosi);
    etaLight.hidden = etaDense.hidden = !vanKeken;
  };
  // Both the advanced list and the simple "how the rock behaves" control
  // above land on this one function, so neither can apply a law the other
  // doesn't also know about.
  law.on("change", (e) => applyViscosity(e.value as ViscosityName));
  // Both contrasts re-invert the preconditioner in f64, so they fire on release
  // rather than while dragging, and each sends both values: the rebuild is one
  // job over μ̄(r), which is a function of γ *and* c, so there is nothing for a
  // per-slider callback to do differently. n is a plain uniform and the two
  // counts are only loop bounds, so those take effect as they are dragged. The
  // equation's symbols follow the slider either way — redrawing costs nothing,
  // and a legend that lagged the handle would be worse than none.
  const applyContrast = (e: { last: boolean }) => {
    eq.redraw();
    if (e.last) hooks.onContrast(state.logContrast, state.logDepthContrast);
  };
  contrast.on("change", applyContrast);
  depth.on("change", applyContrast);
  nExp.on("change", (e) => { eq.redraw(); hooks.onExponent(e.value); });
  iters.on("change", (e) => hooks.onIters(e.value));
  picard.on("change", (e) => hooks.onPicard(e.value));
  sigmaY.on("change", (e) => { eq.redraw(); hooks.onSigmaY(e.value); });
  sigmaB.on("change", (e) => { eq.redraw(); hooks.onSigmaB(e.value); });
  etaStar.on("change", (e) => { eq.redraw(); hooks.onEtaStar(e.value); });
  // On release, like the two contrasts: both feed μ̄(r), so both re-invert
  // the preconditioner in f64 — see `GpuSimulation.setViscosity`.
  const applyVanKekenViscosity = (e: { last: boolean }) => {
    eq.redraw();
    if (e.last) hooks.onEtaVanKeken(state.etaLight, state.etaDense);
  };
  etaLight.on("change", applyVanKekenViscosity);
  etaDense.on("change", applyVanKekenViscosity);
  enable(state.viscosity);

  // Initial condition: which seed mode a fresh run starts from. Split out
  // from the numerics folder above because it is a decision about the
  // *starting picture*, not about how accurately the solve tracks it once
  // running. Read by "restart simulation" and "seed disturbance" in the
  // simulation folder — no button of its own here, since a second restart
  // button would do exactly what that one does.
  const initial = pane.addFolder({ title: "initial condition" });
  advancedFolders.push(initial);
  const seedMode = initial.addBinding(state, "wavenumber",
    { min: 1, max: 12, step: 1, label: "seed mode" });

  // The streamline density, the mesh overlay and the line width both draw
  // with — colour map lives in the simulation folder, next to the legend it
  // has always sat beside, so it is not repeated here.
  const view = pane.addFolder({ title: "view" });
  advancedFolders.push(view);
  const density = view.addBinding(state, "contours",
    { min: 0, max: 60, step: 2, label: "streamline density" });
  density.on("change", (e) => {
    hooks.onStreamlines(e.value, state.lineWidth);
    // The simple "show flow lines" switch reads this same field, so a
    // density dragged to (or off) zero here has to be reflected there too.
    simpleFlow.on = e.value > 0;
    pane.refresh();
  });
  const mesh = view.addBinding(state, "mesh", { options: nameOptions(MESH), label: "mesh overlay" });
  mesh.on("change", (e) => hooks.onMesh(e.value as MeshName));
  const lineWidth = view.addBinding(state, "lineWidth",
    { min: 0.5, max: 3, step: 0.1, label: "line width" });
  lineWidth.on("change", (e) => hooks.onStreamlines(state.contours, e.value));
  // How much of the run the two corner plots show — Nusselt number and RMS
  // velocity share this one control (see `presets.ts`). Costs nothing: both
  // traces keep every sample either way, so this re-scales an existing
  // buffer and does not begin collecting again.
  const plotWindow = view.addBinding(state, "nuWindow", { options: NU_WINDOWS, label: "plot window" });
  plotWindow.on("change", (e) => hooks.onNuWindow(e.value));

  // The tracer overlay in full: the three-way mode list (the simple "show
  // tracers" switch only ever reaches "off" and "visual" — "chemical" lives
  // here), plus every control that only means something once a cloud
  // exists. Structured like the viscosity folder above: one list decides
  // which controls beneath it mean anything, and those are hidden rather
  // than disabled. "Tracer" throughout, the word the simple controls use.
  const trace = pane.addFolder({ title: "tracers" });
  advancedFolders.push(trace);
  const mode = trace.addBinding(state, "particles",
    { options: nameOptions(PARTICLES), label: "tracer overlay" });
  const count = trace.addBinding(state, "particleCount",
    { options: PARTICLE_COUNTS, label: "tracer count" });
  const tint = trace.addBinding(state, "particleTint",
    { options: nameOptions(PARTICLE_TINT), label: "colour by" });
  const pcbar = colorbarBlock(state.particleColormap);
  const size = trace.addBinding(state, "particleSize", {
    min: PARTICLE_SIZE.min, max: PARTICLE_SIZE.max, step: PARTICLE_SIZE.step,
    label: "tracer size",
  });
  const opacity = trace.addBinding(state, "particleOpacity", {
    min: PARTICLE_OPACITY.min, max: PARTICLE_OPACITY.max, step: PARTICLE_OPACITY.step,
    label: "tracer opacity",
  });
  // Chemical-only: the initial composition profile means nothing to a
  // purely visual cloud, Rb has no effect on one (`Rb` stays at 0 regardless
  // of the slider — see `PARTICLES`), and the dense-layer thickness/interface
  // height is only ever read while a fresh cloud is being seeded, which
  // nothing but the chemical mode's own initial composition consumes.
  const species = trace.addBinding(state, "particleSpecies",
    { options: nameOptions(SPECIES_CONDITIONS), label: "composition" });
  const rb = trace.addBinding(state, "logRb", { ...LOG_RB, label: "log₁₀ Rb" });
  const layer = trace.addBinding(state, "layerDepth", {
    min: LAYER_DEPTH.min, max: LAYER_DEPTH.max, step: LAYER_DEPTH.step,
    label: "layer depth",
  });
  const reseedTracers = trace.addButton({ title: "reseed tracers" });
  reseedTracers.on("click", () => hooks.onReseedParticles());

  const enableParticles = (m: ParticlesName): void => {
    const { attached, coupled } = PARTICLES[m];
    count.hidden = tint.hidden = size.hidden = opacity.hidden = !attached;
    pcbar.el.hidden = !attached;
    species.hidden = rb.hidden = layer.hidden = !coupled;
    // The simple "colour tracers by" reads the same condition as this
    // folder's `tint`, combined with which view is showing.
    syncVisibility();
  };
  mode.on("change", (e) => {
    const m = e.value as ParticlesName;
    enableParticles(m);
    // The simple "show tracers" switch reads this same field as an on/off,
    // so a mode picked here — including "chemical", which that switch never
    // reaches on its own — has to be reflected there too.
    simpleParticles.on = m !== "off";
    hooks.onParticles(m);
    pane.refresh();
  });
  count.on("change", () => hooks.onParticleCount());
  // Both the advanced list and the simple "colour tracers by" control above
  // land on this one function, so neither can apply a mode the other
  // doesn't also know about.
  tint.on("change", (e) => applyTint(e.value as TintMode));
  const applyStyle = (): void => hooks.onParticleStyle(state.particleSize, state.particleOpacity);
  size.on("change", applyStyle);
  opacity.on("change", applyStyle);
  // A pure uniform write (see `Rb`'s own header in `gpu/wgsl.ts`), so this
  // takes effect while dragging, like Ra.
  rb.on("change", (e) => hooks.onRb(10 ** e.value));
  // On release only, like the box width and the two contrasts: both reseed
  // the whole cloud with a differently painted initial composition, not a
  // value the running push kernel reads every step.
  species.on("change", () => hooks.onParticleSpecies());
  layer.on("change", (e) => { if (e.last) hooks.onLayerDepth(); });
  enableParticles(state.particles);

  // The one control here that is about the *page* rather than the physics —
  // last, since a first-time reader has the least use for it.
  const dbg = pane.addFolder({ title: "debug" });
  advancedFolders.push(dbg);
  dbg.addBinding(state, "debug", { label: "debug mode" })
    .on("change", (e) => hooks.onDebug(e.value));

  // Simple by default: every folder just built starts hidden, and the
  // "advanced controls" binding above flips the two views together.
  syncVisibility();

  // The equation goes between the law and the knobs, because that is what it
  // connects: the law the list just selected, and the sliders below named
  // against the symbols they set. It is inserted *last* — a rack re-appends its
  // blades' elements as they are added, so a foreign node placed mid-folder
  // drifts to the bottom of it as the rest of the folder is built. The colour
  // bar is the same trick for the same reason: it belongs right after the
  // colour-map list, and `trace` gains several more blades after `tint`.
  law.element.after(eq.el);
  cmap.element.after(cbar.el);
  tint.element.after(pcbar.el);

  // A "?" on every advanced folder but debug, which is about the page rather
  // than the physics. Opened without `startTutorial`'s `setAdvanced(false)`:
  // the help is about the advanced view, and would hide its own folder.
  const sections: [FolderApi, SectionName][] = [
    [dom, "domain"], [numerics, "numerics"], [rheo, "viscosity"],
    [initial, "initial condition"], [view, "view"], [trace, "tracers"],
  ];
  for (const [folder, name] of sections) {
    sectionHelpButton(folder, () => hooks.onSectionHelp(name));
  }

  return {
    pane,
    view3d,
    // Assembled here, at the end, rather than accumulated as each control is
    // built: the compiler checks this object against `PaneTargetName` in one
    // place, so a name added to `TOUR_TARGETS` (`tours.ts`) fails the build
    // right here until the blade behind it is actually wired, rather than
    // failing in front of a reader half-way through a tour.
    targets: {
      planet: planetSelect.element,
      preset: preset.element,
      vigour: vigour.element,
      seed: seed.element,
      rock: rock.element,
      paused: paused.element,
      speed: speed.element,
      flow: flow.element,
      tracers: tracers.element,
      colormap: cmap.element,
      restart: restart.element,
      resetView: resetView.element,
      view3d: view3d.element,
      advanced: advanced.element,
      geometry: geom.element,
      boxWidth: len.element,
      walls: walls.element,
      radialWalls: radialWalls.element,
      resolution: resolution.element,
      isothermal: iso.element,
      courant: courant.element,
      dtMax: dtMax.element,
      dtInitial: dtInitial.element,
      law: law.element,
      equation: eq.el,
      contrast: contrast.element,
      depthContrast: depth.element,
      powerLawN: nExp.element,
      cgIterations: iters.element,
      picardSweeps: picard.element,
      yieldStress: sigmaY.element,
      yieldGradient: sigmaB.element,
      etaStar: etaStar.element,
      etaLight: etaLight.element,
      etaDense: etaDense.element,
      seedMode: seedMode.element,
      streamlineDensity: density.element,
      meshOverlay: mesh.element,
      lineWidth: lineWidth.element,
      plotWindow: plotWindow.element,
      tracerOverlay: mode.element,
      tracerCount: count.element,
      tracerColour: tint.element,
      tracerSize: size.element,
      tracerOpacity: opacity.element,
      composition: species.element,
      logRb: rb.element,
      layerDepth: layer.element,
      reseedTracers: reseedTracers.element,
    },
    applyPatch,
    setPlanetTraveling: (traveling) => {
      pane.element.classList.toggle("planet-traveling", traveling);
    },
    selectPlanet,
    set: {
      // Not `applyPatch({ logRa: v })`: that refreshes the whole pane, and a
      // ramp calls this on every frame of a two-second drag. The three lines
      // it does keep are the ones `vigour`'s own `"change"` handler runs.
      logRa: (v) => {
        state.logRa = v;
        if (!state.isothermal) hooks.onRa(10 ** v);
        vigour.refresh();
      },
      // The refresh fires the binding's "change", which repaints the number
      // and moves the log slider — see `logSlider`.
      courant: (v) => {
        state.courant = v;
        courant.refresh();
      },
    },
  };
}
