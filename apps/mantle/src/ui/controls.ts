/**
 * Tweakpane controls.
 *
 * Two layers of grouping sit on top of each other here. The one a reader
 * meets first is *audience*: plain-language controls live at the pane's
 * root, visible always — try an example, convective vigour, how the rock
 * behaves, playback, show flow lines, show tracers (and, once that is
 * checked, colour tracers by), temperature colour map, restart simulation,
 * reset view —
 * and everything else (still every field the solver reads; nothing below is
 * removed) sits in folders hidden behind the "advanced controls" toggle,
 * named for what they let a reader who already knows the physics reach. The
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

import { Pane, type FolderApi } from "tweakpane";
import { COLORMAPS, type ColormapName } from "../colormaps";
import { boundaryNames } from "../geometry";
import {
  PARTICLE_TINT, SIMPLE_PARTICLE_TINT, SPECIES_CONDITIONS, type TintMode,
} from "../particles";
import { colorbarBlock } from "./colorbar";
import { EQUATION, parseFormula } from "./equation";
import { applyOptgroups, deriveGroups } from "./preset-optgroups";
import {
  BENCHMARKS, BOX_LENGTH, CONTRAST, DEPTH_CONTRAST, ETA_VAN_KEKEN, GEOMETRY,
  LABELS, LAYER_DEPTH, LOG_RB, MESH, NU_WINDOWS, PARTICLE_COUNTS,
  PARTICLE_OPACITY, PARTICLE_SIZE, PARTICLES, PRESETS, QUICK_STARTS,
  RADIAL_WALLS, SIMPLE_VISCOSITY, SPEEDS, VISCOSITY, WALLS, type BenchmarkName,
  type GeometryName, type MeshName, type ParticlesName, type PresetName,
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
  /** A benchmark case has just written its fields onto `state`; rebuild from it. */
  onBenchmark(): void;
  onRa(v: number): void;
  /** `Ra` forced to 0 regardless of the slider, or released back to it — a pure uniform write either way. */
  onIsothermal(v: boolean): void;
  onStreamlines(levels: number, lineW: number): void;
  onMesh(m: MeshName): void;
  onColormap(v: ColormapName): void;
  onNuWindow(steps: number): void;
  onReseed(): void;
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
  /** Draw a fresh cloud at the current settings, without touching T — see the plain "reseed" button for the T + particles combination. */
  onReseedParticles(): void;
}

/** Tweakpane list options want `{ label: value }`. */
const nameOptions = <T extends object>(o: T) =>
  Object.fromEntries(Object.keys(o).map((k) => [k, k]));

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

export function buildPane(state: State, hooks: Hooks): Pane {
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
  const preset = pane.addBinding(presetState, "preset", {
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
  preset.on("change", (e) => {
    const name = e.value;
    if (name === CUSTOM) return;
    Object.assign(state, presetTable[name]);
    // Geometry, viscosity and their dependent visibility all just moved, so
    // the same housekeeping their own change handlers do below has to run
    // here too — those handlers fire from pointer/list events on the pane,
    // none of which this write goes through. `enable`/`eq`/`enableBox`/
    // `enableRa`/`enableParticles`/`pcbar` are all defined further down this
    // function; this callback only ever runs after the whole pane (and so
    // every one of those consts) exists.
    enableBox(state.geometry);
    enableRa(state.isothermal);
    enable(state.viscosity);
    enableParticles(state.particles);
    pcbar.setColormap(state.particleColormap);
    eq.redraw();
    syncSimpleControls();
    presetState.preset = CUSTOM;
    pane.refresh();
    hooks.onBenchmark();
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
  const vigour = pane.addBinding(state, "logRa", {
    min: 0, max: 7, step: 0.05, label: "convective vigour",
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

  // ---- how the rock behaves ----
  //
  // Three of `VISCOSITY`'s seven laws, under `SIMPLE_VISCOSITY`'s plain
  // names. Bound to its own object rather than `state.viscosity` directly —
  // the same reason `presetState` is its own object rather than a `State`
  // field: the four laws this list does not offer (Tackley, Tosi,
  // Blankenbach, van Keken) are still reachable under "law" in the advanced
  // viscosity folder, and a plain binding on `state.viscosity` would have
  // nothing sane to display here the moment one of those is picked there.
  // `applyViscosity` is the one place either list's change lands, so the two
  // can never disagree about what selecting a law costs.
  const isSimpleLaw = (v: ViscosityName): boolean =>
    (Object.values(SIMPLE_VISCOSITY) as ViscosityName[]).includes(v);
  const simpleLaw: { law: ViscosityName } =
    { law: isSimpleLaw(state.viscosity) ? state.viscosity : "constant" };
  const rock = pane.addBinding(simpleLaw, "law", {
    options: SIMPLE_VISCOSITY, label: "how the rock behaves",
  });
  // `enable` and `eq` are defined in the advanced viscosity folder below;
  // referenced here only inside a callback, which never runs before the
  // whole pane (and so both consts) exists.
  const applyViscosity = (v: ViscosityName): void => {
    state.viscosity = v;
    enable(v);
    eq.redraw();
    if (isSimpleLaw(v)) simpleLaw.law = v;
    hooks.onViscosity(v);
    pane.refresh();
  };
  rock.on("change", (e) => applyViscosity(e.value as ViscosityName));

  // ---- playback ----
  pane.addBinding(state, "paused");
  // A list rather than a slider: the useful settings span 1/16 to 16 steps per
  // frame, and the labels say what happens far better than a number would.
  pane.addBinding(state, "speed", { options: SPEEDS });

  // ---- show flow lines ----
  //
  // Stands in for the density slider (`contours`, 0–60) and the mesh/line-
  // width pair beneath it in the advanced view folder: a first visit needs
  // to know streamlines exist, not how many. `24` is an arbitrary but
  // reasonable mid-ladder density — the exact count is exactly what the
  // advanced slider is for.
  const SIMPLE_STREAMLINE_DENSITY = 24;
  const simpleFlow = { on: state.contours > 0 };
  // Preserves a nonzero `contours` a benchmark set rather than snapping it to
  // `SIMPLE_STREAMLINE_DENSITY` — same reason `simpleParticles`'s own "on"
  // handler now guards `state.particles`: `pane.refresh()` fires this "change"
  // event on any programmatic flip of `on`, not only a real click.
  pane.addBinding(simpleFlow, "on", { label: "show flow lines" }).on("change", (e) => {
    state.contours = e.value ? (state.contours > 0 ? state.contours : SIMPLE_STREAMLINE_DENSITY) : 0;
    hooks.onStreamlines(state.contours, state.lineWidth);
    pane.refresh();
  });

  // ---- show tracers ----
  //
  // Off ↔ "visual" — the picture worth a first look. "chemical" (the
  // buoyancy-coupled mode) stays reachable only from the full three-way list
  // in the advanced particles folder: turning tracers off here always lands
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
  pane.addBinding(simpleParticles, "on", { label: "show tracers" }).on("change", (e) => {
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
  // same condition; see `enableParticles`. Bound to its own proxy rather than
  // `state.particleTint` directly, the same reason "how the rock behaves" is:
  // `SIMPLE_PARTICLE_TINT` offers two of the full list's seven rows, and a
  // plain binding on `state.particleTint` would have nothing sane to show
  // the moment one of the other five is picked in the advanced folder.
  // `applyTint` is the one place either list's change lands, so the two can
  // never disagree about what colouring a tracer by X means.
  const isSimpleTint = (t: TintMode): boolean =>
    (Object.values(SIMPLE_PARTICLE_TINT) as TintMode[]).includes(t);
  const simpleTintState: { tint: TintMode } =
    { tint: isSimpleTint(state.particleTint) ? state.particleTint : "initial depth" };
  const simpleTint = pane.addBinding(simpleTintState, "tint", {
    options: SIMPLE_PARTICLE_TINT, label: "colour tracers by",
  });
  simpleTint.hidden = !simpleParticles.on;
  const applyTint = (t: TintMode): void => {
    state.particleTint = t;
    state.particleColormap = PARTICLE_TINT[t].colormap;
    pcbar.setColormap(state.particleColormap);
    if (isSimpleTint(t)) simpleTintState.tint = t;
    hooks.onParticleTint();
    pane.refresh();
  };
  simpleTint.on("change", (e) => applyTint(e.value as TintMode));

  // ---- colour map ----
  //
  // Already plain — a swatch, not a formula — so it stays at the root next
  // to the legend it always sat beside.
  const cbar = colorbarBlock(state.colormap, ["0 (cold)", "1 (hot)"]);
  const cmap = pane.addBinding(state, "colormap",
    { options: nameOptions(COLORMAPS), label: "temperature colour map" });
  cmap.on("change", (e) => {
    cbar.setColormap(e.value as ColormapName);
    hooks.onColormap(e.value as ColormapName);
  });

  // ---- restart simulation ----
  //
  // Re-seeds T and re-solves Stokes from it, and — if a tracer cloud is
  // attached — redraws that too, so every field the picture shows starts
  // over together rather than restarting T and leaving a stale cloud
  // behind. What it restarts *into* (seed mode, composition, species) is
  // still set from the advanced seeding and particles folders; this is the
  // one-click "start over with what's already set" a first-time reader
  // reaches for without touching either.
  pane.addButton({ title: "restart simulation" }).on("click", () => {
    hooks.onReseed();
    if (PARTICLES[state.particles].attached) hooks.onReseedParticles();
  });

  // Scroll to zoom, drag to pan (see main.ts) — this is the way back from
  // either with no pointer precision required.
  pane.addButton({ title: "reset view" }).on("click", () => hooks.onResetView());

  // ---- advanced controls ----
  //
  // Built here, right after the last plain-language control and before any
  // advanced folder — not after them — so this stays put in the rack
  // regardless of whether it is checked. Placed *after* the folders (as
  // "built last" once was), checking it un-hides several screens of content
  // that were sitting between this checkbox and "reset view", which shoves
  // the checkbox itself far down the pane the instant it is clicked — the
  // opposite of what a fixed anchor is for. Here, every advanced folder
  // renders *below* this line whether hidden or not, so this is always the
  // last thing directly under "reset view".
  const ui = { advanced: false };
  pane.addBinding(ui, "advanced", { label: "advanced controls" })
    .on("change", (e) => {
      for (const f of advancedFolders) f.hidden = !e.value;
      // The convective-vigour/log₁₀-Ra control's other face — see its own
      // note above on why this is one binding rather than two.
      vigour.label = e.value ? "log₁₀ Ra" : "convective vigour";
      if (vigourNumber) vigourNumber.style.display = e.value ? "" : "none";
    });

  /**
   * Re-reads `state` into the four simple proxies above, for whichever of
   * them a change made elsewhere (a preset, or an advanced control) may have
   * moved without going through the simple control itself. Cheap and always
   * safe to call — each line is a no-op unless the two actually disagree.
   */
  const syncSimpleControls = (): void => {
    if (isSimpleLaw(state.viscosity)) simpleLaw.law = state.viscosity;
    simpleFlow.on = state.contours > 0;
    simpleParticles.on = state.particles !== "off";
    simpleTint.hidden = !simpleParticles.on;
    if (isSimpleTint(state.particleTint)) simpleTintState.tint = state.particleTint;
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
  const walls = dom.addBinding(state, "walls",
    { options: nameOptions(WALLS), label: "left / right" });
  // Unlike `walls`, legal on *both* geometries — a no-slip radial condition
  // means the same thing on an annulus (inner/outer) as on a box (top/
  // bottom), so this is never disabled the way `len`/`walls` are, only
  // relabelled to say which pair of boundaries it closes. See
  // `boundaryNames` in `geometry.ts`, the same names the Nusselt readout uses.
  const radialWalls = dom.addBinding(state, "radialWalls",
    { options: nameOptions(RADIAL_WALLS), label: "top / bottom" });
  const enableBox = (g: GeometryName) => {
    len.disabled = walls.disabled = GEOMETRY[g] !== "box";
    const bn = boundaryNames(GEOMETRY[g]);
    radialWalls.label = `${bn.inner} / ${bn.outer}`;
  };
  geom.on("change", (e) => {
    enableBox(e.value as GeometryName);
    hooks.onGeometry();
  });
  // On release only. The width changes the azimuthal knot vector, so it is the
  // same second-or-two rebuild the resolution list is; firing it per pointer
  // move would queue one for every pixel dragged. The list has no drag to wait
  // for, so it fires on change like every other list in the pane.
  len.on("change", (e) => { if (e.last) hooks.onGeometry(); });
  walls.on("change", () => hooks.onGeometry());
  radialWalls.on("change", () => hooks.onGeometry());
  enableBox(state.geometry);
  dom.addBinding(state, "resolution", { options: nameOptions(PRESETS) })
    .on("change", (e) => hooks.onResolution(e.value as PresetName));

  // What actually drives the step, plus the ceiling it is held under, and the
  // isothermal override — everything about the solve that isn't "how vigorous"
  // or "which law", both of which moved to the simple root above.
  const numerics = pane.addFolder({ title: "numerics" });
  advancedFolders.push(numerics);
  // Forces Ra = 0 regardless of the convection-vigour slider above — the
  // purely compositional (isothermal) buoyancy the van Keken Rayleigh–Taylor
  // benchmark needs (see `isothermal`'s own header in presets.ts on why this
  // is a checkbox rather than a widened `logRa` floor). Hides the vigour
  // slider rather than disabling it: while this is checked, `logRa`'s value
  // is not what is being solved with, and a slider that is still draggable
  // but silently ignored is worse than one that is briefly not there.
  const iso = numerics.addBinding(state, "isothermal", { label: "isothermal (Ra = 0)" });
  const enableRa = (isothermal: boolean): void => { vigour.hidden = isothermal; };
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
  // a couple of pixels wide. There's no `log` option on a plain binding, so
  // the built-in slider strip (`.tp-sldv`, inside the `.tp-sldtxtv_s` half of
  // this composite view) is hidden and a native `<input type="range">` —
  // dragged in log₁₀ space, from -1 to 2 — stands in for it. The number field
  // stays Tweakpane's own and keeps reading/accepting real Courant values;
  // only the drag mapping is replaced, via the same "reach into our own DOM"
  // move as the colour above, so a Tweakpane upgrade that renames these
  // classes loses the slider, not the control.
  const courantSliderWrap =
    courant.element.querySelector<HTMLElement>(".tp-sldtxtv_s");
  const courantNativeSlider =
    courantSliderWrap?.querySelector<HTMLElement>(".tp-sldv");
  if (courantNativeSlider) courantNativeSlider.style.display = "none";
  const courantLogSlider = document.createElement("input");
  courantLogSlider.type = "range";
  courantLogSlider.className = "co-log-slider";
  courantLogSlider.min = "-1";
  courantLogSlider.max = "2";
  courantLogSlider.step = "0.001";
  courantLogSlider.value = String(Math.log10(state.courant));
  courantSliderWrap?.appendChild(courantLogSlider);
  courantLogSlider.addEventListener("input", () => {
    state.courant =
      Math.min(100, Math.max(0.1, 10 ** courantLogSlider.valueAsNumber));
    pane.refresh();
    paintCourant(state.courant);
  });
  // Fires on every drag tick and on committing a typed value alike (see the
  // Ra binding above), so typing a number directly into the field also drags
  // the log slider's handle to match.
  courant.on("change", (e) => {
    paintCourant(e.value);
    courantLogSlider.value = String(Math.log10(e.value));
  });
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
  // Same treatment as the Courant slider immediately above, and for the same
  // reason: no `log` option on a plain binding, so Tweakpane's own linear
  // strip is hidden and a native `<input type="range">` dragged in log₁₀
  // space (−4 to 3) stands in for it, sharing that slider's `co-log-slider`
  // styling rather than a second copy of it.
  const dtMaxSliderWrap = dtMax.element.querySelector<HTMLElement>(".tp-sldtxtv_s");
  const dtMaxNativeSlider = dtMaxSliderWrap?.querySelector<HTMLElement>(".tp-sldv");
  if (dtMaxNativeSlider) dtMaxNativeSlider.style.display = "none";
  const dtMaxLogSlider = document.createElement("input");
  dtMaxLogSlider.type = "range";
  dtMaxLogSlider.className = "co-log-slider";
  dtMaxLogSlider.min = "-4";
  dtMaxLogSlider.max = "3";
  dtMaxLogSlider.step = "0.001";
  dtMaxLogSlider.value = String(Math.log10(state.dtMax));
  dtMaxSliderWrap?.appendChild(dtMaxLogSlider);
  dtMaxLogSlider.addEventListener("input", () => {
    state.dtMax = Math.min(1e3, Math.max(1e-4, 10 ** dtMaxLogSlider.valueAsNumber));
    pane.refresh();
  });
  // Fires on every drag tick and on committing a typed value alike (see the
  // Ra binding's own note on this), so typing a number directly into the
  // field also drags the log slider's handle to match.
  dtMax.on("change", (e) => { dtMaxLogSlider.value = String(Math.log10(e.value)); });
  // 1e-6 to 1e3: the step `GpuSimulation.create` is seeded with, before the
  // first `pollStats` readback gives `adaptiveDt` a CFL-implied value to work
  // from (see `dtInitial` on `State`). Only read at build time — unlike the
  // cap immediately above, changing it has no effect on a solver already
  // running, only the next one built.
  const dtInitial = numerics.addBinding(state, "dtInitial", {
    min: 1e-6, max: 1e3, step: 1e-6, label: "dt initial",
    format: (v) => v.toExponential(1),
  });
  // Same "hide Tweakpane's own linear strip, drag a log-space native slider
  // in its place" treatment as the cap immediately above.
  const dtInitialSliderWrap =
    dtInitial.element.querySelector<HTMLElement>(".tp-sldtxtv_s");
  const dtInitialNativeSlider =
    dtInitialSliderWrap?.querySelector<HTMLElement>(".tp-sldv");
  if (dtInitialNativeSlider) dtInitialNativeSlider.style.display = "none";
  const dtInitialLogSlider = document.createElement("input");
  dtInitialLogSlider.type = "range";
  dtInitialLogSlider.className = "co-log-slider";
  dtInitialLogSlider.min = "-6";
  dtInitialLogSlider.max = "3";
  dtInitialLogSlider.step = "0.001";
  dtInitialLogSlider.value = String(Math.log10(state.dtInitial));
  dtInitialSliderWrap?.appendChild(dtInitialLogSlider);
  dtInitialLogSlider.addEventListener("input", () => {
    state.dtInitial =
      Math.min(1e3, Math.max(1e-6, 10 ** dtInitialLogSlider.valueAsNumber));
    pane.refresh();
  });
  dtInitial.on("change", (e) => {
    dtInitialLogSlider.value = String(Math.log10(e.value));
  });

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

  // Initial condition: which seed mode a fresh run starts from, and the
  // button that redraws one. Split out from the numerics folder above
  // because these are decisions about the *starting picture*, not about how
  // accurately the solve tracks it once running.
  const seeding = pane.addFolder({ title: "seeding" });
  advancedFolders.push(seeding);
  seeding.addBinding(state, "wavenumber", { min: 1, max: 12, step: 1, label: "seed mode" });
  seeding.addButton({ title: "reseed" }).on("click", () => hooks.onReseed());

  // The streamline density, the mesh overlay and the line width both draw
  // with — colour map lives at the simple root now, next to the legend it
  // has always sat beside, so it is not repeated here.
  const view = pane.addFolder({ title: "view detail" });
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
  view.addBinding(state, "mesh", { options: nameOptions(MESH) })
    .on("change", (e) => hooks.onMesh(e.value as MeshName));
  view.addBinding(state, "lineWidth", { min: 0.5, max: 3, step: 0.1, label: "line width" })
    .on("change", (e) => hooks.onStreamlines(state.contours, e.value));
  // How much of the run the two corner plots show — Nusselt number and RMS
  // velocity share this one control (see `presets.ts`). Costs nothing: both
  // traces keep every sample either way, so this re-scales an existing
  // buffer and does not begin collecting again.
  view.addBinding(state, "nuWindow", { options: NU_WINDOWS, label: "Nu window" })
    .on("change", (e) => hooks.onNuWindow(e.value));

  // The tracer overlay in full: the three-way mode list (the simple "show
  // tracers" switch only ever reaches "off" and "visual" — "chemical" lives
  // here), plus every control that only means something once a cloud
  // exists. Structured like the viscosity folder above: one list decides
  // which controls beneath it mean anything, and those are hidden rather
  // than disabled.
  const trace = pane.addFolder({ title: "particles detail" });
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
    label: "dot size",
  });
  const opacity = trace.addBinding(state, "particleOpacity", {
    min: PARTICLE_OPACITY.min, max: PARTICLE_OPACITY.max, step: PARTICLE_OPACITY.step,
    label: "dot opacity",
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
  trace.addButton({ title: "reseed particles" }).on("click", () => hooks.onReseedParticles());

  const enableParticles = (m: ParticlesName): void => {
    const { attached, coupled } = PARTICLES[m];
    count.hidden = tint.hidden = size.hidden = opacity.hidden = !attached;
    pcbar.el.hidden = !attached;
    species.hidden = rb.hidden = layer.hidden = !coupled;
    // The simple root's own "colour tracers by" reads the same condition as
    // this folder's `tint` — one function deciding it for both, so the two
    // can never disagree about when a cloud exists to colour.
    simpleTint.hidden = !attached;
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
  // "advanced controls" binding above flips all of them together.
  for (const f of advancedFolders) f.hidden = true;

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

  return pane;
}
