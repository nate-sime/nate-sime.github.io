/**
 * Tweakpane controls.
 *
 * The controls are grouped by what they *cost*, because that is the honest
 * distinction here and it is invisible from the labels:
 *
 *   Ra, contours, line width,   — a 128-byte uniform write; next frame.
 *   mesh, n
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
 *                                 the one profile, so either fires the one job.
 *   viscosity tier, resolution, — rebuilds every table and pipeline, 1.3–2.7 s,
 *   geometry, box length          and the page says so rather than appearing to
 *                                 hang. Only *entering or leaving* the Krylov
 *                                 tier does this; μ(T, d) ↔ μ(T, d, ε̇) is a
 *                                 uniform.
 *                                 Geometry is a rebuild because the metric is
 *                                 compiled into the shaders and the box length
 *                                 reaches the knot vector — see `presets.ts`.
 *
 * This module owns no simulation state: it mutates `State` and calls back. That
 * keeps the pane replaceable (and absent, in tests) without the solver noticing.
 */

import { Pane } from "tweakpane";
import { COLORMAPS, type ColormapName } from "../colormaps";
import { colorbarBlock } from "./colorbar";
import { EQUATION, parseFormula } from "./equation";
import {
  BOX_LENGTH, GEOMETRY, LABELS, MESH, NU_WINDOWS, PRESETS, SPEEDS, VISCOSITY,
  WALLS, type GeometryName, type MeshName, type PresetName, type State,
  type ViscosityName, type WallsName,
} from "./presets";

export type {
  GeometryName, MeshName, PresetName, State, ViscosityName, WallsName,
};
export {
  GEOMETRY, MESH, NU_WINDOWS, PRESETS, SPEEDS, VISCOSITY, WALLS,
  defaultState, geometryFor,
} from "./presets";

export interface Hooks {
  onRa(v: number): void;
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
  onResetView(): void;
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

  // *What* is being solved, above everything about how. Both controls in here
  // rebuild every table and pipeline (see `presets.ts`), so both are announced;
  // the length is disabled rather than hidden on the annulus, so selecting a
  // geometry does not move the rest of the pane out from under the pointer.
  const dom = pane.addFolder({ title: "domain" });
  const geom = dom.addBinding(state, "geometry",
    { options: nameOptions(GEOMETRY), label: "geometry" });
  const len = dom.addBinding(state, "boxLength", {
    min: BOX_LENGTH.min, max: BOX_LENGTH.max, step: BOX_LENGTH.step,
    label: "box width",
  });
  // Below the width, because it is a statement about the domain's *edges* and
  // reads as one only once there is a width for them to be the edges of.
  const walls = dom.addBinding(state, "walls",
    { options: nameOptions(WALLS), label: "left / right" });
  const enableBox = (g: GeometryName) => {
    len.disabled = walls.disabled = GEOMETRY[g] !== "box";
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
  enableBox(state.geometry);

  const flow = pane.addFolder({ title: "flow" });
  // Ra spans decades and the interesting behaviour (onset, then plume count) is
  // logarithmic in it, so a linear slider would waste most of its travel.
  flow.addBinding(state, "logRa", { min: 3, max: 7, step: 0.05, label: "log₁₀ Ra" })
    .on("change", (e) => hooks.onRa(10 ** e.value));
  // What actually drives the step, plus the ceiling it is held under. Neither
  // is a factorisation itself — `main.ts` reads both every poll and calls
  // `setDt` only when the CFL-implied step has moved past `adaptiveDt`'s
  // hysteresis band — so, unlike the old single `dt` slider, both take effect
  // while dragging rather than needing a release guard.
  flow.addBinding(state, "courant", { min: 0.1, max: 3, step: 0.1, label: "Courant number" });
  flow.addBinding(state, "dtMax", { min: 2e-5, max: 5e-4, step: 1e-5, label: "dt cap" });

  // Viscosity: the law list picks the rheology, and the knobs below it only mean
  // anything for some of them — so they are disabled rather than hidden, which
  // would make the pane jump and lose the reader's place. Two levels of that:
  // contrast and the CG budget need the Krylov tier, n and the Picard sweeps
  // need the power law on top of it.
  const rheo = pane.addFolder({ title: "viscosity" });
  const law = rheo.addBinding(state, "viscosity",
    { options: nameOptions(VISCOSITY), label: "law" });
  const eq = equationBlock(state);
  const contrast = rheo.addBinding(state, "logContrast",
    { min: 0, max: 5, step: 0.25, label: LABELS.contrast });
  // Directly below the thermal contrast, because they are the same kind of
  // number — a log₁₀ ratio across the layer — and reading them as a pair is what
  // says the total contrast is their product. Its floor is 0 (no depth
  // dependence, the law the app opens with) and its ceiling is lower than the
  // thermal one's: the two multiply inside one clamp, and 10⁵ of each is a
  // contrast no fixed Krylov budget is going to hold.
  const depth = rheo.addBinding(state, "logDepthContrast",
    { min: 0, max: 3, step: 0.25, label: LABELS.depth });
  const nExp = rheo.addBinding(state, "n",
    { min: 1, max: 5, step: 0.25, label: LABELS.n });
  const iters = rheo.addBinding(state, "iters",
    { min: 1, max: 40, step: 1, label: "CG iterations" });
  const picard = rheo.addBinding(state, "picard",
    { min: 1, max: 3, step: 1, label: "Picard sweeps" });

  const enable = (v: ViscosityName) => {
    const { variable, strainRate } = VISCOSITY[v];
    contrast.disabled = depth.disabled = iters.disabled = !variable;
    nExp.disabled = picard.disabled = !strainRate;
  };
  law.on("change", (e) => {
    enable(e.value as ViscosityName);
    eq.redraw();
    hooks.onViscosity(e.value as ViscosityName);
  });
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
  enable(state.viscosity);

  const run = pane.addFolder({ title: "run" });
  run.addBinding(state, "paused");
  // A list rather than a slider: the useful settings span 1/16 to 16 steps per
  // frame, and the labels say what happens far better than a number would.
  run.addBinding(state, "speed", { options: SPEEDS });
  run.addBinding(state, "wavenumber", { min: 1, max: 12, step: 1, label: "seed mode" });
  run.addButton({ title: "reseed" }).on("click", () => hooks.onReseed());

  // Both overlays start off and both are one uniform write, so nothing here is
  // announced. The width serves both — it is the only line weight in the render
  // pass — which is why it sits below the two things it applies to rather than
  // under the streamlines alone.
  const view = pane.addFolder({ title: "view" });
  // The map and its legend come first in the folder — what the colour means,
  // ahead of what is drawn over it.
  const cbar = colorbarBlock(state.colormap);
  const cmap = view.addBinding(state, "colormap",
    { options: nameOptions(COLORMAPS), label: "colour map" });
  cmap.on("change", (e) => {
    cbar.setColormap(e.value as ColormapName);
    hooks.onColormap(e.value as ColormapName);
  });
  view.addBinding(state, "contours", { min: 0, max: 60, step: 2, label: "streamlines" })
    .on("change", (e) => hooks.onStreamlines(e.value, state.lineWidth));
  view.addBinding(state, "mesh", { options: nameOptions(MESH) })
    .on("change", (e) => hooks.onMesh(e.value as MeshName));
  view.addBinding(state, "lineWidth", { min: 0.5, max: 3, step: 0.1, label: "line width" })
    .on("change", (e) => hooks.onStreamlines(state.contours, e.value));
  // How much of the run the two corner plots show — Nusselt number and RMS
  // velocity share this one control (see `presets.ts`). It belongs beside the
  // overlays because it is the same kind of control — what is drawn, not what
  // is solved — and it costs the same nothing: both traces keep every sample
  // either way, so this re-scales an existing buffer and does not begin
  // collecting again.
  view.addBinding(state, "nuWindow", { options: NU_WINDOWS, label: "Nu window" })
    .on("change", (e) => hooks.onNuWindow(e.value));
  // Scroll to zoom, drag to pan (see main.ts) — this is the way back from
  // either with no pointer precision required.
  view.addButton({ title: "reset view" }).on("click", () => hooks.onResetView());

  pane.addBinding(state, "resolution", { options: nameOptions(PRESETS) })
    .on("change", (e) => hooks.onResolution(e.value as PresetName));

  // The equation goes between the law and the knobs, because that is what it
  // connects: the law the list just selected, and the sliders below named
  // against the symbols they set. It is inserted *last* — a rack re-appends its
  // blades' elements as they are added, so a foreign node placed mid-folder
  // drifts to the bottom of it as the rest of the folder is built. The colour
  // bar is the same trick for the same reason: it belongs right after the
  // colour-map list, and `view` gains three more blades after that one.
  law.element.after(eq.el);
  cmap.element.after(cbar.el);

  return pane;
}
