/**
 * The guided tour: the "Tutorial" chip in `#hud`, and the overlay it opens.
 *
 * The screen dims, one control at a time is left lit, and a card beside it
 * says what that control does and what it means for the physics. The tour
 * drives the model as it goes — loading a scene, dragging the vigour slider,
 * flying the camera onto a boundary layer — so a reader watching it happen
 * has seen the connection between the control and the picture before being
 * asked to believe in it. The steps themselves are data (`tours.ts`); this
 * file is only the machinery that makes one true.
 *
 * **Nothing ever dims the field.** The whole point of a step is that the
 * reader watches the canvas while a control is explained, so the simulation
 * is at full brightness for every step of the tour. What goes down is the
 * *chrome* — the four containers `.chrome-hidden` names in index.html — and
 * the one holding the lit control is exempted from that and gets the
 * spotlight treatment instead.
 *
 * **That spotlight is four rectangles, not one element with a hole cut in
 * it.** A `clip-path` would be less markup, but the lit control has to stay
 * both visible *and* clickable at its own stacking level, and the pane it
 * lives in (`#pane`) declares no `z-index` at all — raising it above a
 * covering overlay to un-dim one blade would un-dim the whole panel. Four
 * rectangles tiled around the target simply never cover it: nothing is on top
 * to see through or click past, and each one animates its own
 * `top/left/width/height` rather than asking two `path()`s to interpolate.
 * They are tiled between the lit control and its own container, not the
 * viewport, which is what keeps them off the canvas; they take pointer
 * events, which is what leaves only that one control reachable inside it.
 *
 * **Nothing is added to the frame loop.** `main.ts`'s loop is the solver's;
 * this file runs its own `requestAnimationFrame` while the overlay is open,
 * which re-measures the target, advances a slider ramp, and fills the dwell
 * bar. Re-measuring every frame rather than subscribing to anything is what
 * makes the spotlight track a `--ui-scale` drag, a window resize and a scroll
 * inside `#pane` without this file knowing any of those exist.
 *
 * **A dwell does not advance the step.** It fills a bar and then puts a
 * highlight on "next" — a nudge, not a gate. The tour is meant to be walked
 * with two buttons, so a step that moved on by itself would be taking the one
 * decision the reader was left, usually part-way through the sentence
 * explaining what they were supposed to be watching for.
 */

import { BENCHMARKS, QUICK_STARTS, type State } from "./presets";
import {
  DEFAULT_TOUR, TOURS, type TourDwell, type TourStep, type TourTargetName,
} from "./tours";

/**
 * What the tour needs from the rest of the app. Everything here already
 * exists somewhere — the pane's own `applyPatch` and `set` (`controls.ts`),
 * the 2-D camera and `Globe3D` (`main.ts`) — and is passed in rather than
 * imported so this file holds no reference to the solver, the pane or the
 * device, and can be built before any of the three has settled.
 */
export interface TourActions {
  /** The element to light for a step's `target`; `null` if it isn't on screen. */
  element(name: TourTargetName): HTMLElement | null;
  /** `PaneHandle.applyPatch` — the path a preset selection takes. */
  applyPatch(patch: Partial<State>): void;
  /** `PaneHandle.set.logRa`, called every frame of a ramp. */
  setLogRa(v: number): void;
  /** The live `State`, read for a ramp's starting point and for the snapshot. */
  readState(): Readonly<State>;
  /** `Globe3D.viewMode`, so a step's `view` is only acted on when the two differ. */
  viewMode(): "2d" | "3d";
  toggle3D(): void;
  /** Fly the flat view's camera; world coordinates, as `TourFocus` states them. */
  focusOn(zoom: number, x: number, y: number, ms: number): void;
  resetFocus(): void;
  /** Steps the solver has taken, which is what a `{ steps }` dwell counts. */
  steps(): number;
}

const PRESETS_BY_NAME: Record<string, Partial<State>> = { ...QUICK_STARTS, ...BENCHMARKS };

/** Breathing room between the lit control and the hole's edge. */
const PAD = 6;
/** Between the hole and the card beside it. */
const GAP = 14;
/** Between the card and the edge of the window. */
const MARGIN = 12;
/** Below this the card goes bottom-centre regardless — the same width index.html's one media query turns on. */
const NARROW = 600;

interface Box { x: number; y: number; w: number; h: number }

/**
 * The chrome containers a lit control can live in — the same four
 * `.chrome-hidden` names (index.html), which between them hold every target
 * in `TOUR_TARGETS` that is not the canvas. The spotlight is tiled inside
 * whichever one of these contains the step's target, and that one alone is
 * exempted from the wholesale dimming the other three get.
 */
const HOSTS = "#hud, #traces, #scale, #pane";

const same = (a: Box, b: Box): boolean =>
  a.x === b.x && a.y === b.y && a.w === b.w && a.h === b.h;

const px = (v: number): string => `${Math.round(v)}px`;

const setBox = (el: HTMLElement, b: Box): void => {
  el.style.left = px(b.x);
  el.style.top = px(b.y);
  el.style.width = px(Math.max(0, b.w));
  el.style.height = px(Math.max(0, b.h));
};

/** Smoothstep, the same easing curve `Globe3D`'s own transition uses. */
const ease = (t: number): number => t * t * (3 - 2 * t);

const div = (cls: string, parent: HTMLElement): HTMLDivElement => {
  const e = document.createElement("div");
  e.className = cls;
  parent.append(e);
  return e;
};

const button = (label: string, parent: HTMLElement): HTMLButtonElement => {
  const b = document.createElement("button");
  b.type = "button";
  b.className = "tour-btn";
  b.textContent = label;
  parent.append(b);
  return b;
};

/**
 * Wires the static `#tour-chip`/`#tour` pair (index.html) together, the same
 * shape `buildAcknowledgements` uses for its own chip and popover — but built
 * from inside `main()` rather than beside it, because unlike that list this
 * needs a live solver and a built pane to have anything to drive.
 */
export function buildTour(chip: HTMLElement, root: HTMLElement, actions: TourActions): void {
  const steps: readonly TourStep[] = TOURS[DEFAULT_TOUR];

  // Motion is the tour's main device — a slider seen to travel, a camera seen
  // to fly — so this is the one place in the app that has to ask. Honoured
  // here for the ramp and the camera; index.html's own media query handles
  // the shades and the ring. Read once rather than watched: a reader changing
  // this system-wide preference mid-tour is not a case worth the listener.
  const calm = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  // ---- the overlay ------------------------------------------------------
  //
  // Built once, at construction, and hidden rather than torn down: the tour
  // is opened and closed repeatedly against the same page, and the shades'
  // transitions need something already painted to animate from (the same
  // sequencing `acknowledgements.ts` documents for its own popover).
  const shades = ["top", "bottom", "left", "right"]
    .map((side) => div(`tour-shade tour-${side}`, root));
  const ring = div("tour-ring", root);
  const card = div("tour-card", root);

  const count = div("tour-count", card);
  const title = document.createElement("h2");
  title.className = "tour-title";
  card.append(title);
  const bodyEl = div("tour-body", card);
  const watch = div("tour-watch", card);
  const bar = div("tour-bar", card);
  const barFill = div("tour-bar-fill", bar);
  // Its own row above the navigation rather than a fourth button in it: at
  // the card's width four buttons wrap, and this one is not a way *through*
  // the tour like the other three — it undoes everything the tour did.
  const restoreRow = div("tour-restore", card);
  const restore = button("restore the run I had", restoreRow);
  const nav = div("tour-nav", card);
  const end = button("end tour", nav);
  end.classList.add("tour-btn-quiet");
  const spacer = div("tour-spacer", nav);
  spacer.setAttribute("aria-hidden", "true");
  const back = button("back", nav);
  const next = button("next", nav);
  next.classList.add("tour-btn-primary");

  card.setAttribute("role", "dialog");
  card.setAttribute("aria-live", "polite");

  // ---- state ------------------------------------------------------------
  let at = -1;                       // -1 while closed
  let raf = 0;
  // Both last written, so the loop below only touches the DOM when one moves.
  let host: Box | null = null;
  let hole: Box | null = null;
  let ramp: { from: number; to: number; t0: number; ms: number } | null = null;
  let dwell: { kind: TourDwell; t0: number; base: number } | null = null;
  let dwellDone = false;
  /**
   * The whole `State` as the reader left it, taken at open. The tour is not
   * put back automatically on the way out — it ends in a configured, running
   * model, and that is the useful place to be left — but the last card offers
   * this, applied through the same `applyPatch` a preset goes through.
   */
  let snapshot: State | null = null;

  const open = (): boolean => at >= 0;

  // ---- geometry ---------------------------------------------------------

  const viewport = (): Box => ({ x: 0, y: 0, w: window.innerWidth, h: window.innerHeight });

  /**
   * The element a step lights, and the chrome container it sits in.
   *
   * A step pointing at the canvas — or at nothing — has no host: the field is
   * never what gets dimmed (see this file's header), so there is no spotlight
   * to tile and the chrome simply all goes down together.
   */
  const targetOf = (step: TourStep): { el: HTMLElement; host: HTMLElement } | null => {
    const el = step.target === null || step.target === "canvas"
      ? null : actions.element(step.target);
    const host = el?.closest<HTMLElement>(HOSTS);
    return el && host ? { el, host } : null;
  };

  /**
   * The hole, in viewport coordinates, clipped to its host.
   * `getBoundingClientRect` is deliberately the only measurement taken: it
   * already has `--ui-scale`'s transform on `#pane`/`#hud`/`#traces` folded
   * in, so the spotlight tracks that slider without this file reading the
   * custom property at all. The clip is what stops a blade scrolled half out
   * of `#pane`'s own scroll box from punching a hole through the panel's edge.
   */
  const holeIn = (el: HTMLElement, host: Box): Box => {
    const r = el.getBoundingClientRect();
    const x0 = Math.max(host.x, r.left - PAD);
    const y0 = Math.max(host.y, r.top - PAD);
    const x1 = Math.min(host.x + host.w, r.right + PAD);
    const y1 = Math.min(host.y + host.h, r.bottom + PAD);
    return { x: x0, y: y0, w: Math.max(0, x1 - x0), h: Math.max(0, y1 - y0) };
  };

  /**
   * Tile the four shades in the gap between `host` and `hole`, and put the
   * ring on the hole.
   *
   * With no target the two are the same rectangle — the whole viewport — and
   * every one of the four comes out zero-area against an edge. That is not a
   * special case bolted on: it is what "the hole is everything" means, and it
   * makes the shades slide off the sides of the screen on the way into a
   * canvas step rather than collapsing toward a corner.
   */
  const paint = (host: Box, hole: Box, outline: boolean): void => {
    const [top, bottom, left, right] = shades;
    setBox(top, { x: host.x, y: host.y, w: host.w, h: hole.y - host.y });
    setBox(bottom, {
      x: host.x, y: hole.y + hole.h,
      w: host.w, h: host.y + host.h - hole.y - hole.h,
    });
    setBox(left, { x: host.x, y: hole.y, w: hole.x - host.x, h: hole.h });
    setBox(right, {
      x: hole.x + hole.w, y: hole.y,
      w: host.x + host.w - hole.x - hole.w, h: hole.h,
    });
    setBox(ring, hole);
    // Not `host !== hole`: a step can point at a whole container (the corner
    // traces are one target, not a control inside one), where the two are the
    // same rectangle and the shades come out empty — but the outline saying
    // "this panel" is exactly as wanted there as it is around a single blade.
    ring.style.opacity = outline ? "1" : "0";
  };

  /**
   * Put the card next to the hole, on whichever side has room for it, then
   * clamp it into the window. Preference order is horizontal first: the
   * controls this tour points at are a vertical strip down the right-hand
   * pane, so beside is nearly always both available and the reading order
   * a reader expects.
   */
  const place = (b: Box, beside: boolean): void => {
    const vw = window.innerWidth, vh = window.innerHeight;
    const { width: cw, height: ch } = card.getBoundingClientRect();
    let x: number, y: number;
    // Nothing to sit beside on a canvas step, so the card goes where it
    // covers least of what the step is talking about. Bottom-centre clears
    // the HUD's corner, the corner traces, and the top of the annulus, which
    // is where the one zooming step is looking.
    if (!beside || vw < NARROW) {
      x = (vw - cw) / 2;
      y = vh - ch - MARGIN;
    } else if (b.x - GAP - cw >= MARGIN) {
      x = b.x - GAP - cw;
      y = b.y;
    } else if (b.x + b.w + GAP + cw <= vw - MARGIN) {
      x = b.x + b.w + GAP;
      y = b.y;
    } else if (b.y + b.h + GAP + ch <= vh - MARGIN) {
      x = b.x + b.w / 2 - cw / 2;
      y = b.y + b.h + GAP;
    } else {
      x = b.x + b.w / 2 - cw / 2;
      y = b.y - GAP - ch;
    }
    card.style.left = px(Math.min(vw - cw - MARGIN, Math.max(MARGIN, x)));
    card.style.top = px(Math.min(vh - ch - MARGIN, Math.max(MARGIN, y)));
  };

  // ---- the loop ---------------------------------------------------------

  const tick = (now: number): void => {
    if (!open()) return;
    const step = steps[at];

    // Both measured every frame and written only when either moves: that is
    // what tracks a `--ui-scale` drag, a window resize and a scroll inside
    // `#pane` without this file subscribing to any of the three.
    const found = targetOf(step);
    const h = found ? found.host.getBoundingClientRect() : null;
    const nextHost: Box = h ? { x: h.left, y: h.top, w: h.width, h: h.height } : viewport();
    const nextHole = found ? holeIn(found.el, nextHost) : nextHost;
    if (!host || !hole || !same(host, nextHost) || !same(hole, nextHole)) {
      host = nextHost;
      hole = nextHole;
      paint(nextHost, nextHole, found !== null);
      place(nextHole, found !== null);
    }

    if (ramp) {
      const t = ramp.ms <= 0 ? 1 : Math.min(1, (now - ramp.t0) / ramp.ms);
      actions.setLogRa(ramp.from + (ramp.to - ramp.from) * ease(t));
      if (t >= 1) ramp = null;
    }

    if (dwell) {
      const p = "steps" in dwell.kind
        ? (actions.steps() - dwell.base) / dwell.kind.steps
        : (now - dwell.t0) / dwell.kind.ms;
      const clamped = Math.min(1, Math.max(0, p));
      barFill.style.width = `${(clamped * 100).toFixed(1)}%`;
      if (clamped >= 1 && !dwellDone) {
        dwellDone = true;
        // The nudge: the bar is full, so "next" stops being just the button
        // that happens to be there and starts asking to be pressed.
        next.classList.add("tour-btn-ready");
      }
    }

    raf = requestAnimationFrame(tick);
  };

  // ---- steps ------------------------------------------------------------

  const render = (step: TourStep): void => {
    count.textContent = `${at + 1} / ${steps.length}`;
    title.textContent = step.title;
    bodyEl.replaceChildren(...step.body.map((t) => {
      const p = document.createElement("p");
      p.textContent = t;
      return p;
    }));
    watch.textContent = step.watch ?? "";
    watch.style.display = step.watch ? "" : "none";
    bar.style.display = step.dwell ? "" : "none";
    barFill.style.width = "0%";
    back.disabled = at === 0;
    const last = at === steps.length - 1;
    next.textContent = last ? "finish" : "next";
    next.classList.toggle("tour-btn-ready", last);
    // Only offered where it makes sense to take it: at the end, next to the
    // sentence that says the run has been left where the tour finished.
    restoreRow.style.display = last && snapshot ? "" : "none";
  };

  const enter = (i: number): void => {
    at = i;
    const step = steps[i];

    // The view first: a `focus` below is the flat camera's, and there is no
    // point flying it while the globe is what is actually being drawn.
    // Reconciled rather than asserted, so stepping backwards doesn't toggle
    // the view on a step that never asked about it.
    if (step.view && actions.viewMode() !== step.view) actions.toggle3D();

    // Preset before patch: a step may load a scene and then correct one field
    // of it (the opening example turns off the isothermal override and
    // un-pauses, neither of which a `QUICK_STARTS` entry states).
    if (step.preset) actions.applyPatch(PRESETS_BY_NAME[step.preset]);
    if (step.patch) actions.applyPatch(step.patch);

    if (step.focus === "reset") actions.resetFocus();
    else if (step.focus) {
      const { zoom, x, y, ms } = step.focus;
      actions.focusOn(zoom, x, y, calm ? 0 : ms);
    }

    // After the patch, so the ramp starts from whatever that left behind
    // rather than from a value it is about to overwrite.
    if (step.ramp) {
      const from = actions.readState().logRa;
      if (calm) { actions.setLogRa(step.ramp.to); ramp = null; }
      else ramp = { from, to: step.ramp.to, t0: performance.now(), ms: step.ramp.ms };
    } else {
      ramp = null;
    }

    dwellDone = false;
    dwell = step.dwell
      ? { kind: step.dwell, t0: performance.now(), base: actions.steps() }
      : null;

    // Scrolled into view before the first measurement rather than after, so
    // the hole does not animate from wherever the blade happened to be
    // sitting in a scrolled pane to where it ends up.
    if (step.target) {
      actions.element(step.target)?.scrollIntoView({
        block: "nearest", behavior: calm ? "auto" : "smooth",
      });
    }

    // The one container exempted from the wholesale chrome dim, because it
    // is the one getting the spotlight instead. Set as a class rather than
    // inline, so index.html keeps both halves of the rule together.
    for (const el of document.querySelectorAll(HOSTS)) el.classList.remove("tour-lit-host");
    targetOf(step)?.host.classList.add("tour-lit-host");

    render(step);
    // Forced, rather than left to the loop's own change test: `render` has
    // just resized the card, and the placement it was last given belongs to
    // the previous step's hole.
    host = hole = null;
  };

  const go = (delta: number): void => {
    const i = at + delta;
    if (i < 0) return;
    if (i >= steps.length) return close();
    enter(i);
  };

  // ---- opening and closing ----------------------------------------------

  const onKey = (e: KeyboardEvent): void => {
    if (!open()) return;
    const k = e.key;
    // `h` is not the tour's key, but it is `main.ts`'s: it toggles the chrome
    // this overlay is measuring against, and a pane that vanishes mid-step
    // takes the lit control with it. Swallowed here rather than guarded
    // there, so the tour owns the exception it needs and that handler stays
    // the one place "toggle" is decided.
    if (k === "h" || k === "H") { e.stopPropagation(); return; }
    if (k === "Escape") { e.stopPropagation(); close(); return; }
    if (k === "ArrowRight" || k === " " || k === "Enter") { e.stopPropagation(); go(1); return; }
    if (k === "ArrowLeft") { e.stopPropagation(); go(-1); }
  };

  const start = (): void => {
    if (open()) return;
    // The four chrome containers are what a step points at, and a hidden one
    // measures zero — so the tour brings them back rather than lighting a
    // rectangle with nothing in it. See `.chrome-hidden` in index.html.
    document.documentElement.classList.remove("chrome-hidden");
    snapshot = { ...actions.readState() };
    // On for the whole tour, not per step: the field is never dimmed, so this
    // is the only thing that ever is, and turning it on and off between steps
    // would be the chrome flashing rather than the spotlight moving.
    document.documentElement.classList.add("tour-dim-chrome");
    root.classList.add("tour-open");
    // Same two-step as `acknowledgements.ts`: display first, so the next
    // frame's class change has something painted to transition from.
    requestAnimationFrame(() => root.classList.add("tour-lit"));
    window.addEventListener("keydown", onKey, { capture: true });
    enter(0);
    raf = requestAnimationFrame(tick);
  };

  function close(): void {
    if (!open()) return;
    at = -1;
    ramp = null;
    dwell = null;
    host = hole = null;
    cancelAnimationFrame(raf);
    window.removeEventListener("keydown", onKey, { capture: true });
    root.classList.remove("tour-lit");
    document.documentElement.classList.remove("tour-dim-chrome");
    for (const el of document.querySelectorAll(HOSTS)) el.classList.remove("tour-lit-host");
    const done = (e: TransitionEvent): void => {
      if (e.propertyName !== "opacity") return;
      root.classList.remove("tour-open");
      root.removeEventListener("transitionend", done);
    };
    root.addEventListener("transitionend", done);
  }

  chip.addEventListener("click", start);
  back.addEventListener("click", () => go(-1));
  next.addEventListener("click", () => go(1));
  end.addEventListener("click", close);
  restore.addEventListener("click", () => {
    if (snapshot) actions.applyPatch(snapshot);
    actions.resetFocus();
    close();
  });
}
