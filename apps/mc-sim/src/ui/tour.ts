// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The guided tour overlay, opened from the pane's "guided tours" folder. The
 * steps are data (`tours.ts`); this file makes one true and draws it.
 *
 * Ported from the mantle app's tour, without its solver-specific parts (the
 * Rayleigh-number ramps, the camera, the globe), and keeping its rules:
 *
 * - **The figure and the readout are never dimmed.** They are what a step asks
 *   the reader to watch. Only the pane goes down — and in it, the one control
 *   a step is about stays lit and live.
 * - **The spotlight is four rectangles tiled between the lit control and the
 *   pane**, not one element with a hole cut in it: `#pane` has no stacking
 *   context of its own, so raising one blade above a covering sheet would
 *   raise the whole pane. Four shades around the control never cover it, and
 *   they take pointer events, so inside the pane only that control answers.
 * - **A dwell nudges, it never advances.** It fills a bar — over wall time, or
 *   until the active run has so many samples — and then lights "next". The
 *   reader moves on with two buttons, or the arrow keys; Escape leaves.
 * - **The tour does not undo itself.** It ends where its last step left the
 *   app, which is the useful place to be; the last card offers the settings
 *   the reader had before it began.
 *
 * Everything is re-measured every frame while the tour is open, which is what
 * keeps the spotlight on its control through a resize or a scroll of the pane
 * without subscribing to either.
 */

import type { State } from "./state";
import { TOURS, type TourDwell, type TourName, type TourStep } from "./tours";
import type { ControlName } from "./visibility";

/** What the tour needs from the rest of the app — passed in, so this file holds no reference to the pane or the pool. */
export interface TourActions {
  element(name: ControlName): HTMLElement;
  /** Set fields of `State`, refresh the pane, redraw. */
  applyPatch(patch: Partial<State>): void;
  readState(): Readonly<State>;
  /** Samples the active Monte Carlo run has folded in, over all its levels; 0 with none. */
  samples(): number;
}

const PAD = 6;
const GAP = 14;
const MARGIN = 12;
const NARROW = 640;

interface Box { x: number; y: number; w: number; h: number }
const same = (a: Box, b: Box) => a.x === b.x && a.y === b.y && a.w === b.w && a.h === b.h;
const px = (v: number) => `${Math.round(v)}px`;
const setBox = (el: HTMLElement, b: Box) => {
  el.style.left = px(b.x); el.style.top = px(b.y);
  el.style.width = px(Math.max(0, b.w)); el.style.height = px(Math.max(0, b.h));
};
const div = (cls: string, parent: HTMLElement) => {
  const e = document.createElement("div");
  e.className = cls;
  parent.append(e);
  return e;
};
const button = (label: string, parent: HTMLElement) => {
  const b = document.createElement("button");
  b.type = "button";
  b.className = "tour-btn";
  b.textContent = label;
  parent.append(b);
  return b;
};

/** Builds the overlay into `#tour` and returns the function that starts a tour. */
export function buildTour(root: HTMLElement, pane: HTMLElement, actions: TourActions): (name: TourName) => void {
  const shades = ["top", "bottom", "left", "right"].map((side) => div(`tour-shade tour-${side}`, root));
  const ring = div("tour-ring", root);
  const card = div("tour-card", root);
  card.setAttribute("role", "dialog");
  card.setAttribute("aria-live", "polite");
  const count = div("tour-count", card);
  const title = document.createElement("h2");
  title.className = "tour-title";
  card.append(title);
  const body = div("tour-body", card);
  const watch = div("tour-watch", card);
  const bar = div("tour-bar", card);
  const fill = div("tour-bar-fill", bar);
  const restoreRow = div("tour-restore", card);
  const restore = button("restore the settings I had", restoreRow);
  const nav = div("tour-nav", card);
  const end = button("end tour", nav);
  end.classList.add("tour-btn-quiet");
  div("tour-spacer", nav);
  const back = button("back", nav);
  const next = button("next", nav);
  next.classList.add("tour-btn-primary");

  let steps: readonly TourStep[] = [];
  let at = -1;
  let raf = 0;
  let host: Box | null = null, hole: Box | null = null;
  let dwell: { kind: TourDwell; t0: number } | null = null;
  let snapshot: State | null = null;
  const open = () => at >= 0;

  const lit = (step: TourStep): HTMLElement | null => {
    if (!step.target) return null;
    const el = actions.element(step.target);
    // A hidden blade measures as an empty box at the corner; no spotlight then.
    return el.getClientRects().length ? el : null;
  };

  const paint = (h: Box, o: Box, outline: boolean) => {
    const [top, bottom, left, right] = shades;
    setBox(top, { x: h.x, y: h.y, w: h.w, h: o.y - h.y });
    setBox(bottom, { x: h.x, y: o.y + o.h, w: h.w, h: h.y + h.h - o.y - o.h });
    setBox(left, { x: h.x, y: o.y, w: o.x - h.x, h: o.h });
    setBox(right, { x: o.x + o.w, y: o.y, w: h.x + h.w - o.x - o.w, h: o.h });
    setBox(ring, o);
    ring.style.opacity = outline ? "1" : "0";
  };

  /**
   * Beside the lit control, on the side with room — the pane is a strip down
   * the right, so nearly always to its left. With nothing lit, bottom right of
   * the figure's column: over the readout's empty right-hand end rather than
   * the plots or the numbers, which are left-aligned.
   */
  const place = (o: Box | null) => {
    const vw = window.innerWidth, vh = window.innerHeight, { width: cw, height: ch } = card.getBoundingClientRect();
    const paneLeft = pane.getBoundingClientRect().left;
    let x: number, y: number;
    if (vw < NARROW) { x = (vw - cw) / 2; y = vh - ch - MARGIN; }
    else if (!o) { x = paneLeft - cw - MARGIN; y = vh - ch - MARGIN; }
    else if (o.x - GAP - cw >= MARGIN) { x = o.x - GAP - cw; y = o.y - ch / 3; }
    else { x = o.x + o.w / 2 - cw / 2; y = o.y + o.h + GAP; }
    card.style.left = px(Math.min(vw - cw - MARGIN, Math.max(MARGIN, x)));
    card.style.top = px(Math.min(vh - ch - MARGIN, Math.max(MARGIN, y)));
  };

  const tick = (now: number) => {
    if (!open()) return;
    const step = steps[at], el = lit(step);
    const r = pane.getBoundingClientRect();
    const nextHost: Box = { x: r.left, y: r.top, w: r.width, h: r.height };
    let nextHole: Box = { x: nextHost.x, y: nextHost.y, w: 0, h: 0 };
    if (el) {
      const b = el.getBoundingClientRect();
      const x0 = Math.max(r.left, b.left - PAD), y0 = Math.max(r.top, b.top - PAD);
      const x1 = Math.min(r.right, b.right + PAD), y1 = Math.min(r.bottom, b.bottom + PAD);
      nextHole = { x: x0, y: y0, w: Math.max(0, x1 - x0), h: Math.max(0, y1 - y0) };
    }
    if (!host || !hole || !same(host, nextHost) || !same(hole, nextHole)) {
      host = nextHost;
      hole = nextHole;
      // Nothing lit: the hole collapses to the pane's corner and all four shades cover it.
      paint(nextHost, nextHole, el !== null);
      place(el ? nextHole : null);
    }
    if (dwell) {
      const p = "samples" in dwell.kind ? actions.samples() / dwell.kind.samples : (now - dwell.t0) / dwell.kind.ms;
      const c = Math.min(1, Math.max(0, p));
      fill.style.width = `${(100 * c).toFixed(1)}%`;
      if (c >= 1) { next.classList.add("tour-btn-ready"); dwell = null; }
    }
    raf = requestAnimationFrame(tick);
  };

  const render = (step: TourStep) => {
    count.textContent = `${at + 1} / ${steps.length}`;
    title.textContent = step.title;
    body.replaceChildren(...step.body.map((t) => {
      if (typeof t === "string") {
        const p = document.createElement("p");
        p.textContent = t;
        return p;
      }
      const ul = document.createElement("ul");
      ul.className = "tour-list";
      ul.append(...t.items.map((item) => {
        const li = document.createElement("li");
        li.textContent = item;
        return li;
      }));
      return ul;
    }));
    watch.textContent = step.watch ?? "";
    watch.style.display = step.watch ? "" : "none";
    bar.style.display = step.dwell ? "" : "none";
    fill.style.width = "0%";
    back.disabled = at === 0;
    const last = at === steps.length - 1;
    next.textContent = last ? "finish" : "next";
    next.classList.toggle("tour-btn-ready", last);
    restoreRow.style.display = last && snapshot ? "" : "none";
  };

  const enter = (i: number) => {
    at = i;
    const step = steps[i];
    if (step.patch) actions.applyPatch(step.patch);
    dwell = step.dwell ? { kind: step.dwell, t0: performance.now() } : null;
    const el = lit(step);
    el?.scrollIntoView({ block: "nearest" });
    pane.classList.toggle("tour-lit-host", el !== null);
    render(step);
    host = hole = null;
  };

  /** Walking backwards replays every patch up to the step, so a step always runs in the state it was written for. */
  const go = (delta: number) => {
    const i = at + delta;
    if (i < 0) return;
    if (i >= steps.length) { close(); return; }
    if (delta < 0) for (let k = 0; k < i; k++) if (steps[k].patch) actions.applyPatch(steps[k].patch!);
    enter(i);
  };

  const onKey = (e: KeyboardEvent) => {
    if (!open()) return;
    const k = e.key;
    if (k === "Escape") { e.stopPropagation(); close(); }
    else if (k === "ArrowRight") { e.stopPropagation(); go(1); }
    else if (k === "ArrowLeft") { e.stopPropagation(); go(-1); }
  };

  function close(): void {
    if (!open()) return;
    at = -1;
    dwell = null;
    host = hole = null;
    cancelAnimationFrame(raf);
    window.removeEventListener("keydown", onKey, { capture: true });
    root.classList.remove("tour-lit");
    document.documentElement.classList.remove("tour-dim-chrome");
    pane.classList.remove("tour-lit-host");
    const done = (e: TransitionEvent) => {
      if (e.propertyName !== "opacity") return;
      root.classList.remove("tour-open");
      root.removeEventListener("transitionend", done);
    };
    root.addEventListener("transitionend", done);
  }

  back.addEventListener("click", () => go(-1));
  next.addEventListener("click", () => go(1));
  end.addEventListener("click", close);
  restore.addEventListener("click", () => {
    if (snapshot) actions.applyPatch(snapshot);
    close();
  });

  return (name: TourName) => {
    if (open()) return;
    steps = TOURS[name];
    snapshot = { ...actions.readState() };
    document.documentElement.classList.add("tour-dim-chrome");
    root.classList.add("tour-open");
    requestAnimationFrame(() => root.classList.add("tour-lit"));
    window.addEventListener("keydown", onKey, { capture: true });
    enter(0);
    raf = requestAnimationFrame(tick);
  };
}
