// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The text under the figure. A view hands back plain preformatted lines, or a
 * `Rich` readout: an optional run-status strip, then sections — side by side
 * when there are two, as the beam's beside the plate's — each a run of
 * blocks: a lead line, stat tiles, tables, warnings, muted notes.
 *
 * Views redraw every frame while something moves, so the rich readout is
 * patched into the page node by node rather than rebuilt: a tile's hover text
 * (what the number is) stays open while its number changes under it, and an
 * unchanged table is not touched at all.
 */

export type Tone = "good" | "warn" | "muted";

export interface Tile {
  readonly label: string;
  /** The number, large. */
  readonly value: string;
  /** Its interval, or a second number that qualifies it. */
  readonly detail?: string;
  /** What the number says, with a coloured dot. */
  readonly verdict?: { readonly text: string; readonly tone: Tone };
  /** What the number is, on hover. */
  readonly help?: string;
}

export interface Table {
  readonly caption?: string;
  /** Column heads; the first column is the row's name, the rest are numbers. */
  readonly head: readonly string[];
  readonly rows: readonly (readonly string[])[];
  /** A row to pick out: the mode on screen. */
  readonly mark?: number;
}

export type Block =
  | { readonly kind: "lead" | "note" | "warn"; readonly text: string }
  | { readonly kind: "tiles"; readonly tiles: readonly Tile[]; readonly cols: number }
  | { readonly kind: "table"; readonly table: Table };

export interface Section {
  /** Small capitals above it: BEAM, PLATE. */
  readonly title?: string;
  readonly blocks: readonly Block[];
}

export interface Status {
  /** running, done, paused, waiting, stopped: … */
  readonly state: string;
  readonly done: number;
  readonly total: number;
  /** Counts and rates, after the bar. */
  readonly text: string;
}

export interface Rich {
  readonly status?: Status;
  readonly sections: readonly Section[];
}

/** Blocks, by what they say: the setup, an aside, a caveat, numbers, a table. */
export const ro = {
  lead: (text: string): Block => ({ kind: "lead", text }),
  note: (text: string): Block => ({ kind: "note", text }),
  warn: (text: string): Block => ({ kind: "warn", text }),
  tiles: (tiles: readonly Tile[], cols = 3): Block => ({ kind: "tiles", tiles, cols }),
  table: (head: readonly string[], rows: readonly (readonly string[])[], o: { caption?: string; mark?: number } = {}): Block =>
    ({ kind: "table", table: { head, rows, ...o } }),
};

/** A view's readout as one titled section: a rich one's blocks, or plain lines as notes (a view that cannot draw says why). */
export function section(title: string, r: string | Rich): Section {
  if (typeof r === "string") return { title, blocks: r.split("\n").filter((l) => l).map(ro.note) };
  return { title, blocks: r.sections.flatMap((s) => s.blocks) };
}

// ---- drawing: a small virtual tree, patched into the page ----

interface V {
  readonly t: string;
  readonly c?: string;
  readonly a?: Readonly<Record<string, string>>;
  readonly k?: readonly Kid[];
}
type Kid = V | string;

const h = (t: string, c: string | undefined, k: readonly Kid[] = [], a?: Record<string, string>): V => ({ t, c, k, a });

function tileV(t: Tile): V {
  return h("div", "ro-tile", [
    h("div", "ro-label", [t.label]),
    h("div", "ro-value", [t.value]),
    h("div", "ro-detail", [t.detail ?? ""]),
    h("div", "ro-verdict", [t.verdict?.text ?? ""], t.verdict ? { "data-tone": t.verdict.tone } : undefined),
  ], t.help ? { title: t.help } : undefined);
}

function tableV(t: Table): V {
  const table = h("table", undefined, [
    h("thead", undefined, [h("tr", undefined, t.head.map((s) => h("th", undefined, [s])))]),
    h("tbody", undefined, t.rows.map((r, i) => h("tr", i === t.mark ? "mark" : undefined, r.map((s) => h("td", undefined, [s]))))),
  ]);
  return h("div", "ro-table", t.caption ? [h("div", "ro-caption", [t.caption]), table] : [table]);
}

function blockV(b: Block): V {
  switch (b.kind) {
    case "tiles": return h("div", "ro-tiles", b.tiles.map(tileV), { style: `--tile-cols: ${b.cols}` });
    case "table": return tableV(b.table);
    default: return h("div", `ro-${b.kind}`, [b.text]);
  }
}

function richV(r: Rich): Kid[] {
  const out: Kid[] = [];
  if (r.status) {
    const f = r.status.total > 0 ? Math.min(1, r.status.done / r.status.total) : 0;
    out.push(h("div", "ro-status", [
      h("div", "ro-state", [r.status.state], { "data-state": r.status.state.split(/[ :]/)[0] }),
      h("div", "ro-bar", [h("span", undefined, [], { style: `width: ${(100 * f).toFixed(1)}%` })]),
      h("div", "ro-text", [r.status.text]),
    ]));
  }
  out.push(h("div", "ro-sections", r.sections.map((s) =>
    h("div", "ro-section", [...(s.title ? [h("div", "ro-title", [s.title])] : []), ...s.blocks.map(blockV)])),
  { style: `--sections: ${r.sections.length}` }));
  return out;
}

const put = (parent: Node, node: Node, old: ChildNode | undefined) => (old ? parent.replaceChild(node, old) : parent.appendChild(node));

/** Make `parent`'s children match `kids`, touching only what differs. */
function patch(parent: Element, kids: readonly Kid[]): void {
  kids.forEach((k, i) => {
    const old = parent.childNodes[i];
    if (typeof k === "string") {
      if (old?.nodeType === Node.TEXT_NODE) { if (old.nodeValue !== k) old.nodeValue = k; }
      else put(parent, document.createTextNode(k), old);
      return;
    }
    let e = old instanceof Element && old.localName === k.t ? old : null;
    if (!e) put(parent, (e = document.createElement(k.t)), old);
    const want: Record<string, string> = { ...k.a, ...(k.c ? { class: k.c } : {}) };
    for (const n of e.getAttributeNames()) if (!(n in want)) e.removeAttribute(n);
    for (const [n, v] of Object.entries(want)) if (e.getAttribute(n) !== v) e.setAttribute(n, v);
    patch(e, k.k ?? []);
  });
  while (parent.childNodes.length > kids.length) parent.lastChild!.remove();
}

export function showReadout(el: HTMLElement, r: string | Rich): void {
  if (typeof r === "string") {
    el.classList.remove("rich");
    if (el.firstElementChild || el.textContent !== r) el.textContent = r;
    return;
  }
  el.classList.add("rich");
  patch(el, richV(r));
}
