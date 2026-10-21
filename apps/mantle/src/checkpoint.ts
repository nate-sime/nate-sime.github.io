// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Checkpoints: a running state, tagged under a name, kept in a file, and
 * started from again later.
 *
 * The motivation is internal heating. A heated run has to warm its whole
 * interior until the top boundary carries the heat generated inside on top of
 * what comes in at the bottom, and that is a thermal time — thousands of steps
 * in which only the horizontal mean is doing anything. A checkpoint lets that
 * be paid once.
 *
 * **What a checkpoint holds is T, and ψ beside it.** At infinite Prandtl number
 * T is the only field with a history: ψ is re-solved from it every step. ψ is
 * kept anyway because the Krylov tier warm-starts from it, and
 * `GpuSimulation.writeStateFlat` already takes the pair. Next to the fields go
 * the grid, the domain they were solved on, the whole `State` the pane held,
 * the clock, and the diagnostics of the last poll — among them the heat
 * balance, which is what says whether the tag had come to temperature.
 * Tracers are not in it: a coupled run's composition is state too, and a
 * checkpoint does not yet carry it.
 *
 * Nothing here touches the DOM or the GPU, so the format, the naming, the
 * compatibility rules and the resampling are all tested directly
 * (`tests/checkpoint.test.ts`). `ui/checkpoints.ts` is the library dialog,
 * and `main.ts` reads the fields off the solver and writes them back.
 *
 * ## The file
 *
 *     bytes 0–7   magic "MANTLECK"
 *     u32 LE      format version
 *     u32 LE      header length in bytes
 *     header      UTF-8 JSON: { app, format, entries: [ { …metadata,
 *                   fields: { T: { offset, length }, psi: … | null } } ] }
 *     pad         to a 4-byte boundary
 *     payload     little-endian f32 arrays; offsets in bytes from its start
 *
 * One file holds many tags. The header is plain JSON so that a file can be
 * inspected with `head -c` and nothing else; the fields are raw f32 because
 * that is what the GPU holds — a tag reloads bit for bit.
 */

import type { Stop } from "./colormaps";
import type { Geometry, GeometryKind, RadialWalls, Walls } from "./geometry";
import { PLANETS, planetFor, type PlanetId } from "./planets";
import { cubic } from "./solver/temperature";
import {
  GEOMETRY, PRESETS, RADIAL_WALLS, VISCOSITY, WALLS, defaultState, geometryFor,
  type CustomSurfaceSource, type State,
} from "./ui/presets";

export const CHECKPOINT_EXTENSION = ".mantle";
/** Bumped when a reader of the old format could misread the new one. */
export const CHECKPOINT_FORMAT = 1;
const MAGIC = "MANTLECK";
/** Magic, format version, header length. */
const PREAMBLE = 16;

/**
 * What a field's meaning depends on: the axes it was sampled on. Two runs on
 * the same domain can take each other's T at any resolution; on different
 * domains the numbers are on different points and there is nothing to carry.
 */
export interface Domain {
  readonly kind: GeometryKind;
  readonly walls: Walls;
  readonly radialWalls: RadialWalls;
  readonly lo: number;
  readonly hi: number;
  /** The period φ is solved on — `2L` in a walled box, see `Geometry.span`. */
  readonly span: number;
}

export interface Grid {
  /** ψ spline coefficients, radial × azimuthal. */
  readonly nr: number;
  readonly na: number;
  /** T grid points, radial × azimuthal. */
  readonly gnr: number;
  readonly gna: number;
}

export interface Clock {
  readonly time: number;
  readonly steps: number;
  readonly dt: number;
}

/**
 * The last poll's diagnostics when the tag was taken — a frame or so behind
 * the fields, as everything the readout shows is. NaN where none had landed.
 */
export interface Stats {
  readonly nuInner: number;
  readonly nuOuter: number;
  readonly vrms: number;
  readonly tMax: number;
  /**
   * `Nu_outer − Nu_inner − heatingGap·H`: the residual that goes to zero at a
   * steady state, with internal heating or without (see `heatingGap`).
   */
  readonly heatBalance: number;
}

/**
 * The pane's state as it was tagged. Typed loosely on purpose: a file is
 * input, and `restorePatch` is where any of it is trusted again.
 */
export type TaggedState = Readonly<Record<string, unknown>>;

export interface Checkpoint {
  /** Renamable; unique within a library. */
  name: string;
  /** ISO 8601. */
  readonly created: string;
  /** The build that wrote it — for the reader, not checked. */
  readonly appVersion: string;
  readonly grid: Grid;
  readonly domain: Domain;
  readonly state: TaggedState;
  readonly clock: Clock;
  readonly stats: Stats;
  /** `gnr × gna`, row-major by radius — the layout `GpuSimulation.read("T")` returns. */
  readonly T: Float32Array;
  /** `nr × na`, likewise, or null if the file carried none. */
  readonly psi: Float32Array | null;
}

/** A file that is not a checkpoint library, or not one this build can read. */
export class CheckpointError extends Error {
  override name = "CheckpointError";
}

export const domainOf = (g: Geometry): Domain => ({
  kind: g.kind, walls: g.walls, radialWalls: g.radialWalls, lo: g.lo, hi: g.hi, span: g.span,
});

const near = (a: number, b: number): boolean =>
  Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(a), Math.abs(b));

const fmt = (v: number): string => String(Number(v.toPrecision(3)));

/** A box's drawn width — half the solved period when it is walled. */
const widthOf = (d: Domain): number => (d.walls === "free-slip" ? d.span / 2 : d.span);

const article = (k: GeometryKind): string => (k === "annulus" ? "an annulus" : "a box");
const sides = (w: Walls): string => (w === "free-slip" ? "walled" : "periodic");

/**
 * Why a field solved on `from` cannot seed a run on `to`, or null when it can
 * — phrased to follow "this tag is …".
 *
 * The radial walls are left out on purpose. They constrain ψ — which trial
 * functions exist at the two ends — and not T, whose boundary rows are the
 * Dirichlet data either way; ψ is dropped instead when they differ (see
 * `main.ts`).
 */
export function domainMismatch(from: Domain, to: Domain): string | null {
  if (from.kind !== to.kind) return `on ${article(from.kind)}, not ${article(to.kind)}`;
  if (from.kind === "annulus") {
    if (!near(from.lo, to.lo) || !near(from.hi, to.hi))
      return `on an annulus r ${fmt(from.lo)} … ${fmt(from.hi)}, not r ${fmt(to.lo)} … ${fmt(to.hi)}`;
    return null;
  }
  if (from.walls !== to.walls) return `in a ${sides(from.walls)} box, not a ${sides(to.walls)} one`;
  if (!near(from.span, to.span) || !near(from.lo, to.lo) || !near(from.hi, to.hi))
    return `in a box ${fmt(widthOf(from))} wide, not ${fmt(widthOf(to))}`;
  return null;
}

// ---- reading a tagged state ------------------------------------------------

const str = (s: TaggedState, k: string): string | null =>
  typeof s[k] === "string" ? (s[k] as string) : null;
const num = (s: TaggedState, k: string): number | null =>
  typeof s[k] === "number" && Number.isFinite(s[k]) ? (s[k] as number) : null;
const isPlanetId = (v: unknown): v is PlanetId =>
  typeof v === "string" && Object.hasOwn(PLANETS, v);

/**
 * What to call the domain in a name: the planet whose radii it has, or what
 * it is. The box's width comes from the domain rather than the tagged state,
 * since it is the domain the field is actually on.
 */
export function domainLabel(state: TaggedState, domain: Domain): string {
  if (domain.kind === "box") return `box ×${fmt(widthOf(domain))}`;
  const custom = state.customPlanet;
  if (custom && typeof custom === "object" && typeof (custom as { name?: unknown }).name === "string")
    return (custom as { name: string }).name;
  return isPlanetId(state.activePlanet) ? planetFor(state.activePlanet).label : "annulus";
}

/** `1.0e6`, not `1.0e+6`: a name is read, not parsed. */
const exp1 = (v: number): string => v.toExponential(1).replace("e+", "e");

/**
 * The name a tag gets until it is renamed: where, how vigorous, how heated,
 * and when — `Earth · Ra 2.0e4 · H 20.0 · t 0.153`. H is left out when it is
 * zero; a name is short or it is not read.
 */
export function defaultName(
  state: TaggedState, domain: Domain, time: number, taken: Iterable<string>,
): string {
  const logRa = num(state, "logRa");
  const Ra = state.isothermal === true ? 0 : logRa === null ? NaN : 10 ** logRa;
  const H = num(state, "H") ?? 0;
  const parts = [
    domainLabel(state, domain),
    `Ra ${Ra === 0 ? "0" : exp1(Ra)}`,
    ...(H > 0 ? [`H ${H.toFixed(1)}`] : []),
    `t ${fmt(time)}`,
  ];
  return uniqueName(parts.join(" · "), taken);
}

/** `base`, or `base (2)`, `base (3)`, … — whichever is not already taken. */
export function uniqueName(base: string, taken: Iterable<string>): string {
  const used = new Set(taken);
  if (!used.has(base)) return base;
  for (let n = 2; ; n++) if (!used.has(`${base} (${n})`)) return `${base} (${n})`;
}

/**
 * The one line under a tag's name: what the name leaves out. The heat
 * balance is last and signed — it is the number that says whether the tag
 * had come to temperature.
 */
export function describe(c: Checkpoint): string {
  const law = str(c.state, "viscosity") ?? "unknown viscosity";
  const R = c.stats.heatBalance;
  return `ψ ${c.grid.nr}×${c.grid.na} · ${law} · step ${c.clock.steps.toLocaleString("en-US")} · ` +
    `heat balance ${Number.isNaN(R) ? "—" : (R >= 0 ? "+" : "−") + Math.abs(R).toPrecision(2)}`;
}

// ---- restoring the tagged settings -----------------------------------------

/**
 * The `State` fields `restore` brings back: what the run solved, and how. The
 * view (colour map, speed, overlays, the sun) is the reader's own and stays.
 * The tracer settings stay too, since the tracers themselves are not in a
 * checkpoint — restoring a chemical run's Rb without its composition would
 * claim a coupling the field never had.
 */
export const RESTORE_KEYS = [
  "activePlanet", "customPlanet", "geometry", "boxLength", "walls", "radialWalls",
  "resolution", "logRa", "isothermal", "H", "dtMax", "courant",
  "viscosity", "logContrast", "logDepthContrast", "iters", "n", "picard",
  "sigmaY", "sigmaB", "etaStar", "etaLight", "etaDense",
] as const satisfies readonly (keyof State)[];
export type RestoreKey = typeof RESTORE_KEYS[number];

const isName = (table: object) => (v: unknown): boolean =>
  typeof v === "string" && Object.hasOwn(table, v);
const isNum = (v: unknown): boolean => typeof v === "number" && Number.isFinite(v);
const isPos = (v: unknown): boolean => isNum(v) && (v as number) > 0;
const SURFACES = ["procedural", "earth-daymap", "venus-magellan", "mars-viking"] as const satisfies
  readonly CustomSurfaceSource[];
const SURFACE_NUMBERS = [
  "seed", "rockiness", "terrainScale", "oceanCoverage", "plantLife", "iceCaps",
  "cloudCover", "atmosphereDensity",
] as const;

const isCustomPlanet = (v: unknown): boolean => {
  if (!v || typeof v !== "object") return false;
  const p = v as Record<string, unknown>;
  const s = p.surface as Record<string, unknown> | null;
  return typeof p.name === "string" && isPos(p.innerRadiusKm) && isPos(p.outerRadiusKm)
    && (p.outerRadiusKm as number) > (p.innerRadiusKm as number)
    && (SURFACES as readonly unknown[]).includes(p.surfaceSource)
    && !!s && typeof s === "object" && s.kind === "procedural"
    && SURFACE_NUMBERS.every((k) => isNum(s[k])) && typeof s.atmosphereHue === "string";
};

const VALID: Record<RestoreKey, (v: unknown) => boolean> = {
  activePlanet: (v) => v === null || isPlanetId(v),
  customPlanet: (v) => v === null || isCustomPlanet(v),
  geometry: isName(GEOMETRY),
  boxLength: isPos,
  walls: isName(WALLS),
  radialWalls: isName(RADIAL_WALLS),
  resolution: isName(PRESETS),
  logRa: isNum,
  isothermal: (v) => typeof v === "boolean",
  H: isNum,
  dtMax: isPos,
  courant: isPos,
  viscosity: isName(VISCOSITY),
  logContrast: isNum,
  logDepthContrast: isNum,
  iters: isPos,
  n: isPos,
  picard: isPos,
  sigmaY: isNum,
  sigmaB: isNum,
  etaStar: isPos,
  etaLight: isPos,
  etaDense: isPos,
};

/**
 * The settings `restore` would apply, or why it cannot.
 *
 * A field the tag predates takes its default rather than whatever the pane
 * holds now: it did not exist when the run was tagged, so the run had the
 * behaviour its default stands for (internal heating, say, was 0 before there
 * was a slider for it).
 *
 * The settings must also reproduce the domain the field was solved on. They
 * normally do by construction, but a planet's radii or a preset could change
 * between the build that wrote a file and the one reading it, and restoring
 * the settings onto a domain the field does not fit would be a restore in
 * name only.
 */
export function restorePatch(
  c: Checkpoint, current: State,
): { patch: Partial<State> } | { problem: string } {
  const defaults = defaultState();
  const patch: Record<string, unknown> = {};
  for (const k of RESTORE_KEYS) {
    const v = k in c.state ? c.state[k] : defaults[k];
    if (!VALID[k](v)) return { problem: `its ${k} (${JSON.stringify(v)}) is not one this version offers` };
    patch[k] = v !== null && typeof v === "object" ? structuredClone(v) : v;
  }
  let geom: Geometry;
  try { geom = geometryFor({ ...current, ...patch }); }
  catch (e) { return { problem: e instanceof Error ? e.message : String(e) }; }
  const why = domainMismatch(c.domain, domainOf(geom));
  if (why) return { problem: `its settings no longer give the domain it was solved on: it is ${why}` };
  return { patch: patch as Partial<State> };
}

/** The keys of `patch` whose values differ from `current` — compared as data, not by identity. */
export function changedKeys(patch: Partial<State>, current: State): (keyof State)[] {
  return (Object.keys(patch) as (keyof State)[])
    .filter((k) => JSON.stringify(patch[k]) !== JSON.stringify(current[k]));
}

// ---- resampling --------------------------------------------------------------

/**
 * Carry a T grid onto another resolution of the same domain — what lets a run
 * come to temperature on the coarse grid, at a quarter of the finest grid's
 * dt cap and a sixteenth of its cells, and be reseeded onto a fine one.
 *
 * The stencil is the advection's own — Catmull–Rom in each direction,
 * clamped to the two nodes it sits between (`cubic`, as `Temperature.sample`
 * uses it), periodic in φ — with one change at the radial ends. The sampler
 * repeats the boundary row as the row beyond it, which halves the slope there
 * and makes the first cell first-order accurate; that is the boundary layer,
 * the one place a coarse-to-fine reseed most needs right. Here the row beyond
 * is extrapolated linearly instead. Every sample still lies between the two
 * real nodes either side of it, so the result stays within the source's
 * range — no new extrema, which is the property the advection's limiter would
 * otherwise have to undo.
 *
 * It is written in index space rather than through `sample`'s physical
 * coordinates so that a node coinciding with a source node lands on it
 * exactly: the boundary rows come through as the Dirichlet data they were,
 * and a same-sized grid comes back bit for bit.
 */
export function resample(
  src: Float32Array, srcNr: number, srcNa: number, nr: number, na: number,
): Float32Array {
  if (srcNr === nr && srcNa === na) return src.slice();
  const out = new Float32Array(nr * na);
  const node = (i: number, j: number): number => src[i * srcNa + (((j % srcNa) + srcNa) % srcNa)];
  const at = (i: number, j: number): number =>
    i < 0 ? 2 * node(0, j) - node(1, j)
      : i > srcNr - 1 ? 2 * node(srcNr - 1, j) - node(srcNr - 2, j)
        : node(i, j);
  const col = new Float64Array(4);
  for (let i = 0; i < nr; i++) {
    // Integer numerators: exact, so a shared node is hit exactly.
    const x = nr > 1 ? (i * (srcNr - 1)) / (nr - 1) : 0;
    const i0 = Math.floor(x), tr = x - i0;
    for (let j = 0; j < na; j++) {
      const y = (j * srcNa) / na;
      const j0 = Math.floor(y), tp = y - j0;
      for (let m = -1; m <= 2; m++)
        col[m + 1] = cubic(at(i0 + m, j0 - 1), at(i0 + m, j0), at(i0 + m, j0 + 1), at(i0 + m, j0 + 2), tp);
      out[i * na + j] = cubic(col[0], col[1], col[2], col[3], tr);
    }
  }
  return out;
}

// ---- thumbnails --------------------------------------------------------------

/**
 * A `size × size` RGBA picture of a tag's T, for the library's list: the
 * annulus or the box as the canvas draws it (hot boundary inward or at the
 * bottom, φ anticlockwise from +x), coloured over [0, max T] as the canvas
 * colours it, transparent outside the domain. Nearest-neighbour — it is a
 * way to tell tags apart, not a figure.
 */
export function thumbnail(
  c: Pick<Checkpoint, "T" | "grid" | "domain">, size: number, stops: readonly Stop[],
): Uint8ClampedArray<ArrayBuffer> {
  const { T, grid: { gnr, gna }, domain: d } = c;
  let top = 1;
  for (let k = 0; k < T.length; k++) if (T[k] > top) top = T[k];
  const px = new Uint8ClampedArray(size * size * 4);
  const paint = (p: number, v: number): void => {
    const s = Math.min(1, Math.max(0, v / top)) * (stops.length - 1);
    const k = Math.min(stops.length - 2, Math.floor(s)), f = s - k;
    for (let ch = 0; ch < 3; ch++)
      px[p + ch] = 255 * (stops[k][ch] + f * (stops[k + 1][ch] - stops[k][ch]));
    px[p + 3] = 255;
  };
  const depth = d.hi - d.lo;
  const row = (r: number): number => Math.round(((r - d.lo) / depth) * (gnr - 1));
  if (d.kind === "annulus") {
    const h = d.hi;
    for (let y = 0; y < size; y++)
      for (let x = 0; x < size; x++) {
        const wx = ((x + 0.5) / size) * 2 * h - h, wy = h - ((y + 0.5) / size) * 2 * h;
        const r = Math.hypot(wx, wy);
        if (r < d.lo || r > d.hi) continue;
        let phi = Math.atan2(wy, wx);
        if (phi < 0) phi += 2 * Math.PI;
        const j = Math.floor((phi / d.span) * gna) % gna;
        paint(4 * (y * size + x), T[row(r) * gna + j]);
      }
    return px;
  }
  const width = widthOf(d), ext = Math.max(width, depth);
  for (let y = 0; y < size; y++)
    for (let x = 0; x < size; x++) {
      const wx = ((x + 0.5) / size) * ext - (ext - width) / 2;
      const wz = d.lo + ext - ((y + 0.5) / size) * ext - (ext - depth) / 2;
      if (wx < 0 || wx > width || wz < d.lo || wz > d.hi) continue;
      const j = Math.min(gna - 1, Math.floor((wx / d.span) * gna));
      paint(4 * (y * size + x), T[row(wz) * gna + j]);
    }
  return px;
}

// ---- the file ----------------------------------------------------------------

interface FieldRef { offset: number; length: number }

/** Every tag, as one file — see the module header for the layout. */
export function encodeCheckpoints(entries: readonly Checkpoint[]): Uint8Array<ArrayBuffer> {
  let bytes = 0;
  const place = (a: Float32Array): FieldRef => {
    const ref = { offset: bytes, length: a.length };
    bytes += 4 * a.length;
    return ref;
  };
  const meta = entries.map((c) => ({
    name: c.name, created: c.created, appVersion: c.appVersion,
    grid: c.grid, domain: c.domain, state: c.state, clock: c.clock,
    // JSON has no NaN; `null` reads back as one.
    stats: c.stats,
    fields: { T: place(c.T), psi: c.psi ? place(c.psi) : null },
  }));
  const header = new TextEncoder().encode(
    JSON.stringify({ app: "mantle", format: CHECKPOINT_FORMAT, entries: meta }));
  const start = Math.ceil((PREAMBLE + header.length) / 4) * 4;
  const out = new Uint8Array(start + bytes);
  const view = new DataView(out.buffer);
  for (let k = 0; k < MAGIC.length; k++) out[k] = MAGIC.charCodeAt(k);
  view.setUint32(8, CHECKPOINT_FORMAT, true);
  view.setUint32(12, header.length, true);
  out.set(header, PREAMBLE);
  // Through a DataView rather than a Float32Array over the buffer: that would
  // write in the platform's byte order, and the file states little-endian.
  const write = (a: Float32Array, ref: FieldRef): void => {
    for (let k = 0; k < a.length; k++) view.setFloat32(start + ref.offset + 4 * k, a[k], true);
  };
  entries.forEach((c, n) => {
    write(c.T, meta[n].fields.T);
    if (c.psi && meta[n].fields.psi) write(c.psi, meta[n].fields.psi);
  });
  return out;
}

const fail = (msg: string): never => { throw new CheckpointError(msg); };
const obj = (v: unknown, what: string): Record<string, unknown> =>
  v && typeof v === "object" && !Array.isArray(v) ? v as Record<string, unknown> : fail(`${what} is missing`);
const count = (v: unknown, what: string): number =>
  Number.isSafeInteger(v) && (v as number) > 0 ? v as number : fail(`${what} is not a positive whole number`);
const finite = (v: unknown, what: string): number =>
  isNum(v) ? v as number : fail(`${what} is not a number`);
const oneOf = <T extends string>(v: unknown, options: readonly T[], what: string): T =>
  (options as readonly unknown[]).includes(v) ? v as T : fail(`${what} is not one of ${options.join(", ")}`);
const statOf = (v: unknown): number => (isNum(v) ? v as number : NaN);

/**
 * Read a library back. Throws `CheckpointError`, with a message fit to show
 * the reader, for anything that is not a library this build can read — a
 * file is input, and a field of the wrong length or with a NaN in it would
 * otherwise reach the GPU.
 */
export function decodeCheckpoints(data: ArrayBuffer | Uint8Array): Checkpoint[] {
  const bytes = data instanceof Uint8Array ? data : new Uint8Array(data);
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (bytes.length < PREAMBLE || String.fromCharCode(...bytes.subarray(0, MAGIC.length)) !== MAGIC)
    fail("it is not a mantle checkpoint file");
  const format = view.getUint32(8, true);
  if (format > CHECKPOINT_FORMAT)
    fail(`it was written by a newer version of the app (format ${format}; this one reads ${CHECKPOINT_FORMAT})`);
  if (format < 1) fail(`its format version (${format}) is not one this app wrote`);
  const headerLength = view.getUint32(12, true);
  if (PREAMBLE + headerLength > bytes.length) fail("it is truncated");
  let header: unknown;
  try {
    header = JSON.parse(new TextDecoder("utf-8", { fatal: true })
      .decode(bytes.subarray(PREAMBLE, PREAMBLE + headerLength)));
  } catch { fail("its header is not readable"); }
  const entries = obj(header, "the header").entries;
  if (!Array.isArray(entries)) fail("its header lists no tags");
  const start = Math.ceil((PREAMBLE + headerLength) / 4) * 4;

  return (entries as unknown[]).map((raw, n) => {
    const e = obj(raw, `tag ${n + 1}`);
    const name = typeof e.name === "string" && e.name.trim() ? e.name.trim() : fail(`tag ${n + 1} has no name`);
    const where = `tag "${name}"`;
    const g = obj(e.grid, `${where}'s grid`);
    const grid: Grid = {
      nr: count(g.nr, `${where}'s grid`), na: count(g.na, `${where}'s grid`),
      gnr: count(g.gnr, `${where}'s grid`), gna: count(g.gna, `${where}'s grid`),
    };
    const d = obj(e.domain, `${where}'s domain`);
    const domain: Domain = {
      kind: oneOf(d.kind, ["annulus", "box"], `${where}'s domain`),
      walls: oneOf(d.walls, ["periodic", "free-slip"], `${where}'s side walls`),
      radialWalls: oneOf(d.radialWalls, ["free-slip", "no-slip"], `${where}'s radial walls`),
      lo: finite(d.lo, `${where}'s domain`), hi: finite(d.hi, `${where}'s domain`),
      span: finite(d.span, `${where}'s domain`),
    };
    if (!(domain.hi > domain.lo) || !(domain.span > 0)) fail(`${where}'s domain is empty`);
    const k = obj(e.clock, `${where}'s clock`);
    const clock: Clock = {
      time: finite(k.time, `${where}'s time`),
      steps: Number.isSafeInteger(k.steps) && (k.steps as number) >= 0
        ? k.steps as number : fail(`${where}'s step count is not a whole number`),
      dt: isPos(k.dt) ? k.dt as number : fail(`${where}'s dt is not a positive number`),
    };
    const s = e.stats && typeof e.stats === "object" ? e.stats as Record<string, unknown> : {};
    const stats: Stats = {
      nuInner: statOf(s.nuInner), nuOuter: statOf(s.nuOuter), vrms: statOf(s.vrms),
      tMax: statOf(s.tMax), heatBalance: statOf(s.heatBalance),
    };
    const f = obj(e.fields, `${where}'s fields`);
    const field = (ref: unknown, length: number, what: string): Float32Array => {
      const r = obj(ref, `${where}'s ${what}`);
      const offset = r.offset as number;
      if (!Number.isSafeInteger(offset) || offset < 0 || offset % 4 !== 0 || r.length !== length)
        fail(`${where}'s ${what} does not match its grid`);
      if (start + offset + 4 * length > bytes.length) fail("it is truncated");
      const a = new Float32Array(length);
      for (let i = 0; i < length; i++) {
        a[i] = view.getFloat32(start + offset + 4 * i, true);
        if (!Number.isFinite(a[i])) fail(`${where}'s ${what} has a value that is not a number`);
      }
      return a;
    };
    return {
      name,
      created: typeof e.created === "string" ? e.created : "",
      appVersion: typeof e.appVersion === "string" ? e.appVersion : "",
      grid, domain,
      state: obj(e.state, `${where}'s settings`),
      clock, stats,
      T: field(f.T, grid.gnr * grid.gna, "temperature"),
      psi: f.psi === null || f.psi === undefined ? null : field(f.psi, grid.nr * grid.na, "stream function"),
    };
  });
}
