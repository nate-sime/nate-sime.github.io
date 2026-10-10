// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 2 on screen, the left column of view 2 (`structures.ts`): one beam,
 * solved — its static deflection under the chosen load above, and below,
 * swinging, one of its lowest modes or (stage 8) its steady response to the
 * load applied harmonically at Ω.
 *
 * Deflection is drawn growing downward, the way a loaded beam sags; a mode is
 * normalised to unit amplitude and animated as φ(x) cos ωt at one visual rate
 * for every mode (the real rates differ by orders of magnitude, so the
 * readout carries them instead). The control polygon — the coefficients at
 * their Greville abscissae — is the spline's own picture of the curve, and
 * hugs it more tightly as the mesh refines; the mesh, as ticks on the axis,
 * shows where its pieces join.
 *
 * The forced response u(x) is complex: Rayleigh damping puts every point a
 * little out of phase with the load. What swings is Re(u e^{iΩt}) =
 * Re u cos Ωt − Im u sin Ωt, inside the envelope ±|u|. Near a resonance |u|
 * dwarfs the static deflection (drawn dashed for scale), and the shape tends
 * to the resonant mode's.
 */

import { Beam, SUPPORTS, admissible, type End } from "../../beam/beam";
import { exactEigenvalues } from "../../beam/exact";
import { pinnedResponse } from "../../beam/exact";
import { evaluateQoI, exactStaticOf, harmonicOf, loadOf, qoiPoint, sectionOf, type Harmonic } from "../../hierarchy";
import { meshTicks } from "../heatmap";
import { greville } from "../../spline";
import { SLOT, type Axes, type Plot, type Series } from "../plot";
import { beamCaseOf, continuityK, continuityName, displayOf, type State } from "../state";
import { axisUnit, deflectionScale, fmt, plain, referenceNote, type Display } from "../units";
import { PLATE_MODES } from "./plate";
import { ro, type Block } from "../readout";
import { memo, type ViewResult } from "./view";

/** As many modes as the plate offers: one mode slider drives both. */
const MODES = PLATE_MODES;
const PERIOD_MS = 1600;
const EXACT_INK = "rgba(207, 238, 255, 0.85)";

interface Solved {
  beam: Beam;
  c: Float64Array;
  F: Float64Array;
  modes: { values: Float64Array; vectors: Float64Array[] };
  ms: number;
}

const solved = memo<Solved | string>();

interface Response {
  h: Harmonic;
  re: Float64Array;
  im: Float64Array;
}

const responses = memo<Response>();
const MESH_INK = "rgba(207, 238, 255, 0.55)";

/** The beam into two long, flat panels: its static deflection, and its motion. */
export function renderBeam([top, bottom]: readonly [Plot, Plot], st: State, t: number): ViewResult {
  const p = st.p, k = continuityK(st.continuity, p);
  const spec = { p, k, ne: st.ne, supports: SUPPORTS[st.supports], ...sectionOf(st.section) };
  const bc = { supports: st.supports, load: st.load, section: st.section };
  const why = admissible(spec);
  const r = solved(JSON.stringify([p, k, st.ne, st.supports, st.section, st.load]), () => {
    if (why) return why;
    const t0 = performance.now();
    const beam = new Beam(spec);
    const F = beam.load(loadOf(bc));
    const c = beam.solve(F);
    const modes = beam.modes(Math.min(MODES, beam.dofs));
    return { beam, c, F, modes, ms: performance.now() - t0 };
  });
  if (typeof r === "string") {
    for (const plot of [top, bottom]) plot.draw({ xlabel: "x / L", ylabel: "w", xlim: [0, 1], ylim: [-1, 1], series: [] });
    return { readout: `${continuityName(k)}, degree ${p}: not a beam discretisation.\n${r}`, animate: false };
  }

  // This view is the beam's whatever the pane's structure is.
  const d = displayOf(st, "beam"), s = r.beam.space;
  const header = `${continuityName(k)} splines, degree ${p}, ${s.ne} elements: ${s.n} coefficients, ` +
    `${r.beam.dofs} free after the ${st.supports} supports`;
  const below = st.motion === "response"
    ? response(bottom, st, d, r, responses(JSON.stringify([p, k, st.ne, st.supports, st.section, st.load, st.forceRatio, st.zeta]), () => {
      const bcf = beamCaseOf(st), h = harmonicOf(bcf);
      const u = evaluateQoI(r.beam, "response", loadOf(bcf), qoiPoint(st.supports), h);
      return { h, re: u.c, im: u.ci! };
    }), t)
    : modes(bottom, st, d, r, t);
  const blocks: Block[] = [ro.lead(header), ...deflection(top, st, d, r), ...below];
  if (d.dimensional) blocks.push(ro.note(referenceNote(d.ref)));
  return { readout: { sections: [{ blocks }] }, animate: true };
}

const samplesOf = (r: Solved) => Array.from({ length: r.beam.space.ne * 24 + 1 }, (_, i) => i / (r.beam.space.ne * 24));

function deflection(plot: Plot, st: State, d: Display, r: Solved): Block[] {
  const bc = { supports: st.supports, load: st.load, section: st.section };
  const { beam } = r, s = beam.space, L = d.dimensional ? d.ref.L : 1, xs = samplesOf(r);
  const xUnit = d.dimensional ? "x [m]" : "x / L";
  const exact = exactStaticOf(bc);
  const wHat = xs.map((x) => beam.evaluate(r.c, x)[0]);
  const W = d.dimensional ? deflectionScale(d.ref, st.load) : 1;
  const maxAbs = Math.max(...wHat.map(Math.abs)) * W;
  const u = d.dimensional ? axisUnit(maxAbs, "m") : { factor: 1, label: "" };
  const yv = (w: number) => w * W * u.factor;
  const series: Series[] = [{ label: "w_h (spline solution)", x: xs.map((x) => x * L), y: wHat.map(yv), color: SLOT[0], width: 2.5 }];
  if (exact) series.push({
    label: "w (exact)", x: xs.map((x) => x * L), y: xs.map((x) => yv(exact.w(x))), color: EXACT_INK, width: 1.25, dash: [6, 5],
  });
  if (st.polygon) series.push({
    label: "control polygon", x: Array.from(greville(s), (g) => g * L), y: Array.from(r.c, yv), color: SLOT[1], width: 1, markers: true,
  });
  const top = Math.max(...series.flatMap((se) => Array.from(se.y)), 1e-300);
  plot.draw({
    title: `static deflection — ${st.load === "uniform" ? "uniform load" : "point load"}, ${st.supports}`,
    xlabel: xUnit,
    ylabel: d.dimensional ? `w [${u.label}]  (downward)` : `w EI / ${st.load === "uniform" ? "q₀L⁴" : "P₀L³"}  (downward)`,
    xlim: [0, L],
    ylim: [-0.32 * top, 1.12 * top],
    yflip: true,
    series,
    hover: (se, i) => `${se.label}\n${xUnit} = ${se.x[i].toPrecision(4)}\nw = ${se.y[i].toPrecision(5)}`,
    under: (ctx, a) => baseline(ctx, a, L),
    over: (ctx, a) => {
      if (st.mesh) meshTicks(ctx, a, s.breaks.map((x) => x * L), 0, MESH_INK);
      supports(ctx, a, st.supports, L);
      loads(ctx, a, st.load, qoiPoint(st.supports) * L, L, top);
    },
  });
  const xq = qoiPoint(st.supports), w = beam.evaluate(r.c, xq)[0];
  const C = r.F.reduce((acc, f, i) => acc + f * r.c[i], 0);
  const rel = (v: number, e: number | undefined) => (e === undefined ? "—" : Math.abs(v / e - 1).toExponential(2));
  const where = st.supports === "cantilever" ? "tip" : "midspan";
  const blocks = [
    ro.table(["static", "spline", "exact", "rel. error"], [
      [`w at ${where}`, fmt.deflection(d, w, st.load), exact ? fmt.deflection(d, exact.w(xq), st.load) : "—", rel(w, exact?.w(xq))],
      ["compliance ℓ(w)", fmt.compliance(d, C, st.load), exact ? fmt.compliance(d, exact.compliance, st.load) : "—", rel(C, exact?.compliance)],
    ]),
    ro.note(`assembled, solved, and ${r.modes.values.length} modes found in ${r.ms.toFixed(1)} ms`),
  ];
  if (!exact) blocks.push(ro.note("tapered section: no closed form — compare meshes in the hierarchy view."));
  return blocks;
}

function modes(plot: Plot, st: State, d: Display, r: Solved, t: number): Block[] {
  const { beam } = r, s = beam.space, L = d.dimensional ? d.ref.L : 1, xs = samplesOf(r);
  const n = Math.min(st.mode, r.modes.values.length) - 1;
  const phi = r.modes.vectors[n];
  const shape = xs.map((x) => beam.evaluate(phi, x)[0]);
  let peak = 0;
  for (const v of shape) if (Math.abs(v) > Math.abs(peak)) peak = v;
  const amp = Math.cos((2 * Math.PI * t) / PERIOD_MS);
  const X = xs.map((x) => x * L);
  const series: Series[] = [
    { label: "envelope", x: X, y: shape.map((v) => v / peak), color: "rgba(57, 135, 229, 0.30)", width: 1, unlisted: true, inert: true },
    { label: "envelope", x: X, y: shape.map((v) => -v / peak), color: "rgba(57, 135, 229, 0.30)", width: 1, unlisted: true, inert: true },
    { label: `mode ${n + 1}`, x: X, y: shape.map((v) => (amp * v) / peak), color: SLOT[0], width: 2.5, unlisted: true },
  ];
  if (st.polygon) series.push({
    label: "control polygon", x: Array.from(greville(s), (g) => g * L), y: Array.from(phi, (v) => (amp * v) / peak),
    color: SLOT[1], width: 1, markers: true, unlisted: true,
  });
  plot.draw({
    title: `mode ${n + 1} of the ${st.supports} beam`,
    xlabel: d.dimensional ? "x [m]" : "x / L",
    ylabel: "mode shape (unit amplitude)",
    xlim: [0, L],
    ylim: [-1.3, 1.3],
    series,
    under: (ctx, a) => baseline(ctx, a, L),
    over: (ctx, a) => {
      if (st.mesh) meshTicks(ctx, a, s.breaks.map((x) => x * L), 0, MESH_INK);
      supports(ctx, a, st.supports, L);
    },
  });

  const ex = st.section === "uniform" ? exactEigenvalues(SUPPORTS[st.supports], r.modes.values.length) : null;
  const unit = d.dimensional ? "f [Hz]" : "ω̂ = ω √(ρA L⁴/EI)";
  const rows = Array.from(r.modes.values, (lam, i) => {
    const w = Math.sqrt(lam), we = ex ? Math.sqrt(ex[i]) : NaN;
    return [
      String(i + 1),
      fmt.frequency(d, w),
      ex ? fmt.frequency(d, we) : "—",
      ex ? ((w - we) / we).toExponential(2) : "—",
    ];
  });
  const blocks = [ro.table(["mode", "spline", "exact", "rel. error"], rows, { caption: `lowest ${rows.length} natural frequencies, ${unit}`, mark: n })];
  if (ex) blocks.push(ro.note("every spline frequency lies above the exact one — Rayleigh–Ritz in a conforming space only overestimates."));
  else blocks.push(ro.note("tapered section: no closed form for its frequencies."));
  if (!d.dimensional) blocks.push(ro.note(`λ̂ = ω̂² = ${plain(r.modes.values[n])} for the mode shown (highlighted)`));
  return blocks;
}

function response(plot: Plot, st: State, d: Display, r: Solved, u: Response, t: number): Block[] {
  const { beam } = r, s = beam.space, L = d.dimensional ? d.ref.L : 1, xs = samplesOf(r);
  const re = xs.map((x) => beam.evaluate(u.re, x)[0]), im = xs.map((x) => beam.evaluate(u.im, x)[0]);
  const amp = re.map((v, i) => Math.hypot(v, im[i])), w0 = xs.map((x) => beam.evaluate(r.c, x)[0]);
  const W = d.dimensional ? deflectionScale(d.ref, st.load) : 1;
  const top = Math.max(...amp, ...w0.map(Math.abs)) * W;
  const un = d.dimensional ? axisUnit(top, "m") : { factor: 1, label: "" };
  const k = W * un.factor, th = (2 * Math.PI * t) / PERIOD_MS, c = Math.cos(th), sn = Math.sin(th);
  const X = xs.map((x) => x * L);
  plot.draw({
    title: `forced response at Ω = ${st.forceRatio} ω₁ (uniform beam), ζ = ${st.zeta} — Re(u e^{iΩt}), swinging`,
    xlabel: d.dimensional ? "x [m]" : "x / L",
    ylabel: d.dimensional ? `w [${un.label}]  (downward)` : `w EI / ${st.load === "uniform" ? "q₀L⁴" : "P₀L³"}  (downward)`,
    xlim: [0, L], ylim: [-1.25 * top, 1.25 * top], yflip: true,
    series: [
      { label: "envelope ±|u|", x: X, y: amp.map((v) => v * k), color: "rgba(57, 135, 229, 0.35)", width: 1, inert: true },
      { label: "envelope", x: X, y: amp.map((v) => -v * k), color: "rgba(57, 135, 229, 0.35)", width: 1, inert: true, unlisted: true },
      { label: "static deflection, for scale", x: X, y: w0.map((v) => v * k), color: EXACT_INK, width: 1.25, dash: [6, 5], inert: true },
      { label: "Re(u e^{iΩt})", x: X, y: re.map((v, i) => (v * c - im[i] * sn) * k), color: SLOT[0], width: 2.5 },
    ],
    legend: "tl",
    under: (ctx, a) => baseline(ctx, a, L),
    over: (ctx, a) => {
      if (st.mesh) meshTicks(ctx, a, s.breaks.map((x) => x * L), 0, MESH_INK);
      supports(ctx, a, st.supports, L);
    },
    hover: (se, i) => `${se.label}\nx = ${se.x[i].toPrecision(3)}\n|u| = ${(amp[i] * k).toPrecision(4)}`,
  });

  const xq = qoiPoint(st.supports), where = st.supports === "cantilever" ? "tip" : "midspan";
  const uq = [beam.evaluate(u.re, xq)[0], beam.evaluate(u.im, xq)[0]], aq = Math.hypot(uq[0], uq[1]);
  const wq = beam.evaluate(r.c, xq)[0];
  const exact = st.section === "uniform" && st.supports === "pinned–pinned" ? pinnedResponse(st.load, u.h) : null;
  // As a lag in [0°, 360°).
  const lag = ((((Math.atan2(-uq[1], uq[0]) * 180) / Math.PI) % 360) + 360) % 360;
  return [
    ro.lead(`forced at Ω = ${plain(u.h.Omega, 5)} (${st.forceRatio} × ω₁ of the uniform beam), ζ = ${st.zeta}`),
    ro.table(["at the " + where, "spline", exact === null ? "exact" : "modal series", "rel. error"], [
      ["amplitude |u|", fmt.deflection(d, aq, st.load), exact === null ? "—" : fmt.deflection(d, exact, st.load),
        exact === null ? "—" : Math.abs(aq / exact - 1).toExponential(2)],
      ["static w", fmt.deflection(d, wq, st.load), "", ""],
    ]),
    ro.tiles([
      { label: "dynamic amplification", value: plain(aq / Math.abs(wq), 4), detail: "|u| / w_static" },
      { label: "phase lag behind the load", value: `${plain(lag, 3)}°` },
    ], 2),
    ro.note(st.forceRatio < 1
      ? "below the first resonance the beam moves with the load, nearly in phase; raise Ω/ω₁ toward 1 and the amplitude climbs, limited only by ζ."
      : "above the first resonance the first mode's share swings into antiphase with the load, and the modes nearest Ω take over the shape."),
    ro.note(`Rayleigh damping C = aM + bK with a = ${plain(u.h.a, 4)}, b = ${plain(u.h.b, 4)}: ζ = ${st.zeta} on the uniform beam's first two modes`),
  ];
}

function baseline(ctx: CanvasRenderingContext2D, a: Axes, L: number): void {
  ctx.strokeStyle = "rgba(207, 238, 255, 0.35)";
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 4]);
  ctx.beginPath();
  ctx.moveTo(a.sx(0), a.sy(0));
  ctx.lineTo(a.sx(L), a.sy(0));
  ctx.stroke();
  ctx.setLineDash([]);
}

/** Wall with hatching for a clamp, a triangle for a pin, nothing for a free end. */
function supports(ctx: CanvasRenderingContext2D, a: Axes, name: keyof typeof SUPPORTS, L: number): void {
  const { left, right } = SUPPORTS[name];
  const draw = (end: End, x: number, side: -1 | 1) => {
    const X = a.sx(x), Y = a.sy(0);
    ctx.strokeStyle = "rgba(207, 238, 255, 0.85)";
    ctx.fillStyle = "rgba(207, 238, 255, 0.85)";
    ctx.lineWidth = 2;
    if (end === "clamped") {
      ctx.beginPath(); ctx.moveTo(X, Y - 22); ctx.lineTo(X, Y + 22); ctx.stroke();
      ctx.lineWidth = 1;
      for (let i = -22; i < 22; i += 7) {
        ctx.beginPath(); ctx.moveTo(X, Y + i + 7); ctx.lineTo(X + 7 * side, Y + i); ctx.stroke();
      }
    } else if (end === "pinned") {
      ctx.beginPath(); ctx.moveTo(X, Y); ctx.lineTo(X - 8, Y + 14); ctx.lineTo(X + 8, Y + 14); ctx.closePath(); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(X - 12, Y + 18); ctx.lineTo(X + 12, Y + 18); ctx.stroke();
    }
  };
  draw(left, 0, -1);
  draw(right, L, 1);
}

/** Load arrows above the beam, pointing the way w grows. */
function loads(ctx: CanvasRenderingContext2D, a: Axes, kind: "uniform" | "point", xp: number, L: number, top: number): void {
  ctx.strokeStyle = "#c98500";
  ctx.fillStyle = "#c98500";
  ctx.lineWidth = 1.5;
  const arrow = (x: number, len: number) => {
    const X = a.sx(x), Y1 = a.sy(-0.06 * top), Y0 = Y1 - len;
    ctx.beginPath(); ctx.moveTo(X, Y0); ctx.lineTo(X, Y1 - 5); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(X, Y1); ctx.lineTo(X - 4, Y1 - 7); ctx.lineTo(X + 4, Y1 - 7); ctx.closePath(); ctx.fill();
  };
  if (kind === "uniform") {
    for (let i = 0; i <= 16; i++) arrow((i / 16) * L, 16);
    const Y = a.sy(-0.06 * top) - 16;
    ctx.beginPath(); ctx.moveTo(a.sx(0), Y); ctx.lineTo(a.sx(L), Y); ctx.stroke();
  } else arrow(xp, 34);
}
