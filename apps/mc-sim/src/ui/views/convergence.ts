// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 3 on screen: the hierarchy of discretisations, as error against
 * element size on log–log axes, one series per degree — for the beam on the
 * left and the plate on the right, whatever the pane's structure.
 *
 * Solid: the true error |Q_ℓ − Q|, where a closed form gives Q. Dashed: the
 * successive difference |Q_ℓ − Q_{ℓ−1}|, which needs no Q at all and is the
 * quantity multilevel Monte Carlo is built from. Both are relative to |Q|, so
 * they read the same in either unit system. The degree selected in the pane is
 * drawn at full weight and tabulated below; the others stay for comparison.
 *
 * The grey line is that degree's round-off, measured per level (see
 * `Level.noise`): it rises as h falls, because the beam's stiffness is
 * conditioned like h⁻⁴. Where an error curve meets it, refining further buys
 * nothing — the finest level worth paying for, which MLMC will need to know.
 */

import { CLEAR, runHierarchy, theoryRate, type Hierarchy, type QoI } from "../../hierarchy";
import { plateTheoryRate, runPlateHierarchy } from "../../plate/qoi";
import type { Figure } from "../figure";
import { SLOT, decade, type Plot, type Series } from "../plot";
import { beamCaseOf, continuityK, continuityName, displayOf, plateCaseOf, type State } from "../state";
import { fmt, plain, referenceNote, si, type Display } from "../units";
import { PLATE_MAX_NE, admissibleAt, isPlate, levelsWithin, structureText } from "./structure";
import { ro, section, type Block } from "../readout";
import { memo, type ViewResult } from "./view";

/** The degrees drawn: the quadratic and the cubic, the two a beam is usually built from. */
const DEGREES = [2, 3];
const NEUTRAL = "rgba(207, 238, 255, 0.80)";

export const QOI_NAME: Record<QoI, string> = {
  omega1: "first natural frequency ω₁",
  deflection: "deflection at the tip / midspan / centre",
  compliance: "compliance ℓ(w)",
  field: "RMS deflection ‖w‖ (a field)",
  response: "forced response amplitude |w| at Ω",
};

const runs = memo<Map<number, Hierarchy | string>>();

export const valueOf = (d: Display, qoi: QoI, v: number, load: State["load"]) =>
  qoi === "omega1" ? fmt.frequency(d, v)
  : qoi === "compliance" ? fmt.compliance(d, v, load)
  : fmt.deflection(d, v, load);

/** Both hierarchies side by side, the beam's readout above the plate's. */
export function renderConvergence(fig: Figure, st: State): ViewResult {
  const [left, right] = fig.panels(2, { cols: 2 });
  const beam = renderHierarchy(left, { ...st, structure: "beam" });
  const plate = renderHierarchy(right, { ...st, structure: "plate" });
  return { readout: { sections: [section("BEAM", beam.readout), section("PLATE", plate.readout)] }, animate: false };
}

function renderHierarchy(plot: Plot, st: State): ViewResult {
  const bc = beamCaseOf(st), pc = plateCaseOf(st), plate = isPlate(st);
  // A plate's hierarchy runs on this thread, four degrees at once, and stops at PLATE_MAX_NE.
  const levels = plate ? levelsWithin(st.ne0, st.levels, PLATE_MAX_NE) : st.levels;
  const all = runs(JSON.stringify([st.continuity, plate ? pc : bc, st.qoi, st.ne0, levels]), () =>
    new Map(DEGREES.map((p): [number, Hierarchy | string] => {
      const k = continuityK(st.continuity, p);
      const why = admissibleAt(st, st.ne0, p);
      if (why) return [p, why];
      const spec = { p, k, qoi: st.qoi, ne0: st.ne0, levels };
      return [p, plate ? runPlateHierarchy({ ...spec, plate: pc }) : runHierarchy({ ...spec, beam: bc })];
    })));
  const theory = (p: number) => (plate ? plateTheoryRate(st.qoi, p, pc) : theoryRate(st.qoi, p, bc));

  if (st.continuity === "c0") {
    plot.draw({ xlabel: "element size h / L", ylabel: "relative error", xlog: true, ylog: true, series: [] });
    return { readout: `${all.get(3) as string}\nChoose C¹ or maximal continuity to build a hierarchy.`, animate: false };
  }

  const d = displayOf(st), L = d.dimensional ? d.ref.L : 1;
  const series: Series[] = [];
  for (const p of DEGREES) {
    const H = all.get(p)!;
    if (typeof H === "string") continue;
    const scale = Math.abs(H.exact ?? H.levels[H.levels.length - 1].Q) || 1;
    const focus = p === st.p, color = SLOT[p - 2];
    const h = H.levels.map((l) => l.h * L);
    const th = theory(p);
    const alpha = H.alphaExact ?? H.alpha;
    const label = `p = ${p}   α ≈ ${alpha === null ? "—" : alpha.toFixed(2)}${th === null ? "" : `  (theory ${th})`}`;
    const base = { color, alpha: focus ? 1 : 0.55, width: focus ? 2.5 : 1.5, markers: true };
    if (H.exact !== null) {
      series.push({ ...base, label, x: h, y: H.levels.map((l) => Math.max(l.err / scale, 1e-17)) });
      series.push({ ...base, label: `p = ${p} successive differences`, x: h, y: H.levels.map((l) => l.dQ / scale), dash: [6, 4], unlisted: true });
    } else {
      series.push({ ...base, label, x: h, y: H.levels.map((l) => l.dQ / scale), dash: [6, 4] });
    }
  }
  const focusH = all.get(st.p);
  if (focusH && typeof focusH !== "string") {
    const scale = Math.abs(focusH.exact ?? focusH.levels[focusH.levels.length - 1].Q) || 1;
    series.push({
      label: `round-off, p = ${st.p} (measured)`, x: focusH.levels.map((l) => l.h * L),
      y: focusH.levels.map((l) => l.noise / scale), color: "rgba(207, 238, 255, 0.45)", width: 1.25, dash: [2, 3],
    });
  }
  // Legend keys for the two line styles.
  const anyExact = [...all.values()].some((H) => typeof H !== "string" && H.exact !== null);
  if (anyExact) series.push({ label: "— true error |Q_ℓ − Q|/|Q|", x: [], y: [], color: NEUTRAL, width: 1.5, inert: true });
  series.push({ label: "- - successive |Q_ℓ − Q_ℓ₋₁|/|Q|", x: [], y: [], color: NEUTRAL, width: 1.5, dash: [6, 4], inert: true });

  const ys = series.flatMap((s) => Array.from(s.y)).filter((v) => v > 0 && Number.isFinite(v));
  const ylo = Math.max(1e-17, Math.min(...ys) / 3), top = Math.max(...ys) * 3;
  // Headroom for the legend, top left: about a quarter of the axis above the data.
  const yhi = top * (top / ylo) ** 0.3;
  const H0 = all.get(DEGREES[0])!;
  const hs = typeof H0 === "string" ? [1] : H0.levels.map((l) => l.h * L);
  plot.draw({
    title: `${QOI_NAME[st.qoi]} — ${st.continuity === "max" ? "maximal continuity Cᵖ⁻¹" : "C¹"}, ${structureText(st)}`,
    xlabel: d.dimensional ? "element size h [m]" : "element size h / L",
    ylabel: "relative error",
    xlog: true,
    ylog: true,
    // Errors fall toward the bottom-left and round-off rises along the bottom:
    // the top-left corner, above the coarse meshes' errors, stays empty.
    legend: "tl",
    xlim: [Math.min(...hs) / 1.5, Math.max(...hs) * 1.5],
    ylim: [ylo, yhi],
    xfmt: (v) => (d.dimensional ? si(v, "m", 2) : decade(v)),
    series,
    hover: (s, i) => {
      const round = s.label.startsWith("round-off");
      const p = round ? st.p : DEGREES.find((q) => s.label.startsWith(`p = ${q}`))!;
      const H = all.get(p) as Hierarchy, l = H.levels[i];
      const kind = round ? "round-off estimate" : s.dash ? "successive difference" : "true error";
      return `p = ${p}, ${kind}\nlevel ${l.level}: ${l.ne} elements, ${l.dofs} dofs\nrelative ${s.y[i].toExponential(2)}`;
    },
  });

  // The table: the degree selected in the pane, if it is one of those plotted.
  const p = DEGREES.includes(st.p) ? st.p : 3;
  const H = all.get(p)!;
  if (typeof H === "string") return { readout: `p = ${p}: ${H}`, animate: false };
  const scale = Math.abs(H.exact ?? H.levels[H.levels.length - 1].Q) || 1;
  const rows = H.levels.map((l) => [
    String(l.level), String(l.ne), String(l.dofs), valueOf(d, st.qoi, l.Q, st.load),
    Number.isNaN(l.dQ) ? "—" : (l.dQ / scale).toExponential(2),
    Number.isNaN(l.err) ? "—" : (l.err / scale).toExponential(2),
    (l.noise / scale).toExponential(1),
    l.ms.toFixed(2),
  ]);
  const k = continuityK(st.continuity, p), th = theory(p);
  const blocks: Block[] = [
    ro.lead(`p = ${p}, ${continuityName(k)}: levels ℓ = 0…${levels - 1}, ne = ${st.ne0}·2^ℓ${plate ? " per side" : ""}`),
    ro.table(["ℓ", "ne", "dofs", "Q_ℓ", "|ΔQ_ℓ|/|Q|", "|Q_ℓ−Q|/|Q|", "round-off", "ms"], rows),
    ro.tiles([
      {
        label: "exact Q", value: H.exact === null ? "—" : valueOf(d, st.qoi, H.exact, st.load),
        detail: H.exact === null ? "no closed form: only successive differences" : undefined,
      },
      {
        label: "α, convergence rate", value: H.alpha === null ? "—" : plain(H.alpha, 3),
        detail: [H.alphaExact === null ? "" : `true error ${plain(H.alphaExact, 3)}`, th === null ? "" : `theory ${th}`].filter((x) => x).join(" · ") || "from successive differences",
        help: `|Q_ℓ − Q| ~ h^α, fitted to the finest three levels standing ${CLEAR}× clear of round-off: from the successive differences |ΔQ_ℓ|, and from the true error where the exact Q is known.`,
      },
      {
        label: plate ? "γ, work ~ h^−γ" : "γ, dofs ~ h^−γ", value: H.gamma === null ? "—" : plain(H.gamma, 3),
        detail: plate ? "banded solve: dofs × bandwidth² → γ = 4" : "banded solve O(dofs·p²) → γ = 1",
        help: plate
          ? "A banded plate solve costs dofs × bandwidth² ~ h⁻² · h⁻², so γ → 4 in 2D."
          : "Banded solves cost O(dofs·p²), so the work grows like the dofs: γ = 1 in 1D.",
      },
    ]),
  ];
  if (levels < st.levels)
    blocks.push(ro.note(`${st.levels - levels} level${st.levels - levels > 1 ? "s" : ""} left out: this view stops a plate at ${PLATE_MAX_NE} × ${PLATE_MAX_NE}`));
  if (p !== st.p) blocks.push(ro.note(`degree ${st.p} is not drawn here: this view shows p = ${DEGREES.join(" and ")}`));
  if (st.load === "point" && st.qoi !== "omega1")
    blocks.push(ro.warn(plate
      ? "point load: w ~ r² log r under it, so no smooth-data rate holds."
      : "point load: w‴ jumps under it, so the smooth-data rates need not hold — and with the load on a knot the exact solution may lie in the space."));
  if (plate && st.edges !== "SSSS")
    blocks.push(ro.warn(`${st.edges}: where a clamped or free edge meets another, the solution carries a corner singularity r^s that can cap α below the smooth rate — ` +
      "the clamped–free corners of CFFF hold it near 2 at every degree."));
  if (d.dimensional) blocks.push(ro.note(referenceNote(d.ref)));
  return { readout: { sections: [{ blocks }] }, animate: false };
}
