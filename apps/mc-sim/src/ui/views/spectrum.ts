// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The whole discrete spectrum against the exact one: ω_h,n / ω_n for every
 * mode n the mesh carries, plotted against n/N.
 *
 * This is where continuity earns its keep. On the same mesh at the same degree,
 * maximal-continuity splines resolve most of their spectrum to within a few per
 * cent and misbehave only in a handful of "outlier" modes at the very top,
 * while C¹ (Hermite-like) elements split into branches and lose accuracy over
 * a large fraction of it (Cottrell, Reali, Bazilevs & Hughes, CMAME 2006;
 * Hughes, Evans & Reali, CMAME 2014). Every ratio is ≥ 1: Rayleigh–Ritz.
 *
 * A dense eigensolve per curve, so the mesh is capped for this view.
 */

import { Beam, SUPPORTS, admissible } from "../../beam/beam";
import { betaL } from "../../beam/exact";
import { SLOT, type Plot, type Series } from "../plot";
import { continuityName, type State } from "../state";
import { memo, table, type ViewResult } from "./view";

export const SPECTRUM_MAX_NE = 48;
const YMAX = 2;

interface Curve { k: number; ratio: number[]; N: number; }
const curves = memo<Curve[] | string>();

export function renderSpectrum(plot: Plot, st: State): ViewResult {
  const p = st.p, ne = Math.min(st.ne, SPECTRUM_MAX_NE);
  const supports = SUPPORTS[st.supports];
  const ks = [...new Set([p - 1, 1])].filter((k) => k >= 1);
  const result = curves(JSON.stringify([p, ne, st.supports, st.section]), () => {
    if (st.section !== "uniform") return "The exact spectrum is known in closed form only for a uniform section.";
    const why = admissible({ p, k: 1, ne, supports });
    if (why) return why;
    return ks.map((k) => {
      const lam = new Beam({ p, k, ne, supports }).spectrum();
      const ratio = Array.from(lam, (l, i) => Math.sqrt(l) / betaL(supports, i + 1)! ** 2);
      return { k, ratio, N: lam.length };
    });
  });
  if (typeof result === "string") {
    plot.draw({ xlabel: "n / N", ylabel: "ω_h / ω", xlim: [0, 1], ylim: [0.95, YMAX], series: [] });
    return { readout: result, animate: false };
  }

  const series: Series[] = result.map((c, j) => ({
    label: `${continuityName(c.k)}${c.k === p - 1 ? " (maximal)" : ""}: N = ${c.N}`,
    x: c.ratio.map((_, i) => (i + 1) / c.N),
    y: c.ratio.map((r) => Math.min(r, YMAX * 1.5)),
    color: SLOT[j],
    width: 2,
  }));
  plot.draw({
    title: `normalised spectrum, degree ${p}, ${ne} elements, ${st.supports}`,
    xlabel: "mode number / number of modes   n / N",
    ylabel: "ω_h,n / ω_n",
    xlim: [0, 1],
    ylim: [0.98, YMAX],
    series,
    hover: (s, i) => {
      const c = result[series.indexOf(s)];
      return `${continuityName(c.k)}, mode ${i + 1} of ${c.N}\nω_h/ω = ${c.ratio[i].toFixed(5)}`;
    },
  });

  const rows = result.map((c) => {
    const within = (tol: number) => c.ratio.filter((r) => r - 1 <= tol).length;
    const above = c.ratio.filter((r) => r > YMAX).length;
    return [
      continuityName(c.k), String(c.N),
      `${((100 * within(0.01)) / c.N).toFixed(0)}%`,
      `${((100 * within(0.1)) / c.N).toFixed(0)}%`,
      Math.max(...c.ratio).toFixed(2), String(above),
    ];
  });
  const lines = [
    table(["space", "modes", "≤ 1% err", "≤ 10% err", "max ω_h/ω", `off-scale (> ${YMAX})`], rows),
    "",
    "Same mesh, same degree: maximal continuity has fewer coefficients and resolves more of its spectrum.",
  ];
  if (st.ne > SPECTRUM_MAX_NE) lines.push(`(mesh capped at ${SPECTRUM_MAX_NE} elements here: every eigenvalue is a dense solve)`);
  if (p === 2) lines.push("At p = 2, C¹ is already maximal continuity: one curve.");
  return { readout: lines.join("\n"), animate: false };
}
