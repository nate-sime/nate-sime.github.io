// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Stage 1 on screen: every basis function of the chosen space (or its first or
 * second derivative), with the knots and their multiplicities under the axis.
 *
 * Each element is drawn on its own, one-sided at both ends, and the pen lifts
 * between elements. A derivative that jumps at a knot therefore shows as two
 * ends that do not meet — which is what C^k means — rather than as a vertical
 * segment the space does not contain. All functions share one hue: what is on
 * show is their shape and overlap, not which is which (hover names one).
 */

import { basisOnElement, firstDof, splineSpace } from "../../spline";
import { SLOT, type Plot, type Series } from "../plot";
import { continuityK, continuityName, type State } from "../state";
import type { ViewResult } from "./view";

const PER = 32;

export function renderBasis(plot: Plot, st: State): ViewResult {
  const p = Math.max(1, st.p), k = continuityK(st.continuity, p), r = Math.min(st.derivative, p);
  const s = splineSpace(p, st.ne, k);
  const P1 = p + 1;
  const xs: number[][] = Array.from({ length: s.n }, () => []);
  const ys: number[][] = Array.from({ length: s.n }, () => []);
  const sumX: number[] = [], sumY: number[] = [];
  for (let e = 0; e < s.ne; e++) {
    const a = s.breaks[e], b = s.breaks[e + 1], f0 = firstDof(s, e);
    for (let j = 0; j < P1; j++) {
      // A gap from the previous element this function lived on.
      if (xs[f0 + j].length) { xs[f0 + j].push(NaN); ys[f0 + j].push(NaN); }
    }
    sumX.push(NaN); sumY.push(NaN);
    for (let i = 0; i <= PER; i++) {
      const x = a + ((b - a) * i) / PER, N = basisOnElement(s, e, x, r);
      let sum = 0;
      for (let j = 0; j < P1; j++) {
        xs[f0 + j].push(x);
        ys[f0 + j].push(N[r * P1 + j]);
        sum += N[r * P1 + j];
      }
      sumX.push(x); sumY.push(sum);
    }
  }
  const prime = ["", "′", "″"][r];
  const series: Series[] = xs.map((x, i) => ({
    label: `B${sub(i)}${prime}`, x, y: ys[i], color: SLOT[0], width: 1.75, unlisted: true,
  }));
  series.push({
    label: r === 0 ? "Σ Bᵢ = 1" : `Σ Bᵢ${prime} = 0`, x: sumX, y: sumY, color: "rgba(207, 238, 255, 0.55)",
    width: 1.25, dash: [5, 4], unlisted: true, inert: true,
  });

  plot.draw({
    title: `${continuityName(k)} splines of degree ${p} — ${s.n} functions on ${s.ne} element${s.ne > 1 ? "s" : ""}`,
    xlabel: "x / L",
    ylabel: `Bᵢ${prime}(x)`,
    xlim: [0, 1],
    ylim: r === 0 ? [-0.06, 1.1] : undefined,
    series,
    hover: (se, i) => `${se.label}\nx/L = ${se.x[i].toFixed(4)}\nvalue = ${se.y[i].toPrecision(5)}`,
    over: (ctx, a) => {
      // Knot multiplicities, along the foot of the plot.
      ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "bottom";
      for (let i = 0; i <= s.ne; i++) {
        const X = a.sx(s.breaks[i]), mult = i === 0 || i === s.ne ? p + 1 : s.m;
        ctx.fillStyle = "#c98500";
        ctx.beginPath();
        ctx.arc(X, a.box.b, 3.5, 0, 2 * Math.PI);
        ctx.fill();
        if (s.ne <= 24 || i === 0 || i === s.ne) {
          ctx.fillStyle = "rgba(207, 238, 255, 0.70)";
          ctx.fillText(`×${mult}`, X, a.box.b - 6);
        }
      }
    },
  });

  const knots = Array.from(s.U, (u) => +u.toFixed(4));
  const shown = knots.length > 40 ? `${knots.slice(0, 18).join(" ")} … ${knots.slice(-18).join(" ")}` : knots.join(" ");
  const lines = [
    `degree p = ${p}, continuity ${continuityName(k)} across interior knots`,
    `interior knot multiplicity  m = p − k = ${s.m}`,
    `basis functions             n = p + 1 + (ne − 1)·m = ${p + 1} + ${s.ne - 1}·${s.m} = ${s.n}`,
    `nonzero on each element     p + 1 = ${P1}`,
    `knot vector  [${shown}]`,
  ];
  if (st.continuity === "c1" && p < 2) lines.push("(degree 1 cannot be C¹: shown at C⁰)");
  if (st.derivative > p) lines.push(`(degree ${p} has no derivative of order ${st.derivative}: shown at ${p})`);
  if (k === 0) lines.push("C⁰: fine for a string or a bar — not for a beam, whose energy needs w″ (see the beam view).");
  return { readout: lines.join("\n"), animate: false };
}

const SUB = "₀₁₂₃₄₅₆₇₈₉";
const sub = (i: number) => String(i).split("").map((d) => SUB[+d]).join("");
