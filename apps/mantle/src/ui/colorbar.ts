// This file is part of the Mantle app, a WebGPU mantle convection simulator.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * A colour-bar legend: a horizontal gradient sampling the active colour map,
 * labelled with its endpoints. Every quantity this draws for is mapped onto
 * [0, 1] before it is ever coloured, so the endpoints are usually fixed. The
 * one exception is temperature under internal heating: the field's map then
 * spans [0, max T] (`tMaxSource` in `gpu/wgsl.ts`), and `setLabels` is how
 * the main legend follows it.
 *
 * `labels` defaults to the bare `["0", "1"]` because this one block also
 * legends the particle-tint colour bar (`pcbar` in controls.ts), whose
 * quantity changes with `PARTICLE_TINT`'s own selection — speed, age,
 * initial φ, species — and only some of those endpoints are "cold"/"hot".
 * The main temperature legend (`cbar`) is the one call that knows its
 * quantity never changes, so it is the one call that passes its own pair.
 *
 * A CSS gradient rather than a canvas: the same five control points already
 * drive the shader (`colormaps.ts`), and a `linear-gradient` reproduces their
 * piecewise-linear interpolation exactly, without a redraw call of its own
 * when the pane resizes.
 */

import { COLORMAPS, type ColormapName } from "../colormaps";

const rgb = ([r, g, b]: readonly [number, number, number]): string =>
  `rgb(${Math.round(r * 255)} ${Math.round(g * 255)} ${Math.round(b * 255)})`;

export function colorbarBlock(
  initial: ColormapName,
  labels: readonly [string, string] = ["0", "1"],
): {
  el: HTMLElement;
  setColormap: (name: ColormapName) => void;
  setLabels: (labels: readonly [string, string]) => void;
} {
  const el = document.createElement("div");
  el.className = "cbar";

  const bar = document.createElement("div");
  bar.className = "cbar-bar";

  const ticks = document.createElement("div");
  ticks.className = "cbar-ticks";
  const lo = document.createElement("span");
  lo.textContent = labels[0];
  const hi = document.createElement("span");
  hi.textContent = labels[1];
  ticks.append(lo, hi);

  el.append(bar, ticks);

  const setColormap = (name: ColormapName): void => {
    const stops = COLORMAPS[name];
    const step = 100 / (stops.length - 1);
    const css = stops
      .map((s, i) => `${rgb(s)} ${(i * step).toFixed(3)}%`)
      .join(", ");
    bar.style.background = `linear-gradient(to right, ${css})`;
  };
  setColormap(initial);

  // Compared first: the main legend is offered a label every frame, and an
  // unchanged string should not touch the DOM.
  const setLabels = ([a, b]: readonly [string, string]): void => {
    if (lo.textContent !== a) lo.textContent = a;
    if (hi.textContent !== b) hi.textContent = b;
  };

  return { el, setColormap, setLabels };
}
