/**
 * The temperature colour-bar legend: a horizontal gradient sampling the
 * active colour map, labelled 0 and 1 — the field's fixed range, since T is
 * clamped to [0, 1] before it is ever coloured (see `sample_T` in
 * `gpu/wgsl.ts`), so the bar's endpoints never need to track anything.
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
): { el: HTMLElement; setColormap: (name: ColormapName) => void } {
  const el = document.createElement("div");
  el.className = "cbar";

  const bar = document.createElement("div");
  bar.className = "cbar-bar";

  const ticks = document.createElement("div");
  ticks.className = "cbar-ticks";
  const lo = document.createElement("span");
  lo.textContent = "0";
  const hi = document.createElement("span");
  hi.textContent = "1";
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

  return { el, setColormap };
}
