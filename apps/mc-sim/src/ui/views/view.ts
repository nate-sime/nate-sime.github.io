// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/** What a view hands back after drawing: the text under the plot, and whether it wants frames. */
export interface ViewResult {
  readonly readout: string;
  readonly animate: boolean;
}

/** Memoise the last result per key — the pane redraws far more often than inputs change. */
export function memo<T>(): (key: string, make: () => T) => T {
  let lastKey: string | null = null, last: T;
  return (key, make) => {
    if (key !== lastKey) { last = make(); lastKey = key; }
    return last;
  };
}

/** A fixed-width table, right-aligned, for the readout. */
export function table(head: string[], rows: string[][]): string {
  const w = head.map((h, j) => Math.max(h.length, ...rows.map((r) => r[j].length)));
  const line = (r: string[]) => r.map((c, j) => c.padStart(w[j])).join("  ");
  return [line(head), ...rows.map(line)].join("\n");
}
