// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

import type { Rich } from "../readout";

/** What a view hands back after drawing: the text under the plot, and whether it wants frames. */
export interface ViewResult {
  /** Plain lines, or sections of tiles, tables and notes (`ui/readout.ts`). */
  readonly readout: string | Rich;
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
