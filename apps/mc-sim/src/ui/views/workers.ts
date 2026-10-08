// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The one worker pool the sampling views share, and the pane's handle on it.
 * Created on first use, so the views that never sample never start a worker.
 */

import { McPool } from "../../mc/pool";

let pool: McPool | null = null;
let notify: () => void = () => {};

export const workers = {
  /** Give the pool a way to ask for a redraw. */
  connect(request: () => void): void { notify = request; },
  get pool(): McPool { return (pool ??= new McPool(() => notify())); },
  get size(): number { return pool?.size ?? 0; },
  /** Samples the run on screen has folded in, over all its levels; 0 before any pool exists. */
  samples(): number { return pool?.current()?.streams.reduce((n, s) => n + s.acc.n, 0) ?? 0; },
  /** Pause or resume the run on screen. */
  toggle(): void {
    const r = pool?.current();
    if (!r) return;
    if (r.paused) pool!.resume(); else pool!.pause();
    notify();
  },
  /** The sampling view left the screen: stop dispatching, keep every run. */
  idle(): void { pool?.idle(); },
};
