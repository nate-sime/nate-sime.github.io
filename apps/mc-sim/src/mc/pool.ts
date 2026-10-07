// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The pool of Monte Carlo workers, and the one run it is working on.
 *
 * Samples are independent, so the parallelism is across them: the pool hands
 * each idle worker the next range of indices and folds what comes back into an
 * `Accumulator`, which restores sample order. Batch sizes adapt per worker
 * toward ~50 ms of work — small at first, so the first statistics appear at
 * once, larger once the per-sample cost is known, so message overhead stays
 * negligible.
 *
 * A new job (any input changed) bumps the job id; replies to the old one are
 * ignored, and the workers are kept. Raising the target extends a run in place.
 */

import type { KLData } from "../random/kl";
import { FIELD_POINTS, type McSpec } from "./sampler";
import { Accumulator } from "./stats";
import type { JobMessage, Reply, RunMessage } from "./worker";

export interface McJob {
  readonly spec: McSpec;
  readonly kl: KLData;
  /** Elements on the level sampled; its parent has ne/2. */
  readonly ne: number;
  readonly coarse: boolean;
  readonly wantW: boolean;
}

const BATCH_MS = 50;

interface Slot {
  readonly worker: Worker;
  busy: boolean;
  size: number;
}

export class McRun {
  readonly acc = new Accumulator(FIELD_POINTS);
  /** Samples handed out so far. */
  dispatched = 0;
  running = false;
  error: string | null = null;
  /** Wall time spent running, ms (pauses excluded). */
  private wall = 0;
  private since = 0;

  constructor(readonly id: number, readonly key: string, readonly job: McJob, public target: number) {}

  get wallMs(): number {
    return this.wall + (this.running ? performance.now() - this.since : 0);
  }

  get done(): boolean {
    return this.acc.n >= this.target;
  }

  /** @internal */ setRunning(on: boolean): void {
    if (on === this.running) return;
    if (on) this.since = performance.now();
    else this.wall += performance.now() - this.since;
    this.running = on;
  }
}

export class McPool {
  readonly size: number;
  private readonly slots: Slot[];
  private run: McRun | null = null;
  private nextId = 1;

  constructor(private readonly onUpdate: () => void, size?: number) {
    this.size = size ?? Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 4) - 1));
    this.slots = Array.from({ length: this.size }, () => {
      const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
      const slot: Slot = { worker, busy: false, size: 4 };
      worker.onmessage = (e: MessageEvent<Reply>) => this.receive(slot, e.data);
      worker.onerror = (e) => this.fail(e.message || "worker failed to load");
      return slot;
    });
  }

  /**
   * The run for `key`: the current one if the key matches (its target raised
   * if `target` is larger), else a fresh run of `job`. It starts running unless
   * it is already complete or was paused.
   */
  ensure(key: string, job: () => McJob, target: number): McRun {
    let r = this.run;
    if (!r || r.key !== key) {
      if (r) r.setRunning(false);
      r = this.run = new McRun(this.nextId++, key, job(), target);
      const m: JobMessage = { type: "job", id: r.id, ...r.job };
      for (const s of this.slots) { s.busy = false; s.size = 4; s.worker.postMessage(m); }
      r.setRunning(true);
    } else if (target !== r.target) {
      const grew = target > r.target;
      r.target = target;
      if (grew && !r.error) r.setRunning(true);
    }
    this.pump();
    return r;
  }

  current(): McRun | null {
    return this.run;
  }

  pause(): void {
    this.run?.setRunning(false);
  }

  resume(): void {
    const r = this.run;
    if (r && !r.done && !r.error) { r.setRunning(true); this.pump(); }
  }

  private pump(): void {
    const r = this.run;
    if (!r || !r.running) return;
    for (const s of this.slots) {
      if (s.busy || r.dispatched >= r.target) continue;
      const count = Math.min(s.size, r.target - r.dispatched);
      const m: RunMessage = { type: "run", id: r.id, from: r.dispatched, count };
      r.dispatched += count;
      s.busy = true;
      s.worker.postMessage(m);
    }
  }

  private receive(slot: Slot, m: Reply): void {
    const r = this.run;
    if (!r || m.id !== r.id) return; // a superseded job
    slot.busy = false;
    if (m.type === "error") { this.fail(m.message); return; }
    const b = m.batch;
    const perSample = b.ms / b.count;
    slot.size = Math.max(1, Math.min(4096, Math.round(BATCH_MS / Math.max(perSample, 1e-3)), 4 * slot.size));
    r.acc.push(b);
    if (r.done) r.setRunning(false);
    this.pump();
    this.onUpdate();
  }

  private fail(message: string): void {
    const r = this.run;
    if (!r) return;
    r.error = message;
    r.setRunning(false);
    this.onUpdate();
  }
}
