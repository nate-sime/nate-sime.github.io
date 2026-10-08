// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The pool of Monte Carlo workers, and the runs it works on.
 *
 * Samples are independent, so the parallelism is across them: the pool hands
 * each idle worker the next range of indices of some stream and folds what
 * comes back into that stream's `Accumulator`, which restores sample order.
 * A run is one job — a spec and its random field, which every worker holds —
 * and one or more streams: plain Monte Carlo has one (a level and its parent);
 * multilevel Monte Carlo has one per level, each with its own random stream.
 * Batch sizes adapt per stream toward ~50 ms of work — small at first, so the
 * first statistics appear at once, larger once the per-sample cost is known, so
 * message overhead stays negligible.
 *
 * Only the active run is dispatched; a run left (its view hidden, or an input
 * changed) keeps what it has, and a few are kept so that going back resumes
 * rather than restarts. Every reply names its run, and workers answer in the
 * order they were asked, so a range in flight when the job changed still lands
 * in the run it was taken from. Raising a stream's target extends it in place.
 */

import type { KLData } from "../random/kl";
import { FIELD_POINTS, type McSpec } from "./sampler";
import { Accumulator } from "./stats";
import type { JobMessage, Reply, RunMessage, StreamSpec } from "./worker";

export type { StreamSpec } from "./worker";

export interface McJob {
  readonly spec: McSpec;
  readonly kl: KLData;
  readonly klY?: KLData;
}

const BATCH_MS = 50;
/** Runs kept besides the active one. */
const KEEP = 4;

interface Slot {
  readonly worker: Worker;
  busy: boolean;
}

export class McStream {
  readonly acc: Accumulator;
  target = 0;
  /** Samples handed out so far. */
  dispatched = 0;
  /** Samples per message, adapted toward BATCH_MS of work. */
  size = 4;
  /** Worker time per sample, ms; NaN until a batch has come back. */
  msPerSample = NaN;

  constructor(readonly spec: StreamSpec) {
    this.acc = new Accumulator(spec.wantW ? FIELD_POINTS : 0);
  }

  get done(): boolean {
    return this.acc.n >= this.target;
  }
}

export class McRun {
  readonly streams: McStream[] = [];
  running = false;
  /** Paused by the reader, as opposed to merely not shown. */
  paused = false;
  error: string | null = null;
  /** Wall time spent running, ms (pauses excluded). */
  private wall = 0;
  private since = 0;

  constructor(readonly id: number, readonly key: string, readonly job: McJob) {}

  /** The stream solving `spec`, added (with no target yet) if new. */
  stream(spec: StreamSpec): McStream {
    const same = (s: StreamSpec) => s.ne === spec.ne && s.coarse === spec.coarse && s.wantW === spec.wantW && s.stream === spec.stream;
    let s = this.streams.find((t) => same(t.spec));
    if (!s) this.streams.push((s = new McStream(spec)));
    return s;
  }

  /** Plain Monte Carlo's one stream: its statistics. */
  get acc(): Accumulator {
    return this.streams[0].acc;
  }

  get target(): number {
    return this.streams[0]?.target ?? 0;
  }

  get done(): boolean {
    return this.streams.every((s) => s.done);
  }

  get wallMs(): number {
    return this.wall + (this.running ? performance.now() - this.since : 0);
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
  private runs: McRun[] = [];
  private active: McRun | null = null;
  /** Id of the job the workers were last sent. */
  private held = -1;
  private nextId = 1;

  constructor(private readonly onUpdate: () => void, size?: number) {
    this.size = size ?? Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 4) - 1));
    this.slots = Array.from({ length: this.size }, () => {
      const worker = new Worker(new URL("./worker.ts", import.meta.url), { type: "module" });
      const slot: Slot = { worker, busy: false };
      worker.onmessage = (e: MessageEvent<Reply>) => this.receive(slot, e.data);
      worker.onerror = (e) => this.fail(this.active, e.message || "worker failed to load");
      return slot;
    });
  }

  /**
   * The run for `key` — kept from before if there is one, else a fresh run of
   * `job()` — made the active one. It runs whenever it has work, unless paused.
   */
  open(key: string, job: () => McJob): McRun {
    let r = this.runs.find((x) => x.key === key);
    if (!r) {
      r = new McRun(this.nextId++, key, job());
      this.runs.push(r);
      const old = this.runs.filter((x) => x !== r && x !== this.active);
      for (const x of old.slice(0, Math.max(0, old.length - KEEP))) this.runs.splice(this.runs.indexOf(x), 1);
    }
    if (r !== this.active) {
      this.active?.setRunning(false);
      this.active = r;
      if (this.held !== r.id) {
        const m: JobMessage = { type: "job", id: r.id, ...r.job };
        for (const s of this.slots) s.worker.postMessage(m);
        this.held = r.id;
      }
    }
    this.refresh(r);
    return r;
  }

  /** Plain Monte Carlo: the run for `key`, with one stream raised (or lowered) to `target`. */
  ensure(key: string, job: () => McJob, spec: StreamSpec, target: number): McRun {
    const r = this.open(key, job);
    this.demand(r, r.stream(spec), target);
    return r;
  }

  /** Set a stream's target; work starts at once if its run is the active one. */
  demand(r: McRun, s: McStream, target: number): void {
    s.target = target;
    this.refresh(r);
  }

  current(): McRun | null {
    return this.active;
  }

  /** Stop dispatching, keeping every run: their view is hidden. */
  idle(): void {
    this.active?.setRunning(false);
    this.active = null;
  }

  pause(): void {
    if (this.active) { this.active.paused = true; this.refresh(this.active); }
  }

  resume(): void {
    if (this.active) { this.active.paused = false; this.refresh(this.active); }
  }

  private refresh(r: McRun): void {
    r.setRunning(r === this.active && !r.paused && !r.error && !r.done);
    this.pump();
  }

  /**
   * Each idle worker gets a range of the stream with the least work left, so
   * small demands — a survey, a fine level's handful of samples — are met at
   * once and the statistics that wait on them appear, while a coarse level's
   * long run takes every worker it can once they are.
   */
  private pump(): void {
    const r = this.active;
    if (!r || !r.running) return;
    for (const slot of this.slots) {
      if (slot.busy) continue;
      let best: McStream | null = null, least = Infinity;
      for (const s of r.streams) {
        const left = (s.target - s.dispatched) * (Number.isFinite(s.msPerSample) ? s.msPerSample : s.spec.ne);
        if (s.dispatched < s.target && left < least) { best = s; least = left; }
      }
      if (!best) return;
      const count = Math.min(best.size, best.target - best.dispatched);
      const m: RunMessage = { type: "run", id: r.id, tag: r.streams.indexOf(best), from: best.dispatched, count, ...best.spec };
      best.dispatched += count;
      slot.busy = true;
      slot.worker.postMessage(m);
    }
  }

  private receive(slot: Slot, m: Reply): void {
    slot.busy = false;
    const r = this.runs.find((x) => x.id === m.id);
    if (r) {
      if (m.type === "error") { this.fail(r, m.message); return; }
      const s = r.streams[m.tag], b = m.batch, per = b.ms / b.count;
      s.msPerSample = Number.isFinite(s.msPerSample) ? 0.8 * s.msPerSample + 0.2 * per : per;
      s.size = Math.max(1, Math.min(4096, Math.round(BATCH_MS / Math.max(per, 1e-3)), 4 * s.size));
      s.acc.push(b);
      if (r === this.active) this.refresh(r);
    }
    this.pump();
    this.onUpdate();
  }

  private fail(r: McRun | null, message: string): void {
    if (!r) return;
    r.error = message;
    r.setRunning(false);
    this.onUpdate();
  }
}
