// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * A Monte Carlo worker: holds one `Sampler` and runs whatever range of sample
 * indices it is handed. It keeps no randomness of its own — sample i is drawn
 * from its index (`random/philox.ts`) — so the pool may hand any range to any
 * worker. Messages are tagged with a job id; work for a superseded job is
 * dropped unanswered.
 */

import type { KLData } from "../random/kl";
import { FIELD_POINTS, Sampler, type McSpec } from "./sampler";
import type { Batch } from "./stats";

export interface JobMessage {
  readonly type: "job";
  readonly id: number;
  readonly spec: McSpec;
  readonly kl: KLData;
  readonly ne: number;
  readonly coarse: boolean;
  readonly wantW: boolean;
}

export interface RunMessage {
  readonly type: "run";
  readonly id: number;
  readonly from: number;
  readonly count: number;
}

export type Reply =
  | { readonly type: "batch"; readonly id: number; readonly batch: Batch }
  | { readonly type: "error"; readonly id: number; readonly message: string };

// The DOM lib types `self` as a Window; this is the slice of a worker scope used here.
const scope = self as unknown as {
  onmessage: ((e: MessageEvent<JobMessage | RunMessage>) => void) | null;
  postMessage(m: Reply, transfer?: Transferable[]): void;
};

let job: JobMessage | null = null;
let sampler: Sampler | null = null;

scope.onmessage = (e) => {
  const m = e.data;
  if (m.type === "job") {
    job = m;
    sampler = new Sampler(m.spec, m.kl);
    return;
  }
  if (!job || !sampler || m.id !== job.id) return;
  try {
    const t0 = performance.now();
    const Q = new Float64Array(m.count), Qc = new Float64Array(m.count);
    const W = job.wantW ? new Float64Array(m.count * FIELD_POINTS) : null;
    for (let s = 0; s < m.count; s++) {
      const r = sampler.sample(m.from + s, job.ne, job.coarse, job.wantW);
      Q[s] = r.Q;
      Qc[s] = r.Qc;
      if (W && r.w) W.set(r.w, s * FIELD_POINTS);
    }
    const batch: Batch = { from: m.from, count: m.count, Q, Qc, W, ms: performance.now() - t0 };
    scope.postMessage({ type: "batch", id: m.id, batch }, W ? [Q.buffer, Qc.buffer, W.buffer] : [Q.buffer, Qc.buffer]);
  } catch (err) {
    scope.postMessage({ type: "error", id: m.id, message: (err as Error).message });
  }
};
