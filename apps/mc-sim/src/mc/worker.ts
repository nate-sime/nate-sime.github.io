// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * A Monte Carlo worker: holds one `Sampler` and runs whatever range of sample
 * indices it is handed, on whatever level and stream the message names. It
 * keeps no randomness of its own — sample i is drawn from its index
 * (`random/philox.ts`) — so the pool may hand any range to any worker. Messages
 * carry the id of the job they belong to; a range for a job other than the one
 * held is dropped unanswered. Messages arrive in order, so every range sent
 * before a new job is still run under the old one.
 */

import type { KLData } from "../random/kl";
import { FIELD_POINTS, Sampler, type McSpec } from "./sampler";
import type { Batch } from "./stats";

export interface JobMessage {
  readonly type: "job";
  readonly id: number;
  readonly spec: McSpec;
  readonly kl: KLData;
  /** The expansion along a plate's y; absent for a beam or a square plate. */
  readonly klY?: KLData;
}

/** Where a stream's samples are solved: a level, whether its parent too, and which random stream. */
export interface StreamSpec {
  /** Elements on the level; its parent has ne/2. */
  readonly ne: number;
  readonly coarse: boolean;
  readonly wantW: boolean;
  /** Philox stream the samples are drawn from. */
  readonly stream: number;
}

export interface RunMessage extends StreamSpec {
  readonly type: "run";
  readonly id: number;
  /** Index of the stream within its run, echoed in the reply. */
  readonly tag: number;
  readonly from: number;
  readonly count: number;
}

export type Reply =
  | { readonly type: "batch"; readonly id: number; readonly tag: number; readonly batch: Batch }
  | { readonly type: "error"; readonly id: number; readonly message: string };

// The DOM lib types `self` as a Window; this is the slice of a worker scope used here.
const scope = self as unknown as {
  onmessage: ((e: MessageEvent<JobMessage | RunMessage>) => void) | null;
  postMessage(m: Reply, transfer?: Transferable[]): void;
};

let id = -1;
let sampler: Sampler | null = null;

scope.onmessage = (e) => {
  const m = e.data;
  if (m.type === "job") {
    id = m.id;
    sampler = new Sampler(m.spec, m.kl, m.klY);
    return;
  }
  if (!sampler || m.id !== id) return;
  try {
    const t0 = performance.now();
    const Q = new Float64Array(m.count), Qc = new Float64Array(m.count);
    const W = m.wantW ? new Float64Array(m.count * FIELD_POINTS) : null;
    for (let s = 0; s < m.count; s++) {
      const r = sampler.sample(m.from + s, m.ne, m.coarse, m.wantW, m.stream);
      Q[s] = r.Q;
      Qc[s] = r.Qc;
      if (W && r.w) W.set(r.w, s * FIELD_POINTS);
    }
    const batch: Batch = { from: m.from, count: m.count, Q, Qc, W, ms: performance.now() - t0 };
    scope.postMessage({ type: "batch", id: m.id, tag: m.tag, batch }, W ? [Q.buffer, Qc.buffer, W.buffer] : [Q.buffer, Qc.buffer]);
  } catch (err) {
    scope.postMessage({ type: "error", id: m.id, message: (err as Error).message });
  }
};
