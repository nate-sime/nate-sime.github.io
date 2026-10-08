// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * Reproducible randomness, addressed rather than drawn.
 *
 * Philox4x32-10 (Salmon, Moraes, Dror & Shaw, "Parallel random numbers: as easy
 * as 1, 2, 3", SC 2011; the generator of Random123, cuRAND and JAX) is a keyed
 * bijection of a 128-bit counter: ten rounds of multiply, xor and key bump. It
 * keeps no state, so the j-th normal of sample i is a pure function of
 *
 *   (seed, i, channel, j),
 *
 * which is the property multilevel Monte Carlo needs — sample ω must be the
 * same ω on the fine and the coarse solve, in whichever worker either runs, in
 * whatever order — and the one a GPU port needs, one invocation per sample with
 * nothing shared. A stateful generator (PCG, xoshiro) would have to be jumped
 * ahead per sample to the same effect.
 *
 * The counter is (block, i, channel, stream) and the key (seed, 0): each block
 * gives four 32-bit words, hence four normals by Box–Muller. A channel is an
 * independent sequence per random input (stiffness, load), so changing how many
 * normals one input reads never shifts another's. A stream is an independent
 * copy of the whole experiment: multilevel Monte Carlo draws level ℓ's samples
 * from stream ℓ, so the estimators on different levels are independent while
 * the fine and coarse solves within a level share their ω. Plain Monte Carlo is
 * stream 0.
 */

const M0 = 0xd2511f53, M1 = 0xcd9e8d57;
const W0 = 0x9e3779b9, W1 = 0xbb67ae85;

/** High 32 bits of the 64-bit product of two uint32s, in exact double arithmetic. */
function mulhi(a: number, b: number): number {
  const ah = a >>> 16, al = a & 0xffff, bh = b >>> 16, bl = b & 0xffff;
  const ll = al * bl, lh = al * bh, hl = ah * bl;
  const mid = (ll >>> 16) + (lh & 0xffff) + (hl & 0xffff);
  return (ah * bh + (lh >>> 16) + (hl >>> 16) + (mid >>> 16)) >>> 0;
}

/** Philox4x32-10 of counter (c0, c1, c2, c3) under key (k0, k1), into `out`. */
export function philox4x32(
  c0: number, c1: number, c2: number, c3: number, k0: number, k1: number, out: Uint32Array,
): Uint32Array {
  c0 >>>= 0; c1 >>>= 0; c2 >>>= 0; c3 >>>= 0; k0 >>>= 0; k1 >>>= 0;
  for (let r = 0; r < 10; r++) {
    const hi0 = mulhi(M0, c0), lo0 = Math.imul(M0, c0) >>> 0;
    const hi1 = mulhi(M1, c2), lo1 = Math.imul(M1, c2) >>> 0;
    c0 = (hi1 ^ c1 ^ k0) >>> 0;
    c1 = lo1;
    c2 = (hi0 ^ c3 ^ k1) >>> 0;
    c3 = lo0;
    k0 = (k0 + W0) >>> 0;
    k1 = (k1 + W1) >>> 0;
  }
  out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
  return out;
}

/** A uint32 as a uniform in the open interval (0, 1) — never 0, so a logarithm of it is finite. */
export const unit = (u: number): number => (u + 0.5) * 2 ** -32;

/** Independent streams, one per random input. */
export const CHANNEL = { stiffness: 0, load: 1 } as const;

/**
 * The first `n` standard normals of sample `index` on `channel` of `stream`, by
 * Box–Muller on consecutive pairs of uniforms. The j-th is the same whatever
 * `n` is.
 */
export function normals(
  seed: number, index: number, channel: number, n: number, stream = 0, out = new Float64Array(n),
): Float64Array {
  const w = new Uint32Array(4);
  for (let b = 0; 4 * b < n; b++) {
    philox4x32(b, index, channel, stream, seed, 0, w);
    for (let h = 0; h < 2; h++) {
      const r = Math.sqrt(-2 * Math.log(unit(w[2 * h]))), t = 2 * Math.PI * unit(w[2 * h + 1]);
      const j = 4 * b + 2 * h;
      if (j < n) out[j] = r * Math.cos(t);
      if (j + 1 < n) out[j + 1] = r * Math.sin(t);
    }
  }
  return out;
}

/** The first `n` uniforms in (0, 1) of sample `index` on `channel` of `stream`. */
export function uniforms(
  seed: number, index: number, channel: number, n: number, stream = 0, out = new Float64Array(n),
): Float64Array {
  const w = new Uint32Array(4);
  for (let b = 0; 4 * b < n; b++) {
    philox4x32(b, index, channel, stream, seed, 0, w);
    for (let r = 0; r < 4 && 4 * b + r < n; r++) out[4 * b + r] = unit(w[r]);
  }
  return out;
}
