/**
 * RMS velocity against time — the data behind the second corner plot, and the
 * scale it is drawn on.
 *
 * Nu says whether heat transport has settled; it says nothing about the flow
 * carrying it, and two runs can share a Nusselt number while one is a slow
 * single roll and the other a fast tangle of plumes. `v_rms = √⟨|u|²⟩` is that
 * missing half — the area-weighted reduction in `Temperature.rmsVelocity` and
 * `gpu/wgsl.ts`'s `rmsSource`, read as a series rather than an instant.
 *
 * **The area weight is what "accounts for the coordinate system change" in
 * the header this module answers to.** `⟨·⟩` is an average over the domain,
 * and the domain's area element is `h(r) dr dφ` — the ring of circumference
 * `2πr` on the annulus, uniform in a box. A plain grid mean would let a fixed
 * number of azimuthal samples near the inner boundary outweigh the same count
 * near the outer one, on the annulus alone; weighting by `h` is what makes a
 * single number comparable *between* the two geometries as well as *within*
 * either one — the same discipline `Geometry.nuScale` applies to Nu.
 *
 * Two series, unlike the header above once said: `v` is the domain-average
 * flow speed, and `vs` sits beside it as the same reduction taken over the
 * top boundary alone — `outer` in `boundaryNames`, a real mantle's surface
 * and the one place plate motion is actually observed. They are not two
 * readings of the same number the way Nu's inner and outer are, and are not
 * expected to converge; `vs` is usually the smaller of the two, since a
 * no-penetration top drops the radial component of `u` to zero there and
 * leaves only the tangential part `rmsVelocity`'s area weight lets the
 * interior's faster plume cores dominate. Sharing one trace and one axis
 * still follows `NuTrace`'s reasoning rather than departing from it: both are
 * the same quantity, in the same units, from the same poll — worth reading
 * against each other, which a second y scale would defeat. The ring buffer,
 * the window and the axis padding are still exactly `NuTrace`'s, which is why
 * they are imported rather than re-derived.
 *
 * DOM-free for the same reason as `nusselt.ts`: `rmsplot.ts` renders it.
 */

import { axisDecimals, niceAxis, type NuAxis, type NuExtent } from "./nusselt";

export { axisDecimals, niceAxis, type NuAxis, type NuExtent };

/**
 * Categorical slot 3 (aqua) of the site's reference palette, stepped for the
 * dark surface — the next slot after the Nusselt plot's 1 and 2 (`NU_COLOUR`),
 * not a colour picked to taste. Validated as a pair against both Nu series and
 * against this panel's own translucent surface: worst-case CVD ΔE 9.4,
 * normal-vision 26.5, and ≥ 3:1 contrast against the ground the panel's 8%
 * transparency can produce — the same bar `nusselt.ts` holds its pair to.
 */
export const RMS_COLOUR = "#199e70";

/**
 * Categorical slot 4 (violet) — the next slot after `RMS_COLOUR`'s 3, and
 * clear of both Nu series (orange, blue) and the domain curve (aqua) on hue
 * alone, so the top-boundary curve reads as its own series rather than a
 * shade of one already on screen.
 */
export const RMS_SURFACE_COLOUR = "#8b5cf6";

export interface RmsSample {
  /** Simulation time the reduction was taken at. */
  t: number;
  /** Step count at that time — what the display window is measured in. */
  step: number;
  /** Domain-average √⟨|u|²⟩, h-weighted — `Temperature.rmsVelocity`'s twin. */
  v: number;
  /** The same reduction over the top boundary alone — `rmsSurfaceVelocity`'s. */
  vs: number;
}

/**
 * Samples retained, matched to `NU_CAPACITY`: the two traces are polled by the
 * same frame loop at the same rate, so there is no reason for one to remember
 * more history than the other.
 */
export const RMS_CAPACITY = 16_384;

/** A rolling buffer of the last `capacity` polls — `NuTrace` with `v`/`vs` in place of inner/outer. */
export class RmsTrace {
  private readonly ts: Float64Array;
  private readonly st: Float64Array;
  private readonly v: Float32Array;
  private readonly vs: Float32Array;
  private head = 0;
  private count = 0;

  constructor(readonly capacity = RMS_CAPACITY) {
    if (!Number.isInteger(capacity) || capacity < 2)
      throw new Error(`capacity ${capacity} must be an integer ≥ 2`);
    this.ts = new Float64Array(capacity);
    this.st = new Float64Array(capacity);
    this.v = new Float32Array(capacity);
    this.vs = new Float32Array(capacity);
  }

  get length(): number { return this.count; }

  clear(): void { this.head = 0; this.count = 0; }

  /** Record one poll; see `NuTrace.push` for the drop and restart rules. */
  push({ t, step, v, vs }: RmsSample): boolean {
    if (!(Number.isFinite(t) && Number.isFinite(step)
          && Number.isFinite(v) && Number.isFinite(vs)))
      return false;
    if (this.count > 0) {
      const k = (this.head + this.capacity - 1) % this.capacity;
      if (t < this.ts[k] || step < this.st[k]) this.clear();
      else if (t === this.ts[k]) return false;
    }
    this.ts[this.head] = t;
    this.st[this.head] = step;
    this.v[this.head] = v;
    this.vs[this.head] = vs;
    this.head = (this.head + 1) % this.capacity;
    if (this.count < this.capacity) this.count++;
    return true;
  }

  /** Index of the oldest sample lying within `steps` of the newest. */
  first(steps: number): number {
    if (this.count === 0 || !(steps > 0)) return 0;
    if (!Number.isFinite(steps)) return 0;
    const cutoff = this.st[this.index(this.count - 1)] - steps;
    let i = this.count - 1;
    while (i > 0 && this.st[this.index(i - 1)] >= cutoff) i--;
    return i;
  }

  /** Sample `i`, counting from the oldest held. */
  at(i: number): RmsSample {
    const k = this.index(i);
    return { t: this.ts[k], step: this.st[k], v: this.v[k], vs: this.vs[k] };
  }

  get last(): RmsSample | null {
    return this.count === 0 ? null : this.at(this.count - 1);
  }

  /** Bounds of the samples from `from` onwards, over *both* series — they share the axis. */
  extent(from = 0): NuExtent | null {
    if (from >= this.count || from < 0) return null;
    let lo = Infinity, hi = -Infinity;
    for (let i = from; i < this.count; i++) {
      const k = this.index(i);
      lo = Math.min(lo, this.v[k], this.vs[k]);
      hi = Math.max(hi, this.v[k], this.vs[k]);
    }
    return {
      t0: this.ts[this.index(from)],
      t1: this.ts[this.index(this.count - 1)],
      lo, hi,
    };
  }

  private index(i: number): number {
    if (!Number.isInteger(i) || i < 0 || i >= this.count)
      throw new RangeError(`sample ${i} of ${this.count}`);
    return (this.head - this.count + i + this.capacity) % this.capacity;
  }
}
