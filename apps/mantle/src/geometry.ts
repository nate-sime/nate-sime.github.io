/**
 * The two domains, as the handful of metric quantities that distinguish them.
 *
 * Everything downstream — the spline space, the Galerkin operators, the
 * temperature transport, the GPU kernels — is written on one non-periodic axis
 * `r` crossed with one periodic axis `φ`. That is as true of a Cartesian box as
 * it is of an annulus, and the whole difference between them is the metric
 *
 *     ds² = dr² + h(r)² dφ²,
 *
 * with `h = r` on the annulus and `h = 1` in the box. So the geometry is not a
 * second solver: it is `h`, `h′`, the period of φ, and the two boundary
 * quantities that are defined by reference to the conductive state. Nothing else
 * in the pipeline knows which domain it is in.
 *
 * Writing it this way is what keeps the *scheme* the subject. `u = ∇×(ψ ẑ)` is
 * pointwise divergence-free in both, the dissipation form is the same integral
 * with a different Jacobian, and free-slip is natural to it either way — so the
 * box is the same discretisation checked on a second geometry rather than a
 * parallel code path that could drift from it.
 *
 * **The transform is why the box has two wall settings.** The transverse
 * direction is diagonalised by a radix-2 FFT, which *is* the statement that φ
 * wraps, so `periodic` is what the machinery gives for free. Free-slip side
 * walls are reached without giving any of it up, by **mirroring**: solve on a
 * period of `2L` and hold the state even about x = 0.
 *
 * That is not an approximation of a walled box, it is one. An even T gives an
 * odd ψ, hence `u_x = −ψ_z = 0` and `ψ_xx = 0` — impermeable and stress-free —
 * on both x = 0 and x = L, and `∂_x T = 0` there, the insulating sidewall the
 * benchmarks specify. Nothing in [0, L] is constrained by it: the reflection
 * only decides what the *other* half does, and the visible half keeps every
 * degree of freedom it had.
 *
 * The projection onto that subspace costs nothing, because symmetry about x = 0
 * is exactly `Im T̂ = 0`. The diffusion solve is already in mode space, so
 * dropping the imaginary half there *is* the wall — no extra pass, no
 * mirror-indexing anywhere, and the CPU and GPU do the identical thing. What it
 * does cost is a doubled period: for a given width and a given `na`, a walled
 * box resolves half as finely as a periodic one.
 *
 * ψ needs no projection of its own. An even T gives a load that is odd, hence a
 * ψ̂ that is purely imaginary, so oddness is inherited rather than imposed —
 * exactly in f64, and to f32 round-off on the GPU. Since ψ is re-solved from a
 * freshly projected T every step, that round-off cannot accumulate.
 *
 * Depth runs 0 → 1 with the *hot* boundary at 0, matching the convention in the
 * mantle convection literature and matching the annulus, where `lo = r_i` is the
 * core–mantle boundary. `lo` is therefore always the hot side and `hi` the cold
 * one, and buoyancy always points from `lo` towards `hi`.
 */

export type GeometryKind = "annulus" | "box";

/**
 * What closes the transverse direction.
 *
 * `periodic` is the domain the solver natively has. `free-slip` is that domain
 * mirrored — impermeable, stress-free and insulating at x = 0 and x = L — and
 * the annulus, having no ends, is always the former.
 */
export type Walls = "periodic" | "free-slip";

export interface Geometry {
  readonly kind: GeometryKind;
  readonly walls: Walls;
  /** Hot boundary: r_i, or z = 0. */
  readonly lo: number;
  /** Cold boundary: r_o, or z = 1. */
  readonly hi: number;
  /**
   * **Period** of the transverse axis — what the knot vector, the grid and every
   * transform are built on. 2π for the annulus, L for a periodic box, and `2L`
   * for a walled one, whose second half is the mirror image of the first.
   */
  readonly span: number;
  /**
   * **Physical width** — what is drawn, and what an aspect ratio means. Equal to
   * `span` except in a walled box, where it is half of it. The only consumer is
   * the renderer: the solver works on the whole period and does not care that
   * half of it is a reflection.
   */
  readonly width: number;
  /** h(r): arc length per unit of the transverse parameter. */
  h(r: number): number;
  /**
   * h′(r). Constant in both geometries — 1 and 0 — which is why every metric
   * term below is either present or absent rather than merely small, and why the
   * box never divides by a radius it does not have.
   */
  readonly dh: number;
  /** Steady conduction profile, 1 at `lo` and 0 at `hi`. The Nu = 1 reference. */
  conduction(r: number): number;
  /**
   * The factor turning a summed boundary flux into a Nusselt number:
   * `Nu = nuScale(r) · Σ_j q_j · dφ`. It is `h(r)` over the conductive state's
   * total flux through the same boundary, so `Nu = 1` for pure conduction at
   * both ends by construction.
   */
  nuScale(r: number): number;
}

/**
 * Spherical annulus, nondimensionalised by the mantle's own thickness rather
 * than its outer radius: `r_o − r_i = 1` by default, matching the unit depth
 * of the box. The defaults are Earth's core–mantle boundary and surface, 3486
 * and 6371 km, divided through by the 2885 km between them; `ui/dimensional.ts`
 * reads the same choice from the other side.
 */
export const annulus = (ri = 1.208318891, ro = ri + 1): Geometry => {
  const d = Math.log(ro / ri);
  return {
    kind: "annulus", walls: "periodic", lo: ri, hi: ro,
    span: 2 * Math.PI, width: 2 * Math.PI, dh: 1,
    h: (r) => r,
    conduction: (r) => Math.log(ro / r) / d,
    // Conductive flux at radius r is 1/(r ln(r_o/r_i)), over a boundary of
    // length 2πr — so the total is 2π/ln(r_o/r_i) at both radii, as it must be.
    nuScale: (r) => (r * d) / (2 * Math.PI),
  };
};

/**
 * Cartesian box of width `length`, depth 0 → 1: the cell the mantle convection
 * literature states its benchmarks in.
 *
 * `walls` picks what closes it left and right. A walled box is solved on a
 * period of `2·length` and held even about x = 0 — see the module header for why
 * that is the walled problem exactly rather than a stand-in for it.
 *
 * `nuScale = 1/span` needs no case of its own: the flux is summed over the whole
 * period, so `1/span` is its mean there, and the mean over a mirrored domain is
 * the mean over either half.
 */
export const box = (length: number, walls: Walls = "periodic"): Geometry => {
  const span = walls === "free-slip" ? 2 * length : length;
  return {
    kind: "box", walls, lo: 0, hi: 1, span, width: length, dh: 0,
    h: () => 1,
    conduction: (z) => 1 - z,
    nuScale: () => 1 / span,
  };
};

/** The geometry every verification and manufactured solution is written on. */
export const ANNULUS = annulus();

/**
 * What to call the two boundaries on screen.
 *
 * The code names them `inner` and `outer` everywhere, because that is what they
 * are on the annulus and because the series keys of the Nusselt trace are those
 * words. A box has no inner or outer — it has a floor and a ceiling — and a
 * legend reading "inner 4.21" over a rectangle is a label the reader has to
 * translate. So the *names* are a display concern and live here, next to the
 * geometry that decides them, rather than being spelled into either panel.
 */
export const boundaryNames = (
  kind: GeometryKind,
): { inner: string; outer: string } =>
  kind === "annulus"
    ? { inner: "inner", outer: "outer" }
    : { inner: "bottom", outer: "top" };

/**
 * Transverse wavenumber of DFT mode `k`: `∂_φφ → −wavenumber(g, k)²`.
 *
 * On the annulus φ spans 2π and this is just `k`, which is why the existing
 * kernels could write `k²` and be right. In a box φ spans L, so mode k is
 * `exp(2πikx/L)` and the multiplier is `(2πk/L)²`. The *discrete* azimuthal
 * symbols in `solver/operators.ts` need no such correction — they are computed
 * from the knot vector, so they already carry the period.
 */
export const wavenumber = (g: Geometry, k: number): number =>
  (2 * Math.PI * k) / g.span;
