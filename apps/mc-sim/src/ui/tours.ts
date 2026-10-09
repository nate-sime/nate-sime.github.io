// This file is part of the MC-sim app, Monte Carlo for vibrating structures.
// Copyright (c) 2026 Nathan Sime
// SPDX-License-Identifier: MIT

/**
 * The guided tours, as data — the four of MC_PLAN.md's stage 9:
 *
 *   Monte Carlo and the √N rate
 *   discretisation hierarchies (h, p and k refinement)
 *   the bias–variance split
 *   multilevel Monte Carlo and its complexity gain
 *
 * A step names the control it is about, the settings it wants (a patch over
 * whatever the steps before it left), and what the reader should watch for;
 * `tour.ts` makes each of those true. No DOM here, so `tests/tour.test.ts`
 * can check every invariant a tour relies on — that a patch only sets fields
 * `State` has, to values the pane can show, and that every control a step
 * lights is on screen in the state that step leaves — without a browser,
 * rather than part-way through a tour a reader clicks through once.
 *
 * Every number quoted in a card is one measured and recorded in the README —
 * a tour states what the app will show, not what theory hopes.
 */

import type { State } from "./state";
import type { ControlName } from "./visibility";

/** A bulleted list in a card's body. */
export interface TourList {
  readonly items: readonly string[];
}

/**
 * How long a step holds before "next" is offered as the thing to press: wall
 * time for a step that only asks the reader to look, or a sample count for one
 * waiting on a run — the same statistics on a slow machine and a fast one.
 */
export type TourDwell = { readonly ms: number } | { readonly samples: number };

export interface TourStep {
  /** Stable across edits to the prose: what tests and logs name a step by. */
  readonly id: string;
  readonly title: string;
  /** One paragraph or list per entry. */
  readonly body: readonly (string | TourList)[];
  /**
   * The pane control lit while the rest of the pane dims; null lights none.
   * The figure and the readout are never dimmed — they are what every step
   * asks the reader to look at.
   */
  readonly target: ControlName | null;
  /** Settings this step makes, over everything the steps before it made. */
  readonly patch?: Partial<State>;
  /** What to look for, set apart from why. */
  readonly watch?: string;
  readonly dwell?: TourDwell;
}

/**
 * Where every tour starts: the app's defaults for everything a tour's views
 * read, so a tour lands the same whatever the reader had set before it. Each
 * tour's first step applies it, then its own view on top. One exception: the
 * survey is Giles' 2·10³ per level, not the app's 10², since the MLMC tour's
 * figures (plain MC's count from the survey's V[Q₃]) are measured with it.
 */
export const TOUR_BASE: Partial<State> = {
  dimensional: false, structure: "beam",
  p: 3, continuity: "max", ne: 6,
  supports: "cantilever", section: "uniform", load: "uniform", motion: "mode", mesh: false,
  edges: "SSSS", aspect: 1, nu: 0.3,
  qoi: "omega1", ne0: 4, levels: 6, forceRatio: 0.8, zeta: 0.02,
  kernel: "matern32", ell: 0.2, sigma: 0.3, terms: 24, massFollows: false, loadSigma: 0, seed: 1,
  mcLevel: 2, mcSamples: 10000, mlSurvey: 2000, mlEps: 3e-4, liveLevel: 2, cmpTol: 4,
};

const MONTE_CARLO: readonly TourStep[] = [
  {
    id: "mc-intro",
    title: "A random beam",
    body: [
      "The cantilever's stiffness EI(x) is a random field: lognormal, Matérn-correlated over a fifth of the span, varying by about 30%. So its first natural frequency ω₁ is a random number, and the question is its mean.",
      "Monte Carlo answers it the plainest way: draw N independent beams, solve each, and average. The workers are doing that now, on a 16-element cubic spline mesh.",
    ],
    target: null,
    patch: { ...TOUR_BASE, view: "montecarlo", mcSamples: 1000 },
    watch: "Top left: the histogram of ω₁ filling in. Its spread is the input's randomness passed through the beam; the white line is the estimate, the green the frequency of the mean beam.",
    dwell: { samples: 1000 },
  },
  {
    id: "mc-interval",
    title: "The estimate and its error bar",
    body: [
      "The sample mean is itself random. By the central limit theorem it is close to normal, with standard deviation σ/√N — σ the spread of ω₁ — so Q̄_N ± 1.96 σ̂/√N is a 95% interval for the true mean, computed from the samples alone.",
      "This step raises N tenfold, to 10⁴.",
    ],
    target: "mcSamples",
    patch: { mcSamples: 10000 },
    watch: "Top right: the running mean settles and its shaded interval narrows. Ten times the samples makes it √10 ≈ 3.2 times narrower — not ten.",
    dwell: { samples: 10000 },
  },
  {
    id: "mc-rate",
    title: "The √N rate",
    body: [
      "Bottom left is the same thing on log axes: the measured σ̂/√N falls as a straight line of slope −½. Every extra digit of accuracy costs a hundred times the samples.",
      "That rate is the whole price list of Monte Carlo, and the reason everything after this tour exists.",
    ],
    target: null,
    watch: "The blue line's slope: down one decade across two decades of N.",
    dwell: { ms: 4000 },
  },
  {
    id: "mc-seed",
    title: "A new seed is a new experiment",
    body: [
      "Sample i is a pure function of (seed, i), so a run is reproducible to the last bit. Changing the seed draws a different set of beams.",
      "The estimate moves — by about its own error bar. A 95% interval should miss the true mean one time in twenty; the test suite runs 300 independent experiments and checks that the coverage comes out between 88% and 99%.",
    ],
    target: "seed",
    patch: { seed: 2 },
    watch: "The readout's E[Q_ℓ] ± …: compared with seed 1, the two intervals overlap.",
    dwell: { samples: 10000 },
  },
  {
    id: "mc-sigma",
    title: "More randomness, more samples",
    body: [
      "Double the standard deviation of log EI. The spread of ω₁ roughly doubles, and with it the error bar at any N: to win the same accuracy back takes four times the samples.",
      "Error ∝ σ/√N. Monte Carlo can only be made cheaper by making σ smaller or each sample cheaper — the second is what multilevel Monte Carlo does.",
    ],
    target: "sigma",
    patch: { sigma: 0.6 },
    watch: "The histogram widens, and the readout's coefficient of variation roughly doubles.",
    dwell: { samples: 10000 },
  },
  {
    id: "mc-dimension",
    title: "The rate ignores the dimension",
    body: [
      "The random field is 24 independent normals here — a 24-dimensional integral. Take it to 64. Quadrature on a grid would need exponentially more points; Monte Carlo's error is still σ/√N.",
      "That indifference to dimension is why Monte Carlo, slow as it is, is the method for random fields. Here the estimate barely moves either: the 40 new terms hold under a tenth of a per cent of the variance, and each sample's first 24 normals are the ones it drew before — normal j of sample i is addressed, not drawn in sequence.",
    ],
    target: "terms",
    patch: { sigma: 0.3, terms: 64 },
    watch: "The error-against-N line keeps its slope of −½.",
    dwell: { samples: 10000 },
  },
  {
    id: "mc-end",
    title: "What it costs",
    body: [
      "Monte Carlo's estimate is unbiased for the mesh it samples, slow (σ/√N), and blind to the dimension of the input. What it cannot see is the mesh: every sample here was solved on 16 elements, and the average converges to the answer for that mesh, not for the beam.",
      { items: ["the bias–variance tour: how far that mesh is from the beam", "the multilevel tour: how to stop paying the fine mesh's price for every sample"] },
    ],
    target: null,
  },
];

const HIERARCHY: readonly TourStep[] = [
  {
    id: "h-intro",
    title: "A ladder of meshes",
    body: [
      "Before any randomness: one beam, clamped at both ends, solved on meshes of ne₀·2^ℓ elements, ℓ = 0 … 5, at degrees 2 and 3. The quantity is ω₁, against its closed form.",
      "Each line is h-refinement — the same space on a finer mesh. On log–log axes the error falls as h^α, and α is the slope.",
    ],
    target: null,
    patch: { ...TOUR_BASE, view: "convergence", supports: "clamped–clamped", p: 2 },
    watch: "The left panel's legend, the beam's: fitted α against theory 2(p − 1), measured 2.00 and 3.98 for p = 2 and 3.",
    dwell: { ms: 4000 },
  },
  {
    id: "h-p",
    title: "p-refinement",
    body: [
      "Moving between lines at a fixed h is p-refinement: raise the degree on the same mesh. The error drops by orders of magnitude, and the slope steepens with it, since a frequency converges at twice the energy rate, h^(2(p−1)).",
      "The table below the plot follows the degree chosen here.",
    ],
    target: "p",
    patch: { p: 3 },
    watch: "The table: |Q_ℓ − Q|/|Q| down each column, a factor 2⁴ = 16 per level at p = 3, where p = 2 managed 4.",
    dwell: { ms: 3500 },
  },
  {
    id: "h-k",
    title: "k-refinement",
    body: [
      "Continuity is the third axis. C¹ at degree 3 is the Hermite cubic, the classic beam element: two new coefficients per element. Maximal continuity, C², is the smooth spline: one per element. Same degree, same rate, 2(p − 1) = 4.",
      "On the same mesh the two are as accurate as each other — 3.3·10⁻⁷ on 32 elements — but C¹ needs twice the unknowns. Per unknown the smooth spline wins outright: C² on 64 elements has 63 unknowns and an error of 2.1·10⁻⁸; C¹ on 32 elements has 62 and an error of 3.3·10⁻⁷.",
    ],
    target: "continuity",
    patch: { continuity: "c1" },
    watch: "The dofs and |Q_ℓ − Q|/|Q| columns, against the same columns a step ago.",
    dwell: { ms: 3500 },
  },
  {
    id: "h-roundoff",
    title: "Round-off is a level of the hierarchy",
    body: [
      "The same ladder, on a cantilever. The stiffness of a fourth-order problem is conditioned like h⁻⁴, so f64 round-off in ω₁ grows as the mesh refines. The grey line is that round-off, measured per level by re-solving with every entry perturbed by (p + 1)ε.",
      "Over eight levels at p = 3 the error meets it at 64 elements, 4.7·10⁻¹⁰ against a round-off of 1.5·10⁻⁹. Past 128 the measured error is round-off and climbs with it: 10⁻⁷ at 256 elements, 2.5·10⁻⁶ at 512. A finer mesh buys nothing there; the rates are fitted only where they stand clear of it.",
    ],
    target: "supports",
    patch: { view: "convergence", supports: "cantilever", p: 3, continuity: "max", levels: 8 },
    watch: "The left panel's p = 3 line turning back up along the grey round-off line.",
    dwell: { ms: 4000 },
  },
  {
    id: "h-plate",
    title: "The plate: the same rates, a steeper price",
    body: [
      "On the right, the Kirchhoff plate, simply supported. The rates are the beam's — 2.01 and 4.03 measured for p = 2 and 3 — but a banded solve now costs dofs × bandwidth² ∝ h⁻⁴: γ = 4 against the beam's 1.",
      "This view stops a plate at 32 × 32 elements; that is where the cost tells.",
    ],
    target: "edges",
    patch: { edges: "SSSS", p: 3 },
    watch: "The plate readout's γ line: 3.6 on these meshes, rising toward 4 as they refine.",
    dwell: { ms: 5000 },
  },
  {
    id: "h-corners",
    title: "Corners cap the rate",
    body: [
      "A cantilever plate: clamped on one edge, free on three. Where the clamped edge meets a free one the solution carries a corner singularity, and the frequency converges at about h² whatever the degree — measured 1.82 and 2.00 for p = 2 and 3.",
      "There is no closed form, so only the successive differences are drawn: they are what multilevel Monte Carlo measures too.",
    ],
    target: "edges",
    patch: { edges: "CFFF" },
    watch: "Two nearly parallel lines on the right: raising p hardly helps.",
    dwell: { ms: 5000 },
  },
  {
    id: "h-end",
    title: "Three rates",
    body: [
      "A hierarchy is described by how fast the quantity settles (α), how fast its corrections' variance falls once the input is random (β), and how fast a solve's cost grows (γ). The bias–variance tour adds the randomness; the multilevel tour spends all three.",
    ],
    target: null,
  },
];

const BIAS_VARIANCE: readonly TourStep[] = [
  {
    id: "bv-intro",
    title: "Two errors",
    body: [
      "Plain Monte Carlo on a coarse level: quadratic splines on 8 elements. Its estimate carries two errors. The sampling error σ/√N falls with N. The bias — the mesh's own error in E[Q] — does not.",
      "MSE = bias² + σ²/N. Bottom left plots both, relative to E[Q].",
    ],
    target: null,
    patch: { ...TOUR_BASE, view: "montecarlo", p: 2, mcLevel: 1, mcSamples: 100000 },
    watch: "The blue sampling error falls at slope −½, crosses the flat orange bias line early, and carries on below it — while the dashed total error stops falling.",
    dwell: { samples: 20000 },
  },
  {
    id: "bv-indicator",
    title: "Measuring a bias without the answer",
    body: [
      "The true E[Q] is unknown, so the bias cannot be measured directly. But every sample is also solved on the parent mesh, 4 elements, from the same ω. The mean of the difference, E[Q_ℓ − Q_ℓ₋₁], is how far one level moves from the last, and bounds the bias when the error falls geometrically.",
    ],
    target: "mcLevel",
    watch: "The readout: E[Y] ± its interval, resolved from zero.",
    dwell: { samples: 40000 },
  },
  {
    id: "bv-nstar",
    title: "The crossing, N*",
    body: [
      "Where the two lines meet, N* = σ²/bias², the errors are equal: about 70 samples here. Before it, more samples help. Past it, the mesh dominates: at N = 10⁵ the root-mean-square error is the bias alone to within a fraction of a per cent.",
    ],
    target: "mcSamples",
    watch: "The dashed vertical line, and the readout's verdict on which side of it the run is.",
    dwell: { samples: 100000 },
  },
  {
    id: "bv-refine",
    title: "Refine the mesh",
    body: [
      "Two levels finer: 32 elements. The bias indicator falls about twelvefold, from 1.4·10⁻² to 1.2·10⁻³ of E[Q] — asymptotically α = 2 for ω₁ at p = 2, a factor 16 over two levels, a little less on meshes this coarse. N* = σ²/bias² climbs with its square, from about 70 to about 10⁴.",
      "Refining is how to beat the bias. The price is that every one of the N samples is now solved on the finer mesh.",
    ],
    target: "mcLevel",
    patch: { mcLevel: 3 },
    watch: "The orange bias line dropping, and N* moving off to the right.",
    dwell: { samples: 20000 },
  },
  {
    id: "bv-correction",
    title: "The correction is quiet",
    body: [
      "Look at the readout's line V[Y]/V[Q]. The correction Y = Q_ℓ − Q_ℓ₋₁ of two solves from the same ω is a tiny fraction as variable as Q itself: the fine and coarse beams wobble together.",
      "A quantity that varies that little needs few samples to pin down. Estimate E[Q] on a coarse mesh with many cheap samples, and each correction on the finer meshes with few: that is multilevel Monte Carlo.",
    ],
    target: null,
    watch: "V[Y]/V[Q] in the readout: orders of magnitude below 1.",
    dwell: { ms: 4000 },
  },
  {
    id: "bv-end",
    title: "Bias, variance, and what to pay for",
    body: [
      "To reach a root-mean-square error ε, plain Monte Carlo needs a mesh whose bias is below ε and about σ²/ε² samples on it. Both get dearer as ε falls, and they multiply. The multilevel tour shows them added instead.",
    ],
    target: null,
  },
];

const MULTILEVEL: readonly TourStep[] = [
  {
    id: "ml-intro",
    title: "A telescoping sum",
    body: [
      "E[Q_L] = E[Q₀] + Σ_ℓ E[Q_ℓ − Q_ℓ₋₁]. Each term is estimated by its own Monte Carlo average — independent across levels, but with the fine and coarse solve of a sample sharing one ω. Q₀ is cheap and variable; the corrections are dear and quiet.",
      "The run first surveys every level with 2,000 samples, as Giles' mlmc_test does, then runs the adaptive algorithm at five tolerances.",
    ],
    target: null,
    patch: { ...TOUR_BASE, view: "mlmc" },
    watch: "The six panels filling as the survey's levels come in.",
    dwell: { samples: 12000 },
  },
  {
    id: "ml-rates",
    title: "α and β",
    body: [
      "Top row: variance and |mean| per level on log₂ axes, orange squares for Q_ℓ (what standard Monte Carlo sees) and blue circles for the correction Y_ℓ = Q_ℓ − Q_ℓ₋₁ (what MLMC sees). The blue lines fall with slopes −β and −α, given at each panel's top right: the variance each level adds, and the bias it leaves. Measured here: α ≈ 3.5, β ≈ 7.7, with V[Y_ℓ]/V[Q_ℓ] from 7.5·10⁻⁴ at ℓ = 1 to about 5·10⁻¹³ at ℓ = 5.",
      "The cost per sample grows as 2^(γℓ) with γ ≈ 1 for a beam. β > γ is the good case: the variance falls faster than the cost rises.",
      "Middle row: Giles' two checks. The consistency check stays below 1 unless the fine and coarse solves of a sample fail to share their ω; a large kurtosis would mean V[Y_ℓ] rests on a few rare samples.",
    ],
    target: null,
    watch: "Top left: the orange line flat, the blue one diving.",
    dwell: { ms: 4000 },
  },
  {
    id: "ml-allocation",
    title: "Samples where they are cheap",
    body: [
      "Giles' algorithm sets N_ℓ ∝ √(V_ℓ/C_ℓ) — the allocation that minimises cost for a variance budget — and adds a level while the bias estimate is above √θ ε. At ε = 3·10⁻⁴ it settles on four levels with N_ℓ = 236,328 · 4,152 · 585 · 58.",
      "The run is a pure function of the seed: those are the numbers you will get.",
    ],
    target: "mlEps",
    watch: "Bottom left: the blue bars, steeply falling with level, against the one orange bar standard MC would need on the finest.",
    dwell: { samples: 240000 },
  },
  {
    id: "ml-cost",
    title: "The complexity gain",
    body: [
      "Bottom right: cost against ε. The dashed lines are both costs predicted from the survey; the markers are each tolerance's run, standard MC's hollow because it is computed, not run. With β > γ the theorem says multilevel Monte Carlo costs O(ε⁻²) — as if the mesh were free — while standard Monte Carlo on the mesh the same ε needs costs ε^(−2−γ/α). The panel's top right gives the slopes the markers make.",
      "At ε = 3·10⁻⁴ MLMC takes 241,123 samples, and standard MC would need 202,402, every one on level 3. So MLMC takes more samples, yet costs 5.3× less: 98% of its samples are on level 0, where one costs 6.6× less than a level-3 solve. Standard MC's count is computed from the survey's V[Q₃], not run.",
      "Modest: cubic splines converge so fast that only three or four levels are ever needed.",
    ],
    target: null,
    watch: "The readout's saving column, growing as ε falls.",
    dwell: { ms: 4000 },
  },
  {
    id: "ml-live",
    title: "Inside the run, live",
    body: [
      "The live view shows the same run from the solver's side: one sample, its random stiffness, and its solutions on level ℓ and on level ℓ − 1 from the same ω, with every level's mesh. It steps through the level's samples every few seconds. Each pair is re-solved here from (seed, level, i) alone, and the readout checks it against the workers' numbers — they agree to the last bit.",
      "Below, the MLMC view's top row again, from every sample so far and redrawn as batches arrive, and N_ℓ: the samples in hand on each level against what ε = 3·10⁻⁴ asks for, 236,328 · 4,152 · 585 · 58. The dotted line marks the level shown: this pair is one of the 585 on level 2. \"run again from zero\" in the pane throws the run away to watch it arrive.",
    ],
    target: "liveLevel",
    patch: { view: "live", structure: "beam" },
    watch: "The two curves: nearly one, which is why their difference is so quiet.",
    dwell: { ms: 6000 },
  },
  {
    id: "ml-plate",
    title: "Where it pays: the plate",
    body: [
      "On the plate a solve costs h⁻⁴ (γ ≈ 4), so plain Monte Carlo's fine-mesh price is steep. Measured on the simply supported plate: savings of 1.0×, 3.8×, 13× and 38× as ε falls from 4.8·10⁻³ to 6·10⁻⁴.",
      "Each sample here is a 2D solve on up to 32 × 32 elements, so this run takes longer than the beam's.",
    ],
    target: "structure",
    patch: { view: "mlmc", structure: "plate", edges: "SSSS" },
    watch: "Bottom right: the gap between standard MC and MLMC widening as ε falls, and the readout's saving column filling in.",
    dwell: { samples: 30000 },
  },
  {
    id: "ml-end",
    title: "The whole argument",
    body: [
      { items: [
        "Monte Carlo: error σ/√N, whatever the dimension of the input.",
        "Discretisation: bias h^α, at a cost h^(−γ) per sample.",
        "Together, plain Monte Carlo pays both at once.",
        "Multilevel Monte Carlo spends its samples where they are cheap, on corrections that are quiet, and pays ε⁻² when β > γ.",
      ] },
    ],
    target: null,
  },
];

export const TOURS = {
  "Monte Carlo and √N": MONTE_CARLO,
  "discretisation hierarchies": HIERARCHY,
  "the bias–variance split": BIAS_VARIANCE,
  "multilevel Monte Carlo": MULTILEVEL,
} as const satisfies Record<string, readonly TourStep[]>;

export type TourName = keyof typeof TOURS;
export const TOUR_NAMES = Object.keys(TOURS) as TourName[];
