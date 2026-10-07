# Monte Carlo for vibrating structures

A teaching app for Monte Carlo error estimation of functions and functionals of
vibrating beams and plates, discretised with arbitrary-order spline finite
elements: plain Monte Carlo, hierarchies of discretisation error, and multilevel
Monte Carlo. The staged build is laid out in [`PLAN.md`](PLAN.md); stages 0–3
are in — the deterministic foundation every Monte Carlo stage samples.

## Layout

    src/
      main.ts        entry point: one canvas, one readout, one pane
      quad.ts        Gauss–Legendre of any order
      spline.ts      B-spline spaces of any degree p and continuity C^k; flat basis tables
      band.ts        symmetric banded storage and Cholesky (half-bandwidth p)
      eig.ts         dense tred2/tql2; subspace iteration for the lowest modes
      beam/
        beam.ts      Euler–Bernoulli assembly, supports, static and modal solves
        exact.ts     closed forms: βL roots, Macaulay statics
      hierarchy.ts   the level ladder, its rates α and γ, and measured round-off
      ui/
        units.ts     the one place the beam acquires metres and hertz
        plot.ts      a 2-D canvas chart: log axes, legend, hover
        controls.ts  Tweakpane pane; owns no solver state
        views/       one per stage: basis, beam, convergence, spectrum
    tests/           npm test: quadrature, splines, linear algebra, beam, hierarchy

This directory is the development workshop, mirroring `apps/mantle`; it is
excluded from the Jekyll build. `npm run build` emits the static bundle into
`../../assets/mc-sim/`, which the site serves and embeds via an iframe in
`/mc-sim.html`.

## Develop

    npm install
    npm run dev        # hot-reloading dev server
    npm test
    npm run build      # type-check + bundle into ../../assets/mc-sim/

## The discretisation

A space is a degree p, a continuity k ≤ p − 1 across interior knots and a
uniform mesh, realised as an open knot vector whose interior knots repeat
p − k times. That one multiplicity spans h-, p- and k-refinement: p = 3, k = 1
is the Hermite cubic beam element, p = 3, k = 2 the cubic spline space on the
same mesh. The Euler–Bernoulli energy needs w″, so a beam needs C¹ at least;
the app refuses C⁰ for a beam and says why.

Supports are imposed strongly on an open knot vector: w(0) is c₀ alone and
w′(0) involves only c₀, c₁, so a clamp removes the first two coefficients and a
pin the first. The free system is therefore a contiguous block of the full
banded matrices — still half-bandwidth p, no renumbering — and a solve costs
O(n p²).

The solver is nondimensional (L, EI₀, ρA₀ all 1). `ui/units.ts` converts for
display only, with a reference beam that is printed beside every dimensional
number; the toggle never re-solves.

## What the tests hold it to

- **Exactness.** Under a uniform load the solution is a quartic, and every
  space of degree ≥ 4 reproduces it to round-off, for all four support pairs
  and for C¹ as well as maximal continuity. A point force on a knot gives a
  piecewise cubic, reproduced exactly by cubic C¹ and C² spaces.
- **Rates.** A manufactured cantilever — eˣ stiffness, w = 1 − cos x, with the
  end moment and shear it implies — converges at h^(p−1) in the H² seminorm
  and h^min(p+1, 2(p−1)) in L² (Aubin–Nitsche gains only 2(p − 1) at p = 2).
  First eigenvalues converge at h^(2(p−1)) for every support pair, and every
  discrete eigenvalue lies above the exact one.
- **Closed forms.** β_nL against Blevins' tables; statics against the textbook
  deflections; the 1 m, 20 mm square steel cantilever at 16.71 Hz.

Measured in the hierarchy view, maximal continuity, ne = 4 · 2^ℓ: clamped–
clamped ω₁ converges at α = 2.00, 3.98, 6.16, 8.40 for p = 2…5 (theory 2, 4,
6, 8); the cantilever at 2.00, 4.01, 5.94, with p = 5 already at round-off by
its third level.

## Round-off is a level of the hierarchy

The stiffness of a fourth-order problem is conditioned like h⁻⁴, and λ₁ of a
cantilever is small beside its largest discrete eigenvalue (λ_max/λ₁ ≈ 10⁹ at
p = 3 on 64 elements). Any f64 solve — the subspace iteration here and the dense
QL alike — therefore leaves λ₁ with an error near 10⁻¹⁰ that *grows* as the mesh
refines. The hierarchy measures it rather than assuming it: each level is
re-solved with K, M and F perturbed entrywise by (p + 1)ε (a backward-error
model of what Cholesky does anyway), and the change in Q is that level's
round-off. It is drawn in the convergence view, tabulated, and rates are fitted
only to the finest three levels standing ten times clear of it. Where an error
curve meets it, the level above is the finest worth paying for — a number
multilevel Monte Carlo will need.

Two things were learned getting the eigensolver there. A Ritz value converges
as the square of its vector, so stopping on the values alone hands back vectors
good to √tol; `lowestModes` also waits for the vectors, or for them to stop
improving. And the projected stiffness is formed as Yᵀ(MX) (Bathe's form, K Y =
M X by construction) rather than Yᵀ K Y, which keeps K's cancellation out of the
Rayleigh quotient; the cantilever, which had needed up to 17 sweeps, now takes
4–6 like everything else.

## Toward WebGPU

The plan is for sampling to move to the GPU; nothing here does yet. What is
already shaped for it:

- basis tables (`tabulate`) are flat arrays indexed (element, point,
  derivative, local function), and coefficients arrive sampled at quadrature
  points — the form a random-field sample takes, and a storage buffer holds;
- matrices are flat lower bands, so the natural kernel is one invocation per
  sample running a banded Cholesky of half-bandwidth p.

What is not settled is precision. The round-off section above is measured in
f64; in f32 the same conditioning (λ_max/λ₁ of 10⁶–10⁹ already at 64 elements,
depending on the supports) leaves nothing. A GPU path will need mixed precision with refinement in
emulated double, or double-f32 arithmetic throughout — to be decided, with this
CPU reference as the parity target, when stage 5 needs the throughput.
