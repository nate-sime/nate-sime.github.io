# Monte Carlo for vibrating structures

A teaching app for Monte Carlo error estimation of functions and functionals of
vibrating beams and plates, discretised with arbitrary-order spline finite
elements: plain Monte Carlo, hierarchies of discretisation error, and multilevel
Monte Carlo. The staged build is laid out in `MC_PLAN.md`; stages 0–5
are in — the deterministic foundation, the random input, and plain Monte Carlo
on any level of the hierarchy, run in Web Workers.

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
      random/
        philox.ts    counter-based Philox4x32-10: normals addressed by (seed, sample, channel, j)
        kl.ts        Karhunen–Loève by Nyström: exponential, Matérn-3/2, squared-exponential
        field.ts     a sample's stiffness, mass and load at any set of points
      mc/
        sampler.ts   one sample on a level and its parent, same ω — pure, tested directly
        stats.ts     Welford moments merged in sample order; histogram
        worker.ts    a Web Worker around the sampler
        pool.ts      the worker pool and the current run
      ui/
        units.ts     the one place the beam acquires metres and hertz
        plot.ts      a 2-D canvas chart: log axes, legend, hover
        controls.ts  Tweakpane pane; owns no solver state
        views/       one per stage: basis, beam, convergence, spectrum, field, montecarlo
    tests/           npm test: quadrature, splines, linear algebra, beam, hierarchy,
                     random inputs, Monte Carlo

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

## The random input

Stiffness is lognormal about the deterministic section,
e(x, ω) = e₀(x) exp(σ g(x, ω) − ½σ² s_M(x)), with g a Karhunen–Loève field
truncated at M terms and s_M(x) = Σ λ_j φ_j(x)² its pointwise variance.
Subtracting s_M rather than 1 makes E[e(x)] = e₀(x) exactly, for every x and
every M. Optionally the mass follows the stiffness as a random section depth
would (I ∝ d³, A ∝ d, so μ ∝ e^{1/3}), and the load gets an independent Gaussian
part (a field for a distributed load, a magnitude for a point load).

The KL eigenpairs come from a Nyström discretisation on 384 composite Gauss
nodes, and the Nyström interpolant evaluates them anywhere. Each level tabulates
√λ_j φ_j once at its own quadrature points, and a sample is then a matrix–vector
product. Because the normals are addressed rather than drawn, sample ω is the
same function on every level, which is the coupling MLMC needs.

The normals come from Philox4x32-10 (Random123's generator), not from PCG or
xoshiro as MC_PLAN.md first said. It is counter-based: normal j of sample i on a
channel is a pure function of (seed, i, channel, j). Any worker can run any
sample, in any order, with no state to jump ahead. A WebGPU kernel needs exactly
that property. It reproduces Random123's known-answer vectors.

What the tests hold it to:

- **Nyström** matches the exponential kernel's exact eigenvalues (Ghanem &
  Spanos' transcendental roots) to about (j/N)² for mode j, measured 7·10⁻⁴ at
  j = 10 and 3·10⁻² at j = 64. The kink of e^{−r/ℓ} at r = 0 caps the rule at
  second order. The spectra decay as j⁻² (exponential) and j⁻⁴ (Matérn-3/2),
  and the modes reproduce the covariance.
- **Coupling.** The same sample index gives the same e, μ and q at a point
  shared by two point sets, to 10⁻¹³.
- **Moments.** E[e] = e₀ and the claimed lognormal quantiles hold to within
  four standard errors over 2·10⁴ samples.

## Plain Monte Carlo

Each sample is solved on level ℓ (ne₀·2^ℓ elements) and, with the same ω, on
level ℓ − 1. Batches from the workers are folded in sample order, so a run's
statistics after N samples are those of samples 0 … N − 1, bit for bit, whatever
the worker count or batch sizes. The view shows four plots: the histogram
(with Q of the mean beam, to show the Jensen gap), the running mean with its
CLT interval, the RMSE split as bias² + σ²/N with the crossover N* = σ²/E[Y]²,
and the mean deflection field with ±σ and ±2σ bands.

Against closed forms: a perfectly correlated field makes e = exp(σξ − σ²/2) a
single number per sample. K then scales by it exactly, so E[w] = w₀e^{σ²}, and
E[ω₁] = ω₁ₕe^{−σ²/8} (or e^{−σ²/9} when the mass follows the depth), all on the
discrete mesh. Monte Carlo hits each within four standard errors. Over 300
independent runs of N = 60, the 95% interval covers the truth between 88% and
99% of the time. A random load leaves the mean alone and adds exactly its own
variance.

Measured on the coupled levels: Matérn-3/2, ℓ = 0.3, σ = 0.5, cubic C²
cantilever tip deflection, 200 samples per level. V[Y_ℓ]/V[Q_ℓ] falls as
2.7·10⁻⁴, 7.5·10⁻⁶, 1.2·10⁻⁷, 3.9·10⁻¹⁰ for ne = 8 … 64. The variance MLMC
works with is tiny from the start and shrinks by 36–300 times per level. In the
browser (dev build, 8 workers): tip deflection at 16 elements and its parent
runs about 18,500 samples/s; ω₁, a subspace iteration per solve, about 2,400.
With the default inputs, ω₁ on 16 elements has V[Y]/V[Q] = 2.3·10⁻⁵ and
N* ≈ 3·10⁴. Past that many samples the mesh, not the sampling, sets the error.

## Toward WebGPU

The plan is for sampling to move to the GPU. Today it runs in Web Workers. What
is already shaped for the move:

- basis tables (`tabulate`) are flat arrays indexed (element, point,
  derivative, local function), and coefficients arrive sampled at quadrature
  points — the form a random-field sample takes, and a storage buffer holds;
- matrices are flat lower bands, so the natural kernel is one invocation per
  sample running a banded Cholesky of half-bandwidth p;
- the random field reaches a level as one table √λ_j φ_j(x_q), and sample i's
  normals are Philox of (seed, i, channel, j): no state passes between samples.

What is not settled is precision. The round-off section above is measured in
f64; in f32 the same conditioning (λ_max/λ₁ of 10⁶–10⁹ already at 64 elements,
depending on the supports) leaves nothing. A GPU path will need mixed precision with refinement in
emulated double, or double-f32 arithmetic throughout — to be decided, with this
CPU reference as the parity target, when stage 5 needs the throughput.
