# Monte Carlo for vibrating structures

A teaching app for Monte Carlo error estimation of functions and functionals of
vibrating beams and plates, discretised with arbitrary-order spline finite
elements: plain Monte Carlo, hierarchies of discretisation error, and multilevel
Monte Carlo. The staged build is laid out in `MC_PLAN.md`; stages 0–6
are in — the deterministic foundation, the random input, plain Monte Carlo on
any level of the hierarchy, and multilevel Monte Carlo across it, run in Web
Workers.

## Layout

    src/
      main.ts        entry point: one figure of panels, one readout, one pane
      quad.ts        Gauss–Legendre of any order
      spline.ts      B-spline spaces of any degree p and continuity C^k; flat basis tables
      band.ts        symmetric banded storage and Cholesky (half-bandwidth p); complex symmetric LDLᵀ
      eig.ts         dense tred2/tql2; subspace iteration for the lowest modes
      beam/
        beam.ts      Euler–Bernoulli assembly, supports, static and modal solves
        exact.ts     closed forms: βL roots, Macaulay statics, the damped modal series
      hierarchy.ts   the level ladder, its rates α and γ, measured round-off; every QoI
      random/
        philox.ts    counter-based Philox4x32-10: normals addressed by (seed, sample, channel, stream, j)
        kl.ts        Karhunen–Loève by Nyström: exponential, Matérn-3/2, squared-exponential
        field.ts     a sample's stiffness, mass and load at any set of points
      mc/
        sampler.ts   one sample on a level and its parent, same ω — pure, tested directly
        stats.ts     Welford moments merged in sample order, of any prefix; histogram
        mlmc.ts      Giles' adaptive MLMC as a state machine; the survey; a sweep of tolerances
        worker.ts    a Web Worker around the sampler
        pool.ts      the worker pool: runs of one or more streams, kept across views
      ui/
        units.ts     the one place the beam acquires metres and hertz
        plot.ts      a 2-D canvas chart: log axes, legend, hover
        figure.ts    the grid of panels a view draws all its graphs into at once
        controls.ts  Tweakpane pane; owns no solver state
        views/       one per stage: basis, beam, convergence, spectrum, field, montecarlo,
                     mlmc; workers.ts holds the pool they share
    tests/           npm test: quadrature, splines, linear algebra, beam, hierarchy,
                     random inputs, forced response, Monte Carlo, multilevel Monte Carlo

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

## The forced response

The fifth quantity of interest is the steady amplitude |w(x_q)| when the load
is applied harmonically at Ω, with Rayleigh damping C = aM + bK:

    (K − Ω²M + iΩC) u = F.

Ω is set as a fraction of ω₁ of the uniform beam on the same supports, and a, b
so that that beam's first two modes have the damping ratio ζ — the textbook
fit, ζ_n = a/(2ω_n) + bω_n/2, which the tests check against the computed modes.
The defaults are Ω = 0.8 ω₁ and ζ = 0.02. Damping is required, not optional:
undamped, the amplitude is infinite wherever a sample's ω₁ lands on Ω, and so
is the variance (MC_PLAN.md's second pitfall).

The system is complex symmetric, not Hermitian. `complexSolve` factors it as
banded LDLᵀ in complex arithmetic without pivoting: same band, same O(n p²) as
the real Cholesky. That is safe because the imaginary part ΩC is positive
definite. −iA then has a positive definite Hermitian part, which is the
classical condition for Gaussian elimination to need no pivots.

What the tests hold it to:

- **Modal superposition.** Rayleigh damping is diagonal in the undamped modes,
  so on the discrete problem u = Σ φ_n(φ_nᵀF)/(λ_n − Ω² + iΩ(a + bλ_n)) is an
  identity. The banded solve matches it to about 10⁻¹¹ on cantilevers and propped
  beams with non-uniform stiffness and mass, above resonance.
- **A closed form.** For the uniform pinned–pinned beam the same sum over the
  exact modes √2 sin nπx is the exact answer (`pinnedResponse`). It reduces to
  5/384 and 1/48 as Ω → 0. At resonance the first mode amplifies its static
  share by 1/(2ζ).
- **Rates.** Against that series, midspan amplitude under a distributed load
  converges at α = 2.00, 4.00 for p = 2, 3 (as 2(p − 1)). For p = 4 and 5 it
  converges at about 6.2 and 6.3, not 6 and 8. A point value has no single rate,
  and theory gives none. Under a point load at midspan: 2.0, 4.0, 3.0, 2.9.

Near resonance the measured round-off probe can sit an order of magnitude
below the round-off a solve actually makes. So rate fits now also stop where
an error curve turns upward (`fitRate`), whatever the probe says.

As a random quantity it is the hard case. At the defaults (Matérn-3/2,
ℓ = 0.2, σ = 0.3, cantilever), ω₁ varies by about 12%, so some samples sit near
Ω: the coefficient of variation is about 110%, the mean is twice the mean
beam's amplitude, and the corrections have a kurtosis of 24–49.

## Multilevel Monte Carlo

Level ℓ's samples come from Philox stream ℓ: the fourth counter word, unused
until now. Within a level, the fine and coarse solves read the same ω. Across
levels the estimators are independent. The telescoping estimator

    E[Q_L] ≈ Σ_ℓ (1/N_ℓ) Σ_i Y_ℓ⁽ⁱ⁾,   Y₀ = Q₀,  Y_ℓ = Q_ℓ − Q_ℓ₋₁,

has MSE = bias² + Σ V_ℓ/N_ℓ. `Mlmc` is Giles' adaptive algorithm (the 2015
`mlmc.m`), written as a state machine. `wants()` names the cumulative samples
each level needs. The caller produces them, in workers, in any order. Then
`update()` reads the moments of exactly those prefixes and decides the next
round:

- it fits α and β on levels 1 … L;
- it sets N_ℓ = ⌈√(V_ℓ/C_ℓ) Σ√(V_k C_k) / ((1 − θ)ε²)⌉;
- it adds a level while the bias estimate exceeds √θ ε.

θ = ¼, as in Giles' code. A run therefore depends only on its seed, and one
store of samples per level serves a whole sweep of tolerances. Each tolerance
reads its own prefixes and lands exactly where it would have alone; the tests
check this bit for bit. `Accumulator.prefix(n)` makes that cheap: Welford
snapshots every 1024 samples, so any prefix costs at most 1024 updates.

Three departures from Giles' code:

- **Bias test.** It reads corrections only, levels ≥ 1. His reads |E[Q₀]| when
  L = 2, the quantity itself rather than a correction, and so always asks for a
  third level.
- **Cost.** It is the model C_ℓ = dofs_ℓ + dofs_ℓ₋₁ (free coefficients of both
  solves), not wall time. Wall time differs from run to run and would break
  reproducibility. Worker time per sample is measured and shown beside it.
- **Budget.** A tolerance that would need more than 2·10⁶ samples on one level
  stops as "over budget" and reports what it would have needed.

The view follows Giles' `mlmc_test`. A survey takes 2000 samples on every
level; the survey's level 0 also fixes the scale ε is relative to. Five
tolerances, ε_min·{16, 8, 4, 2, 1}, run at once. Four plots: mean and variance
against level, N_ℓ against level, and ε²·cost against ε for MLMC and for plain
Monte Carlo on the finest level each tolerance needed. The survey table also
gives each level's kurtosis and the consistency check |E[Y_ℓ] + E[Q_ℓ₋₁] −
E[Q_ℓ]| / 3σ, which stays below 1 unless the coupling is broken.

What the tests hold it to:

- **β against stage 3.** A perfectly correlated field factors out of every
  solve, so Y_ℓ = ΔQ_ℓ · e^{−σξ+σ²/2} for compliance. That gives
  V[Y_ℓ] = ΔQ_ℓ² e^{2σ²}(e^{σ²} − 1) exactly, with ΔQ_ℓ the stage 3 successive
  difference. The survey matches it level by level within four standard errors,
  and fits β = 2α of the deterministic hierarchy to within 0.15 for p = 2, 3.
- **Unbiased telescoping.** The MLMC estimate lands within four standard
  errors of E[Q_L] in closed form on the finest discrete mesh. The sampling
  variance is within its (1 − θ)ε² budget.
- **RMSE.** Over 32 independent runs at ε = 4·10⁻³ (relative), the RMSE
  against the exact continuum answer is 1.08 ε; the test allows up to 1.3 ε for
  χ² scatter. Every run added a third level on its own.
- **Complexity.** ε²·cost = 1.06, 1.01, 1.02 at ε = 4, 2, 1 ·10⁻³ (p = 2,
  β = 4 > γ = 1): O(ε⁻²), flat as the theorem says.

Measured in the browser (dev build, 8 workers), at the defaults: cubic C²
cantilever, ω₁, ne₀ = 4, six levels, Matérn-3/2 with ℓ = 0.2 and σ = 0.3.

- **Rates.** The survey gives α ≈ 3.5, β ≈ 7.7 and γ = 0.95.
- **Variance ratio.** V[Y_ℓ]/V[Q_ℓ] falls from 7.5·10⁻⁴ at ℓ = 1 to 5·10⁻¹³ at
  ℓ = 5.
- **The sweep.** It takes 248,000 samples and 7.8 s. At ε = 3·10⁻⁴, MLMC uses
  L = 3 with N_ℓ = 236,328 · 4,152 · 585 · 58 and costs 5.3 times less than
  plain Monte Carlo on level 3. The coarser tolerances save 1.1–2.8 times.

The savings are modest because a cubic spline hierarchy converges so fast that
the bias needs only three or four levels. The gain is bounded by about C_L/C₀,
which is small for a short hierarchy. Quadratic splines, a rough (exponential)
field, or the forced response lengthen it. The forced response under the
defaults needs L = 3 already at ε = 2.4·10⁻³ and saves 3.7 times. Tolerances
below 10⁻³ stop over budget: V[Q₀]/E[Q]² ≈ 1.2 alone would need ~10⁷ samples.

## Toward WebGPU

The plan is for sampling to move to the GPU. Today it runs in Web Workers. What
is already shaped for the move:

- basis tables (`tabulate`) are flat arrays indexed (element, point,
  derivative, local function), and coefficients arrive sampled at quadrature
  points — the form a random-field sample takes, and a storage buffer holds;
- matrices are flat lower bands, so the natural kernel is one invocation per
  sample running a banded Cholesky of half-bandwidth p;
- the random field reaches a level as one table √λ_j φ_j(x_q), and sample i's
  normals are Philox of (seed, i, channel, stream, j): no state passes between
  samples;
- a multilevel run is a set of independent streams of independent samples,
  produced in any order — one dispatch per level would do.

What is not settled is precision. The round-off section above is measured in
f64; in f32 the same conditioning (λ_max/λ₁ of 10⁶–10⁹ already at 64 elements,
depending on the supports) leaves nothing. A GPU path will need mixed precision with refinement in
emulated double, or double-f32 arithmetic throughout — to be decided, with this
CPU reference as the parity target, when the throughput is needed. Stages 5 and
6 have not needed it: the default MLMC sweep takes 248,000 samples in 8 s on
workers. The 2D plate of stage 7 may.
