# Monte Carlo for vibrating structures

A teaching app for Monte Carlo error estimation of functions and functionals of
vibrating beams and plates, discretised with arbitrary-order spline finite
elements: plain Monte Carlo, hierarchies of discretisation error, and multilevel
Monte Carlo. The staged build is laid out in `MC_PLAN.md`, and all of it —
stages 0–9 — is in — the deterministic foundation, the random input, plain Monte Carlo on
any level of the hierarchy, and multilevel Monte Carlo across it, run in Web
Workers — for the beam and for the Kirchhoff plate, with a view that shows the
solver and the statistics of one run side by side, and four guided tours. The
formulation is written up on the site page, `/mc-sim.html`.

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
      plate/
        plate.ts     Kirchhoff assembly on tensor-product splines, edges C/S/F, banded solves
        exact.ts     Navier (SSSS) and Lévy (CSCS) closed forms; Leissa's frequencies
        qoi.ts       the plate's quantities of interest and its hierarchy
      hierarchy.ts   the level ladder (any structure), its rates α and γ, measured round-off
      random/
        philox.ts    counter-based Philox4x32-10: normals addressed by (seed, sample, channel, stream, j)
        kl.ts        Karhunen–Loève by Nyström: exponential, Matérn-3/2, squared-exponential
        field.ts     a sample's stiffness, mass and load at any set of points
        field2d.ts   the plate's separable field: products of two 1D expansions, on a grid
      mc/
        sampler.ts   one sample on a level and its parent, same ω — pure, tested directly
        stats.ts     Welford moments merged in sample order, of any prefix; histogram
        mlmc.ts      Giles' adaptive MLMC as a state machine; the survey; a sweep of tolerances
        worker.ts    a Web Worker around the sampler
        pool.ts      the worker pool: runs of one or more streams, kept across views
      ui/
        units.ts     the one place the beam acquires metres and hertz
        plot.ts      a 2-D canvas chart: log axes, legend, hover
        heatmap.ts   a field over a rectangle: colour scales, contours, colour bar, meshes
        figure.ts    the grid of panels a view draws all its graphs into at once
        controls.ts  Tweakpane pane; owns no solver state
        options.ts   what each control offers: lists and slider ranges
        visibility.ts which controls each view shows, as a pure function of State
        tour.ts      the guided-tour overlay (after the mantle app's)
        tours.ts     the four tours, as data
        views/       one per stage: basis, structures (beam.ts and plate.ts side by
                     side), convergence, field, montecarlo, mlmc, live;
                     structure.ts says what each takes from the beam or the plate;
                     workers.ts holds the pool they share
    tests/           npm test: quadrature, splines, linear algebra, beam, hierarchy,
                     random inputs, forced response, Monte Carlo, multilevel Monte Carlo,
                     the plate, the tours

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

Measured with `runHierarchy`, maximal continuity, ne = 4 · 2^ℓ: clamped–
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

The view follows Giles' `mlmc_test`. A survey takes a fixed number of samples
on every level: 10² by default, 2·10³ in Giles' code and in the tours. The
survey's level 0 also fixes the scale ε is relative to. The survey sets a floor
under every level's sample count, but no tolerance reads past its own N_ℓ. A
small survey therefore leaves each tolerance's N_ℓ alone, apart from the shift
in scale, and saves the survey's solves on the fine levels, the dearest ones.
It costs precision in the survey's statistics. At 10², V[Q_L], and so plain
MC's count, is good to about ±30% (95%). Five tolerances,
ε_min·{16, 8, 4, 2, 1}, run at once.

The figure is Giles' `mlmc_plot`: six panels, each headed with what it plots
and, at the top right, the rate it measures. Colour and marker are the
estimator in every panel: orange squares are standard Monte Carlo and what it
sees (Q_ℓ), blue circles are MLMC and what it sees (Y_ℓ). Dashed lines are
predictions; hollow marks are computed or still settling, not measured.
- (a) variance and (b) |mean| of Q_ℓ and of Y_ℓ = Q_ℓ − Q_ℓ₋₁ against level,
  from the survey, on log₂ axes. The blue lines' slopes are −β and −α.
- (c) the consistency check |E[Y_ℓ] + E[Q_ℓ₋₁] − E[Q_ℓ]| / 3σ, which stays
  below 1 unless the coupling is broken, and (d) the kurtosis of Y_ℓ.
- (e) samples per level for the "tolerance shown": MLMC's N_ℓ as bars, beside
  the bar standard Monte Carlo would need on level L. The other tolerances'
  N_ℓ are faint lines; γ is at the top right.
- (f) cost against ε. Dashed: both costs predicted from the survey across the
  sweep and a factor 2 either side, with L from the bias test and the optimal
  N_ℓ unrounded. Markers: each tolerance's MLMC cost and standard MC's for the
  same ε, with the slopes they make at the top right.

The readout is `mlmc_test`'s printout: the survey table, α, β and γ, and per
tolerance the estimate, both costs, the saving and N_ℓ.

Both estimators are held to the same budget: bias from the finest level L the
sweep needed, and variance (1 − θ)ε². Plain MC's count is the one it would need, N = V[Q_L]/((1 − θ)ε²), from the
survey's V[Q_L]. It is computed, not run. At the defaults with a 2·10³ survey,
at ε = 3·10⁻⁴:
- MLMC takes 241,123 samples and plain MC 202,402;
- MLMC is 5.3× cheaper all the same, because 98% of its samples are on level 0,
  where one costs 6.6× less than a level-3 solve.

With the default 10² survey, MLMC takes 237,580 samples and plain MC 170,794.
MLMC is 4.5× cheaper; the difference is the survey's V[Q₃], read from 100
samples rather than 2000.

On the simply supported plate, at ε = 6·10⁻⁴, MLMC also takes more samples
(1.13× as many) and is 38× cheaper.

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

Measured in the browser (dev build, 8 workers), at the defaults with a 2·10³
survey: cubic C² cantilever, ω₁, ne₀ = 4, six levels, Matérn-3/2 with ℓ = 0.2
and σ = 0.3.

- **Rates.** The survey gives α ≈ 3.5, β ≈ 7.7 and γ = 0.95.
- **Variance ratio.** V[Y_ℓ]/V[Q_ℓ] falls from 7.5·10⁻⁴ at ℓ = 1 to 5·10⁻¹³ at
  ℓ = 5.
- **The sweep.** It takes 248,000 samples and 7.8 s. At ε = 3·10⁻⁴, MLMC uses
  L = 3 with N_ℓ = 236,328 · 4,152 · 585 · 58 and costs 5.3 times less than
  plain Monte Carlo on level 3. The coarser tolerances save 1.1–2.8 times.
- **Default survey.** With 10² per level, the sweep takes 238,000 samples and
  36% less work, with N_ℓ = 232,855 · 4,091 · 577 · 57 at ε = 3·10⁻⁴ (sample
  counts and work from a headless replay of the same run).

The savings are modest because a cubic spline hierarchy converges so fast that
the bias needs only three or four levels. The gain is bounded by about C_L/C₀,
which is small for a short hierarchy. Quadratic splines, a rough (exponential)
field, or the forced response lengthen it. The forced response under the
defaults needs L = 3 already at ε = 2.4·10⁻³ and saves 3.7 times. Tolerances
below 10⁻³ stop over budget: V[Q₀]/E[Q]² ≈ 1.2 alone would need ~10⁷ samples.

## The Kirchhoff plate

The plate is a × 1 (lengths in units of the side along y), with bending
stiffness D and mass ρt relative to their reference values, and Poisson's
ratio ν. Its energy needs second derivatives, like the beam's, which is where
splines earn their place: a C¹ triangle is Argyris' 21-dof element, while a
tensor product of two C¹ spline spaces is just two beams' bases multiplied.
The assembly reads 1D tables in x and y and forms w_xx, w_yy, w_xy as products
of them at each Gauss point.

Each edge is clamped (C), simply supported (S) or free (F), named in Leissa's
order: x = 0, y = 0, x = a, y = 1. As on the beam, the essential conditions are
imposed strongly by removing coefficient rows (two for C, one for S), so the
free coefficients form a rectangle of the index grid. The natural ones come
from the weak form: zero moment on S and F, plus the Kirchhoff shear and corner
forces on F, all of them ν-dependent. The free unknowns are numbered x-fastest,
so the matrix is banded with half-width p·nfx + p and the same banded Cholesky
serves.

That costs O(N·bw²) = O(p² ne⁴): γ = 4 against the beam's 1. The hierarchy
fits γ to that work, dofs × bandwidth², and so does the MLMC cost model; the
tests measure γ = 3.8–4 for p = 2 between 8 and 64 elements a side. A free
plate (FFFF) has three rigid modes, so its modes are found by iterating on
K + σM and shifting back; it has no static solution, and the views say so.

The random field is separable, C₁(|Δx|) C₁(|Δy|) with C₁ the beam's kernel. Its
Karhunen–Loève pairs are then products of two 1D Nyström solves,
λ_ij = λ_i λ_j, φ_ij = φ_i(x) φ_j(y), and the M largest products are kept. On
a level's tensor grid of Gauss points a sample is A Ξ Bᵀ: one small product per
distinct j, never a table of NX·NY·M. For the squared exponential the product
is the isotropic kernel itself; for the others it is the standard separable
stand-in. Along x the expansion is the one on [0, 1] for length ℓ/a, read at
x/a. Products of two decaying spectra decay more slowly than either, so the
plate needs more terms than the beam for the same share of the variance: at the
defaults (Matérn-3/2, ℓ = 0.2) the beam's 24 terms keep 99.9%, while 24
products keep 88%.

What the tests hold it to:

- **Exactness.** The biquartic w = x²(1 − x)² y²(1 − y)² is clamped on every
  edge. Under the load ∇⁴w it is reproduced to 10⁻¹¹ by every space of degree ≥ 4,
  for C¹ as well as maximal continuity, at ν = 0 and 0.3.
- **Navier (SSSS).** The series gives 0.00406 q b⁴/D and 0.01160 P b²/D at the
  centre of the square, Timoshenko's coefficients. Measured in the convergence
  view on ne = 4…32, ω₁ converges at α = 2.01, 4.03, 5.96, 8.35 for p = 2…5
  (theory 2(p − 1)). The compliance converges at the same rate, and the
  ‖w‖ field against the series on a Gauss grid at α = 4.06 for p = 3.
- **Lévy (CSCS).** Clamped on x = 0, a and simply supported on y = 0, 1, the
  frequencies come from two transcendental equations, solved by bisection.
  Quintic splines on 16 × 16 match them to 10⁻⁷ at aspect ratios 1 and 2. The
  square's first is ω̂ = 28.95085.
- **Leissa.** Splines reproduce every digit he gives for CCCC (35.985, 73.394,
  108.22) and for the free plate at ν = 0.3 (13.468, 19.596, 24.270). For CFFF they
  converge from above to 3.4710, 8.506 and 21.284 (p = 4 on 48² and 64² agree to
  those digits). That is 0.6% below his 3.4917; a conforming Rayleigh–Ritz
  method can only overestimate, so his values are not the limit.
- **Corners.** Where a clamped or free edge meets another, the solution has a
  corner singularity that caps the rate. The successive differences of ω₁ fall
  at 2.16, 3.68, 4.91, 6.11 for CCCC and at 1.82, 2.00, 2.15, 2.26 for CFFF
  (p = 2…5). The cantilever plate converges at about h² whatever the degree,
  so its hierarchy is long, which is the case MLMC is for.
- **Forced response.** The banded complex solve equals the modal sum over the
  discrete modes to 10⁻¹⁰, for CSCS with non-uniform D and μ. On SSSS it converges
  to the Navier modal series at rate > 3.5 for p = 3, within 10⁻⁵ by 16 × 16.
- **The field.** With 600 products the Matérn covariance between grid points is
  reproduced to 10⁻⁴. Points shared by two grids get the same g, to 10⁻¹³, from the
  same ξ. Along x on a 2 × 1 plate the squared-exponential correlation is that of
  length ℓ, to 10⁻⁶.
- **Monte Carlo.** With a perfectly correlated field, E[w] = w₀e^{σ²} and
  E[ω₁] = ω₁ₕe^{−σ²/8} on the discrete plate, within four standard errors over
  3,000 samples.

Measured in the browser (dev build, 8 workers), at the defaults with the
structure set to the plate: SSSS, cubic C², ω₁, Matérn-3/2 with ℓ = 0.2 and
σ = 0.3, M = 24.

- **Plain Monte Carlo.** On 16 × 16 with its 8 × 8 parent: 540 samples/s,
  14.6 ms per sample per worker. V[Y]/V[Q] = 6.5·10⁻⁷ and N* ≈ 3·10⁵.
- **Survey.** MLMC on ne₀ = 4, four levels to 32 × 32: α ≈ 4.5, β ≈ 9.6,
  γ = 3.6 in the cost model. V[Y_ℓ]/V[Q_ℓ] = 4.6·10⁻⁴, 6.7·10⁻⁷, 7.3·10⁻¹⁰.
- **Savings.** 1.0×, 3.8×, 13× and 38× over plain Monte Carlo on the finest
  level needed, at ε = 4.8, 2.4, 1.2, 0.6 ·10⁻³ (relative), each reaching
  L = 2. The steep growth with 1/ε is γ = 4 at work: plain Monte Carlo pays
  h⁻⁴ for every sample, MLMC for a few.

Worker time grows more slowly than the cost model says: 6, 9, 13, 75 ms per
sample by level, a measured γ of 1.5. At these sizes the eigensolver's sweeps
(N·bw per vector per sweep, nine vectors, about seven sweeps) cost more than the
single factorisation (N·bw²). So the time grows as h⁻³ until the band passes a
few hundred.
The model keeps the factorisation's h⁻⁴, the asymptotic cost, and, unlike wall
time, it is reproducible.

## Visualisation: the solver beside its statistics

Stage 8 adds motion and one view that ties the stages together.

- **Forced response, animated.** The beam and plate view can show the steady
  response to the load applied at Ω instead of a mode: Re(u e^{iΩt}) =
  Re u cos Ωt − Im u sin Ωt, inside its envelope ±|u|, beside the static
  deflection for scale. The readout gives the amplitude at the QoI point
  against the closed form where there is one (pinned–pinned beam, SSSS
  plate), the dynamic amplification and the phase lag. Damping puts the
  points of the structure out of phase with one another, so a forced plate's
  zero contour travels, where a mode's stands still.
- **Meshes.** A "mesh" switch draws the element boundaries over the beam (as
  ticks) and over the plate (as lines).
- **7 · MLMC, live.** The MLMC view's run, from the solver's side. Across the
  top is one sample of level ℓ: the stiffness it drew, and its solution on
  level ℓ and on level ℓ − 1 from the same ω — the coupled pair whose difference
  the estimator averages. The beam also shows a row of mesh ticks for every
  level of the hierarchy, with the pair's two in colour. The shape moves as
  the quantity does: the first mode swinging for ω₁, the forced response for
  the response, still for a static quantity. The view steps through the
  level's first 48 samples, 2.5 s each.

  Below are three of the MLMC view's panels, recomputed from every sample each
  level has so far, so they move as batches arrive (the MLMC view reads the
  survey's fixed prefix instead):
  - (a) variance and (b) |mean| of Q_ℓ and Y_ℓ against level;
  - samples per level: the samples in hand as a line, over the bars of what
    the "tolerance shown" asks for and of what standard MC would need, and
    the survey's floor.

  A dotted line marks the shown level in each. **run again from zero** throws
  the run away so it can be watched arriving. The samples are the same again,
  since a run is a function of its seed: the default beam sweep refills in
  about five seconds, the plate's in about a minute. The readout gives the
  shown pair's Q and Y, the tolerance's N_ℓ, and what standard MC would pay.

Each pair is re-solved on the main thread from (seed, level, i) alone through
`Sampler.inspect`, the code path the workers run, and checked against the Q
and Q_c they returned. The test suite holds this to bit equality for beam and
plate and every quantity, and the live view's readout says, for each pair it
shows, whether the two agree.

`mlmcSession` is the run's setup — open it in the pool, set its demand, build
the sweep and the survey — taken out of the MLMC view so the live view drives
the very same run. Switching between the two views loses no samples, and the
default beam sweep still lands on N_ℓ = 232,855 · 4,091 · 577 · 57 at
ε = 3·10⁻⁴ (236,328 · 4,152 · 585 · 58 with a 2·10³ survey).

## Guided tours

Stage 9's four tours open from the pane's first folder: Monte Carlo and the √N
rate; discretisation hierarchies (h, p and k refinement); the bias–variance
split; and multilevel Monte Carlo with its complexity gain. Each is a list of
steps (`ui/tours.ts`). A step names the control it is about, a patch of
settings over everything earlier steps set, what to watch for, and a dwell:
wall time, or a sample count the active run must reach.

The overlay (`ui/tour.ts`) is the mantle app's, cut down. The pane dims, with
four shades tiled around the one lit control, which stays live. The figure and
the readout are never dimmed. A dwell fills a bar and lights "next" rather than
advancing. Walking back replays every patch up to the step. The last card
offers the settings the reader had before the tour.

Every tour starts from the same base settings, seed 1 included, so the numbers
its cards quote are the ones the reader will see. Each was checked against the
running app: for instance N* ≈ 70 at level 1 of the quadratic cantilever,
rising to about 10⁴ two levels finer, where the bias indicator has fallen
twelvefold (α = 2 asymptotically, a little less on meshes this coarse).

`tests/tour.test.ts` checks, without a browser, what would otherwise only fail
part-way through a tour:
- every patch sets only fields `State` has, to a value its control offers
  (`ui/options.ts`) or within its slider's range;
- the control each step lights is visible in the state that step leaves
  (`ui/visibility.ts`, the pane's show/hide rules as a pure function);
- every mesh a step runs on is one the app will solve;
- every sample dwell can be met by its run.

## Toward WebGPU

The plan is for sampling to move to the GPU. Today it runs in Web Workers. What
is already shaped for the move:

- basis tables (`tabulate`) are flat arrays indexed (element, point,
  derivative, local function), and coefficients arrive sampled at quadrature
  points — the form a random-field sample takes, and a storage buffer holds;
- matrices are flat lower bands, so the natural kernel is one invocation per
  sample running a banded Cholesky of half-bandwidth p (p·nfx + p for a plate);
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
workers. The plate is where it starts to: a sample on 32 × 32 with its parent
takes ~70 ms of worker time against the beam's 0.4–3 ms, so the views stop a
plate at 32 × 32, and a 64 × 64 level (~0.5 s a sample) waits for a faster
solver.
