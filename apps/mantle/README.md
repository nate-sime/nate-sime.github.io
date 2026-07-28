# Mantle convection — WebGPU solver

GPU-accelerated thermal convection in a 2D spherical annulus or a Cartesian box.
The velocity is represented by a stream function `u = ∇×(ψ ẑ)` (pointwise
divergence-free) and obtained from a biharmonic Stokes solve discretised in
high-order (C²) splines; temperature is transported by semi-Lagrangian advection
with implicit diffusion. The write-up of the formulation lives in
`/mantle-convection.html` on the site.

## Layout

    src/
      main.ts        device bootstrap + frame loop
      geometry.ts    the two domains, as the metric that distinguishes them
      spline.ts      B-spline axes, tensor field, u = ∇×ψ
      linalg / quad / dft    dense + tridiagonal f64 kernels, Gauss, reference DFT
      solver/        the CPU reference — f64, and the parity target for the GPU
      gpu/
        wgsl.ts      every compute and render kernel, as source builders
        sim.ts       buffers, pipelines, bind groups, frame encoding
      ui/
        controls.ts  Tweakpane pane; owns no solver state
        equation.ts  the selected law, written out under the list
        nusselt.ts   the Nu time series: rolling buffer + axis scales
        nuplot.ts    that series, drawn into the corner panel
        dimensional.ts  the one place the run acquires physical units
    tests/           npm test: convergence, boundary conditions, GPU parity

This directory is the development workshop; it is excluded from the Jekyll
build. `npm run build` emits the static bundle into `../../assets/mantle/`,
which the site serves and embeds via an iframe.

The CPU solver is not superseded by the GPU one. It is where a scheme is worked
out in f64, and it is what the GPU pipeline is verified against.

## The two geometries

Everything is written on one non-periodic axis `r` crossed with one periodic axis
`φ`, and the whole difference between the annulus and a box is the metric
`ds² = dr² + h(r)² dφ²` — `h = r` in one, `h = 1` in the other. `geometry.ts`
carries `h`, `h′`, the period of φ, the conduction profile and the Nusselt
normalisation, and nothing downstream knows which domain it is in. So the box is
**the same discretisation checked on a second geometry**, not a parallel code
path: `u = ∇×(ψ ẑ)` is pointwise divergence-free in both, the dissipation form is
the same integral with a different Jacobian, and free-slip is natural to it
either way. The buoyancy load needed no change at all — the `1/h` of the
along-gravity velocity cancels the `h` of the area element.

The box closes left and right either way, and the choice is a control:

- **periodic** — what the radix-2 FFT gives natively, since diagonalising the
  transverse direction *is* the statement that φ wraps.
- **free-slip walls** — impermeable, stress-free and insulating at x = 0 and
  x = L, reached by **mirroring**: solve on a period of 2L and hold the state
  even about x = 0. An even T gives an odd ψ, hence `u_x = −ψ_z = 0` and
  `ψ_xx = 0` on both planes, and `∂ₓT = 0` there. That is the walled problem,
  not a stand-in for it, and nothing in [0, L] is constrained — the reflection
  only decides what the invisible half does.

The projection is free, which is the point. Symmetry about x = 0 is exactly
`Im T̂ = 0`, and the diffusion solve is already in mode space, so **dropping the
imaginary half there is the wall** — no extra pass, no mirror-indexing anywhere
in the pipeline, and the CPU and the GPU make the identical projection. Doing it
every step is what makes it a boundary condition rather than a property of the
initial condition: the odd component is annihilated each step, so it can neither
drift in from round-off nor grow out of a symmetry-breaking mode. ψ needs no
projection of its own — oddness is inherited from the load, exactly in f64 and to
f32 round-off on the GPU, and ψ is re-solved from a freshly projected T each step
so that round-off cannot accumulate.

The price is resolution: a walled box of a given width is solved on twice the
period, so at a fixed `na` it resolves half as finely as a periodic one. Periodic
is therefore the default.

Two checks pin the box, and between them they cover the linear and the
finite-amplitude problem. `Ra_c = 27π⁴/4 = 657.5` is the sharpest single
statement about the metric: a box of length `2√2` is one wavelength of the
critical mode, so a single seeded roll decays at Ra = 400 and saturates at
Ra = 1200. A stray `1/h` in the operator, a wrong Jacobian in the assembly or a
transverse wavenumber that still assumed a 2π period all move that number, and
none of them would be visible in a run that merely looked like convection. And
the **Blankenbach 1a benchmark** — unit square, free-slip, Ra = 10⁴, accepted
`Nu = 4.884409` — comes out at **4.8947 on a 33×64 grid, 0.2% high**, and 4.8734
(0.2% low) at 49×96. Four configurations were sampled and all landed inside
0.25%; that is agreement at usable resolutions and not a convergence study — Nu
converges only first order here, and the sampled runs varied dt as well as the
grid. It is run twice, once as a length-2 periodic domain read in halves and once
as the walled unit square the benchmark actually states.

The walls themselves are checked against code that was already verified rather
than against a claim about what a wall should look like: a walled box of width L
*is* a periodic box of width 2L whose state stays symmetric, so the two run
**step for step to round-off** (measured < 1e-12 over 200 steps) and report the
same Nu. Then the walls behaving as walls — an antisymmetric field written in is
gone in one step, and on the settled unit square `u_x` and `∂ₓT` at both walls
measure 1.4e-15 and 1.2e-14 against an interior |u_x| of 61.

**GPU parity is stated in two checkpoints, and they are different claims.** A
mis-ported metric term is systematic: it makes the two paths solve different
equations, so it shows on the first step at a size f32 round-off never reaches —
five steps under a tight bound is what catches it. Twenty-five steps is not that
claim and cannot be, because the box at Ra = 10⁴ convects at max|u| ≈ 56 (nine
times the annulus at the same Ra, whose layer is 0.45 deep against this one's 1),
so it turns over more than once and two trajectories seeded ε apart in f32
separate at the rate the flow mixes. Measured: 1.0e-6 at one step, 8.7e-6 at
five, 1.9e-3 at twenty-five against a field that has changed by O(1). The late
bound says they have not *diverged*. Nu is read at the early checkpoint for the
same reason — it is a boundary gradient, so it divides by dr and amplifies that
separation by ~30.

Depth runs 0 → 1 with the hot boundary at 0, as the mantle convection literature
states it — which is also the annulus' own convention, where `r_i` is the
core–mantle boundary. Because the box's unit length is the layer *depth* and the
annulus' is the outer *radius*, the two are not on the same dimensional clock;
`ui/dimensional.ts` scales each by its own, and the readout names which.

Geometry, and the box's length, are **rebuilds**. The metric is emitted into the
WGSL rather than branched on a uniform — a box would otherwise evaluate `1/r` at
z = 0 — and the length reaches the azimuthal knot vector, so the quadrature
tables, the discrete symbols and the per-mode inverses are all built against it.

## Develop

    npm install
    npm run dev        # hot-reloading dev server (Chrome/Edge: best WebGPU support)
    npm test           # includes a real GPU run via headless Dawn
    npm run build      # type-check + bundle into ../../assets/mantle/

`npm test` uses the `webgpu` package (Dawn as a node addon) so the parity suite
exercises actual compute shaders. Where no adapter is available — most CI
runners — those tests skip; everything they assert is a parity claim about code
the CPU suite covers regardless.

## The pipeline

The whole pipeline runs in WebGPU compute shaders: buoyancy assembly (a
quadrature *gather*, since a scatter would need atomics) → azimuthal FFT
(shared-memory Stockham, one workgroup per row) → per-mode radial solve (a
matvec against dense inverses factorised in f64 at init) → inverse FFT →
semi-Lagrangian/BFECC advection → implicit diffusion → render. Nothing crosses
back to the host in the frame loop; the diagnostics readout is an asynchronous
poll of a GPU-side reduction, off the frame's dependency chain.

The render pass can overlay **ψ isocontours, which are exactly the streamlines**
— `u = ∇×(ψ ẑ)` is tangent to level sets of ψ, so there is no particle tracing and
no second buffer, just one spline evaluation per pixel. Contour spacing follows
`max|ψ|` from a GPU reduction, so the density stays readable across decades of
Ra without the host ever seeing ψ.

Either **mesh** can be overlaid the same way, again as a distance field and not
as geometry: both discretisations are uniform in (r, φ), so a line family is the
distance to the nearest multiple of an element width, and the counts — `nr − 3`
by `na` elements for ψ, `gnr − 1` by `gna` cells for T — are read off the axes
and uploaded with the rest of the uniform block. They are offered separately
because they are two different meshes. Each family fades out as it approaches one
line per pixel, the same Nyquist policy the contours follow. Both overlays start
off.

Both **Nusselt numbers are plotted against time** in the bottom-left corner,
accumulated from the same asynchronous poll the readout uses — so the chart adds
nothing to the frame's dependency chain, and redraws when a sample lands rather
than every frame. It is one axis for both series deliberately: they are the same
quantity at the two boundaries — named *inner*/*outer* on the annulus and
*bottom*/*top* in a box — and the reading is whether they have met. The pair
converging is a genuine global heat balance (the two boundary fluxes are
independent reductions over independent rows); the two instantaneous numbers in
the readout show the balance but not whether the run has reached it, is
oscillating, or is still absorbing a change to Ra. Because they *do* coincide at
steady state to more digits than a pixel can hold, they are drawn nested — the
outer curve wide, the inner narrow over it — so agreement reads as a rim around a
core rather than as a series that failed to draw. Reseeding or a rebuild drops the
trace, since joining two runs with a line would draw a trajectory nothing
followed.

*Nu window* under **view** sets how much of the run is on screen, from the last
500 steps to all of it. The span is in **solver steps** and not samples: the poll
rate is the frame loop's, so a window of "the last 400 polls" would cover a
different stretch of the simulation at every playback speed and would shift under
the user when they moved that list. Narrowing rescales the y axis onto the visible
slice, which is the point of it — a settled band is a thousandth of the initial
transient's height, so on an axis the transient still sets it is a flat line.
Nothing is discarded: the buffer keeps its last 16 384 polls whatever is
displayed, so widening the window again brings the earlier history back rather
than starting over.

Time is reported **both nondimensionally and dimensionally**, as two rows under
the plot's x axis and on the readout's clock. The solver has no physical units in
it and should not acquire any — that is the whole point of the Boussinesq scaling,
where Ra is the only parameter — so the conversion is a display choice confined to
`ui/dimensional.ts` and its assumption is printed in the readout's header rather
than left implicit. One nondimensional time unit is one thermal diffusion time
across whichever length is 1 in code units, and **that is not the same length in
the two geometries**: the annulus scales by the outer radius `R_o` and the box by
its depth `d`, a factor of nearly five in the clock. Both are the same choice
about Earth stated twice — the radius ratio `r_i/r_o = 0.55` is the core–mantle
boundary against the surface, and `d = R_o − R_i = 2891 km` is that same layer.
**Expect
figures in Gyr and Tyr, and read them as a result.** `R_o²/κ ≈ 1.3×10¹² yr` —
diffusion across the mantle is orders of magnitude slower than the age of the
Earth, which is exactly why the real mantle convects; and at Ra = 2×10⁴ this model
is far more viscous than the mantle, whose Ra is 10⁷–10⁹, so it needs many
diffusion times to do what the Earth does in a fraction of one.

Controls (Tweakpane) are grouped by what they cost: Ra, contour count, line
width, the mesh overlay and the power-law index are uniform writes; dt
re-factorises the diffusion
operator in f64; the contrast re-inverts the preconditioner's radial blocks;
reseed re-solves; changing the geometry, the box length, the resolution or the
solver tier rebuilds every table and pipeline and says so. The viscosity folder writes the selected law out beneath
the list that selects it, and names the slider behind each symbol with its
current value — neither symbol is the number on its slider, since `γ =
ln(contrast)` and `n` acts only through the exponent `(1−n)/n`.
The resolution ladder runs ψ 48×128 to ψ 192×512, and playback speed spans one
step per 16 frames to 16 per frame — a coarse mesh is otherwise correct and far
too fast to watch. Slowing down throttles the frame loop rather than shrinking
dt, which would change the trajectory rather than the rate it is shown at.

Verified against the CPU reference: `max|ΔT| = 1.1e-6` over 25 lockstep steps,
`‖∇·u‖/|u| < 1e-5`, and `ψ ∝ Ra` exactly. Runs at 0.98 ms/step for ψ 96×256
with T on 193×256; the streamline overlay costs 0.008 ms at 1024². Startup takes
~2 s to factorise the radial operators and compile every pipeline — announced on
the canvas rather than hidden, and paid once.

**μ(T)** switches the solve from direct to matrix-free preconditioned CG,
with the *same* FFT radial kernel as the preconditioner — on the GPU literally
the same pipelines with other buffers bound, differing in one uniform. The
viscosity is Frank–Kamenetskii, `exp(−γ(T−½))`, centred so the geometric mean
stays 1 and the contrast slider does not quietly rescale the effective Ra. Each
frame warm-starts from the previous ψ, which is what lets a fixed budget work:
four iterations already reach the f32 floor, and the default is twelve.

**μ(T, ε̇)** adds a regularised power law on top, `n ≈ 3` for dislocation
creep, with a viscosity floor and ceiling. It is nearly free: the operator's
first pass already computes `∂_φA[ψ]` and `C[ψ]` at every quadrature point, and
those *are* `ε_rr` and `−2ε_rφ`, so the second invariant is that pass and one
`hypot`. Measured in one run at ψ 96×256: 0.24 ms/step for the direct solve,
3.18 ms for μ(T) at 12 CG iterations, and **3.09 ms for μ(T, ε̇) at n = 3 and the
same budget** — the rheology costs three dispatches against ~7 per iteration.

`n = 1` collapses the power law to the identity *exactly*, so the two variable
laws are one tier, one set of pipelines and one uniform apart; only entering or
leaving the Krylov path rebuilds. The rheology is time-lagged — μ is frozen for
the duration of each solve, which is what keeps the operator CG sees linear and
symmetric — and extra Picard sweeps re-lag against the ψ just computed. That
sweep count, not the iteration budget, is the accuracy knob at n > 1: measured
in the running loop, tripling the CG budget buys ~30% and one extra sweep buys
4×, which is also why a multigrid fallback for the Krylov path was measured out
rather than built.

Note that a residual is *not* a convergence diagnostic here: ψ is stored in f32,
so it sits ε from anything, and the operator amplifies that by κ ~ h⁻⁴. The
measured residual is flat in the iteration count and grows with resolution —
convergence is asserted in f64, in `tests/rheology.test.ts`.

The bundle is embedded in `/mantle-convection.html` on the site, which carries
the write-up of the formulation. Publishing runs from `.github/workflows/`, with
the repository's Pages source set to *GitHub Actions*.
