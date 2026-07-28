---
layout: default
title: Mantle convection
---

Thermal convection in a two-dimensional spherical annulus — or a Cartesian box —
solved live in the browser. Incompressible Stokes flow is coupled to the heat
equation under the Boussinesq approximation, and every part of the solve —
assembly, the Stokes system, temperature transport, and the rendering — runs in
[WebGPU](https://www.w3.org/TR/webgpu/) compute shaders. Nothing returns to the
CPU while it is running.

The point of interest is not the physics, which is textbook, but the
discretisation: the velocity field is **pointwise divergence free, exactly and
by construction**, and the free-slip boundary condition is imposed in the form
the curved geometry actually demands rather than the one inherited from
Cartesian boxes.

<figure class="embed">
  <div id="webgpu-note" class="embed-note" hidden>
    This browser does not appear to support WebGPU, so the frame below will show
    a notice rather than the simulation. Chrome, Edge, or Safari&nbsp;18+ will
    run it.
  </div>
  <iframe src="/assets/mantle/"
          title="Mantle convection in a 2D spherical annulus or Cartesian box, simulated with WebGPU"
          loading="lazy"></iframe>
  <figcaption>
    Temperature, dark to bright: cold to hot. Gravity points from the hot
    boundary towards the cold one; both are free-slip. Under <em>domain</em>,
    the geometry can be switched from the spherical annulus to a Cartesian
    box — depth 0 to 1, hot below, with an adjustable length — which rebuilds
    the solver and takes a second or two. Use the rest of the pane
    to change the Rayleigh number, switch the viscosity between constant,
    temperature-dependent and temperature- and strain-rate-dependent, reseed, or
    pause. Under <em>view</em> are two overlays, both off to begin with: the
    \(\psi\) isocontours, which are the streamlines, and either mesh — the spline
    elements \(\psi\) is solved in, or the grid the temperature is carried on.
    Bottom left, the two Nusselt numbers are plotted against time on a shared
    axis: they disagree while the layer is still storing heat and converge once
    it is not, so the two curves meeting is the run arriving at a steady state.
    <em>Nu window</em> sets how much of the run that panel shows, from the last
    500 steps to all of it; narrowing it rescales the axis onto the settled band,
    which is otherwise a flat line at the top of a plot the initial transient
    sizes. Time is given twice, nondimensionally and in years — see below for
    what fixes the second, and why the figures are so large.
    The first couple of
    seconds are spent factorising the radial operators and compiling pipelines.
    <a href="/assets/mantle/">Open it full screen.</a>
  </figcaption>
</figure>
<script>
  if (!navigator.gpu) document.getElementById("webgpu-note").hidden = false;
</script>

## The model

In non-dimensional Boussinesq form, buoyancy-driven Stokes flow is quasi-static,

$$
-\nabla\cdot\big(2\mu\,\varepsilon(\vec{u})\big) + \nabla p
  = \mathrm{Ra}\;T\,\hat{\vec{e}}_r,
\qquad \nabla\cdot\vec{u} = 0,
$$

coupled to advection–diffusion of temperature,

$$
\frac{\partial T}{\partial t} + \vec{u}\cdot\nabla T = \nabla^2 T.
$$

The Rayleigh number $$\mathrm{Ra}$$ is the primary control: below its critical
value the layer simply conducts, and above it it breaks into
convection cells whose number and vigour grow with $$\mathrm{Ra}$$. The radius
ratio is $$r_i/r_o = 0.55$$, the boundaries are isothermal, and both are free
slip.

The same equations are offered on a **Cartesian box**, depth 0 to 1 with an
adjustable length, which is the cell the mantle convection literature states its
benchmarks in. It is not a second solver. Everything below is written on one
non-periodic axis $$r$$ crossed with one periodic axis $$\varphi$$, and the whole
difference between the two domains is the metric

$$
\mathrm{d}s^2 = \mathrm{d}r^2 + h(r)^2\,\mathrm{d}\varphi^2,
\qquad h = r \ \text{(annulus)}, \qquad h = 1 \ \text{(box)},
$$

so the box is the *same* discretisation evaluated with $$h' = 0$$: the curl
representation, the dissipation form, the natural free-slip condition and the
spectral solve all carry over unchanged, and the curvature terms simply drop.
One consequence is worth stating plainly, because it is what makes the box
useful here — it puts the scheme on a geometry with published benchmark numbers,
which the annulus has not got.

The box is **periodic** in the horizontal, not walled. That is structural rather
than a choice: the transverse direction is diagonalised by a radix-2 FFT, which
is precisely the statement that $$\varphi$$ wraps, and a reflecting sidewall
would need a cosine transform and a different set of radial blocks. It costs less
than it sounds — a periodic box of length $$2L$$ carrying a symmetric pair of
counter-rotating cells has $$u_x = 0$$ and $$\partial_x u_z = 0$$ on the planes
between them, which is exactly free slip, so the half-domain *is* the walled
benchmark cell.

The viscosity has three settings. Constant; temperature-dependent,
$$\mu(T) = \exp(-\gamma(T-\tfrac12))$$, with the contrast $$e^{\gamma}$$ on a
slider; and temperature- and strain-rate-dependent, which multiplies that by a
regularised power law $$\hat s^{(1-n)/n}$$ in the second invariant of the strain
rate, with $$n \approx 3$$ for dislocation creep and a viscosity floor and
ceiling.

Both exponents are *centred*, and for the same reason. The thermal one sits on
$$T = \tfrac12$$ so the geometric mean viscosity stays 1: raising the contrast
then stiffens the cold lid and weakens the hot interior symmetrically rather than
quietly rescaling the effective Rayleigh number. The power-law argument
$$\hat s$$ is the strain rate divided by its own *geometric* mean over the
domain, so that factor has geometric mean 1 too. That second choice is not
cosmetic: the exponent is negative, so viscosity is convex in strain rate, and
the domain mostly deforms more slowly than its root-mean-square — normalising by
the RMS instead stiffens nearly everywhere and, measured, cuts the peak velocity
by a factor of four. The power law is meant to redistribute viscosity into weak
shear zones and strong stagnant regions, not to shift its overall level.

The three settings are not three models so much as two *solvers*, one of them
used two ways — which is the subject of the last section below.

## Divergence free by construction

Writing the velocity as the curl of a stream function,

$$
u_r = \frac{1}{r}\frac{\partial \psi}{\partial \varphi}, \qquad
u_\varphi = -\frac{\partial \psi}{\partial r},
$$

gives $$\nabla\cdot\vec{u} = 0$$ pointwise and *identically*, since
$$\nabla\cdot\nabla\times \equiv 0$$. This is a structural property of the
representation, not of the solution: it holds for any coefficients whatsoever,
so it survives discretisation error, incomplete convergence, and — as it turns
out — the drop to single precision on the GPU. Taking the curl of the momentum
equation eliminates pressure, and for constant viscosity leaves the biharmonic
problem

$$
\mu\,\Delta^2\psi = -\,\frac{\mathrm{Ra}}{r}\,\frac{\partial T}{\partial\varphi}.
$$

Why go to this trouble for a quantity most schemes get to within discretisation
error anyway? Because [some things break on the error being nonzero](research.html),
not on its magnitude. Tracer methods lose coverage and fail to conserve
composition when the velocity carrying them is only approximately solenoidal,
and the defect accumulates over the many thousands of steps a mantle simulation
needs. An exactly divergence-free velocity removes that failure mode rather than
bounding it.

## Free slip on a curved boundary

Impermeability $$\vec{u}\cdot\vec{n} = 0$$ gives $$u_r = 0$$, hence
$$\psi_\varphi = 0$$, hence $$\psi = \text{const}$$ on each radius. Vanishing
tangential traction gives $$\sigma_{r\varphi} = 2\mu\varepsilon_{r\varphi} = 0$$,
where in polar coordinates

$$
\varepsilon_{r\varphi} = \tfrac{1}{2}\Big[\tfrac{1}{r}\partial_\varphi u_r
  + \partial_r u_\varphi - \tfrac{u_\varphi}{r}\Big].
$$

That last term is the curvature contribution, absent in Cartesian geometry.
Since $$u_r \equiv 0$$ along the whole boundary circle, the condition reduces to
$$\partial_r u_\varphi = u_\varphi / r$$, and in stream-function form

$$
\psi_{rr} - \psi_r/r = 0 \qquad \text{on } r = r_i,\ r_o.
$$

**This is not $$\omega = 0$$.** The vorticity condition reads
$$\psi_{rr} + \psi_r/r = 0$$ — the same expression with the sign of the
lower-order term flipped. The two coincide only when $$u_\varphi = 0$$, that is,
for no slip. A rigid rotation $$u_\varphi = \Omega r$$, $$u_r = 0$$ makes the
point: it has $$\varepsilon \equiv 0$$ everywhere, so it is unambiguously
stress free, yet $$\omega = 2\Omega$$. Imposing $$\omega = 0$$ would forbid
rigid rotation at a stress-free boundary. Solving in $$(r,\varphi)$$ does not
remove the curvature; it moves it into the metric.

This matters because the wrong condition is nearly invisible. It differs by one
sign in a lower-order term, and it produces smooth, plausible, entirely
believable velocity fields. Nothing blows up.

The way out is to avoid deriving the boundary term at all. Discretising the
viscous dissipation directly,

$$
a(\psi, v) = \int_\Omega 2\mu\,\varepsilon(\vec{u}[\psi]) :
             \varepsilon(\vec{u}[v])\,\mathrm{d}x,
\qquad \vec{u}[\cdot] = \nabla\times(\cdot\,\hat{\vec{z}}),
$$

has $$\sigma_{r\varphi} = 0$$ as its *natural* boundary condition, exactly, with
the curvature term carried automatically by $$\varepsilon_{r\varphi}$$. Only
$$\psi = \text{const}$$ remains essential. Pressure never appears, because
$$\nabla\cdot\vec{u}[v] \equiv 0$$ for every test function. The naive biharmonic
form $$\int \Delta\psi\,\Delta v$$, whose natural condition is
$$\Delta\psi = 0$$, is *not* the right statement here — it is missing the
curvature term.

Not deriving the boundary term is also what lets the box share this code rather
than fork it. In a box $$h' = 0$$, the offending $$u_\varphi/r$$ term is simply
absent, and $$\psi_{rr} = 0$$ and $$\Delta\psi = 0$$ do coincide — which is why
the distinction above is invisible from Cartesian geometry and easy to inherit
wrongly. The dissipation form is indifferent: it is the same integral against a
different Jacobian, so it gets both right without being told which it is in.

## Discretisation and the fast solve

The stream function is a tensor-product B-spline of degree 3, clamped in $$r$$
and uniformly periodic in $$\varphi$$. Degree 3 is the minimum that is
$$H^2$$-conforming for the fourth-order operator while also giving a continuous
strain rate, which the variable-viscosity extension will need since
$$\dot\varepsilon_{II}$$ appears differentiated there. The formulation is
Galerkin rather than collocation by necessity: cubic splines are $$C^2$$, so
$$\psi''''$$ does not exist in the space and a fourth-order operator cannot be
collocated in it.

Uniform periodic knots in $$\varphi$$ are what make the solve fast. Such a basis
is translation-invariant in the coefficient index, so the azimuthal operator
matrices are **circulant**, and circulant matrices are diagonalised exactly by
the DFT. Every coefficient in the operator — the metric's $$1/h$$ and $$1/h^2$$
— is independent of $$\varphi$$ in both geometries, so this survives the whole
problem: the DFT decouples it into
independent one-dimensional radial problems, one per azimuthal mode, each
factorised once at start-up and merely *applied* thereafter. There is no
iteration.

Two details here are easy to get wrong. The per-mode multiplier is the **symbol
of the discrete B-spline operator** — the DFT of the circulant's first row — and
not the analytic $$-k^2$$; substituting the latter is a silent consistency error
that a convergence study will catch and nothing else will. Being read off the
knot vector, that symbol also carries the period of $$\varphi$$ with it, so a box
of any length needs no correction — which is exactly the sort of place the
analytic multiplier would have gone wrong. And the $$k = 0$$
mode deserves a second look, because the continuous problem has a
two-dimensional kernel there: $$\psi = \text{const}$$, a trivial gauge, and
$$\psi = -\Omega r^2/2$$, rigid rotation, which as above is stress free. (In a
box the second is $$\psi = -Uz$$, uniform horizontal translation; the argument
below is the same one.) That is
the discretisation correctly reporting the physics, and it looks like something
that needs an explicit constraint.

It does not. Free slip makes $$\psi$$ constant on each boundary circle, and
imposing $$\psi = 0$$ on *both* — which is how the essential condition is applied
here, by dropping the two interpolatory end coefficients — removes the whole
kernel by itself, since a parabola vanishing at two distinct radii is zero. The
block is then as well conditioned as any other. That is a gauge choice, and a
consistent one: the difference $$\psi(r_o)-\psi(r_i)$$ is the net azimuthal
volume flux, so fixing it at zero says "no net circulation", and the buoyancy
load is exactly orthogonal to rigid rotation anyway — that mode's velocity is
purely azimuthal, and the load pairs only with the radial component. The physics
never asks for the part being suppressed. For constant viscosity the mode's
forcing vanishes identically in any case, since the source is
$$\propto \partial_\varphi T$$ and
$$\oint \partial_\varphi T \,\mathrm{d}\varphi = 0$$, so it is simply skipped.

**Precision.** The biharmonic operator's condition number grows like
$$h^{-4}$$, which is fatal to single-precision *elimination*: at 128 radial
degrees of freedom, $$\kappa \approx 3\times10^{8}$$ against
$$\varepsilon_{f32} \approx 6\times10^{-8}$$ is a total loss of significance. But
the radial operators are time-invariant, and JavaScript numbers are doubles, so
the factorisations are computed once at start-up in double precision on the CPU
and uploaded as single-precision *inverses*. The GPU then only ever performs a
matrix–vector product, whose error is relative to the result rather than
amplified by the conditioning. The ill-conditioning is paid for once, in double
precision, and never enters the frame loop.

## Variable viscosity: the same solver, twice

Everything above rests on the operator's coefficients being independent of
$$\varphi$$. Make the viscosity depend on temperature and that fails
immediately — $$\mu$$ now varies along $$\varphi$$, the azimuthal modes couple,
and the exact DFT decoupling is gone.

What does *not* change is the variational form. The dissipation form was chosen
in the first place because free slip is natural to it, and $$\mu(T)$$ moves
inside the integral with nothing re-derived:

$$
a(\psi, v) = \int_\Omega 2\mu(T)\,\varepsilon(\vec{u}[\psi]) :
  \varepsilon(\vec{u}[v]) \,\mathrm{d}x .
$$

That claim is worth testing rather than asserting, and it holds: the boundary
shear still converges to zero at the same second order measured for constant
viscosity. Had the curved-boundary condition needed rederiving for variable
coefficients, that is where it would have shown.

So the operator is applied **matrix-free** — evaluate the two strain-rate
components at each quadrature point, weight by $$\mu$$, and gather back to
degrees of freedom, all as a gather rather than a scatter for the same reason the
buoyancy assembly is — and the system is solved by conjugate gradients. The
preconditioner is the fast solve above, run on the azimuthally averaged
$$\bar\mu(r)$$: a radial profile keeps every coefficient $$\varphi$$-independent,
so the circulant structure and the exact DFT decoupling survive there. Constant
viscosity and $$\mu(T)$$ are therefore not two solvers but one, used twice — on
the GPU they are literally the same two pipelines, differing in which buffers are
bound.

Two things make a *fixed*, small iteration budget workable, which matters because
a convergence test would need the residual on the host and so a readback in the
hot loop. Stokes is quasi-static, so each frame starts from the previous frame's
$$\psi$$ and is already nearly right; and the preconditioner absorbs the
$$h^{-4}$$ conditioning, leaving only the viscosity variation for the iteration
to deal with. Warm-started, four iterations already reach the limit of what
single precision can represent.

Making the viscosity depend on the strain rate as well costs almost nothing on
top, because the quantity it needs has already been computed. The matrix-free
operator's first pass evaluates the two strain-rate components at every
quadrature point — and those *are* the invariant, so the power law is that pass
and a square root, three extra passes per solve against roughly seven per
conjugate-gradient iteration.

It does make the problem nonlinear, which is handled by *lagging*: viscosity is
frozen from the previous step's flow for the duration of each solve, so what the
iteration sees is still a symmetric positive-definite linear system. Repeating
the freeze-and-solve is a Picard iteration, and it is that count, not the
Krylov budget, which limits accuracy once the law is nonlinear. Measured in the
running simulation, tripling the number of conjugate-gradient iterations improves
the flow field by about a third; one extra rheology update at the original budget
improves it fourfold. This is why the geometric multigrid fallback originally
planned for high viscosity contrasts was not built: it would accelerate the part
of the calculation that had already converged.

Setting $$n = 1$$ makes the power-law factor exactly one — not nearly one — so
temperature-dependent viscosity is the $$n = 1$$ case of the general law rather
than a separate branch. The two share every buffer and every compiled shader, and
switching between them writes a single number.

There is a temptation to display the residual as a convergence readout. It would
be misleading. Since $$\psi$$ is *stored* in single precision it already sits
$$\varepsilon_{f32}$$ from anything, and the operator amplifies that by
$$\kappa \sim h^{-4}$$ — so the measured residual is flat in the iteration count
and grows with resolution, reporting the width of a mantissa rather than
anything the solver did. It is the precision argument above seen from the other
side: what makes the matrix–vector product accurate is exactly what decouples
the answer's accuracy from the residual's size. Convergence is checked in double
precision, off-line, where the quantity means what it says.

## Temperature and time stepping

Temperature is carried as point values on a uniform $$(r,\varphi)$$ grid rather
than as spline coefficients, because semi-Lagrangian advection, BFECC and
monotone limiting are all point-value operations; coefficients would force a
global interpolation solve every step and introduce damping through repeated
interpolate–evaluate cycles. Advection is a backward RK2 characteristic trace
along the exactly divergence-free velocity, with monotone-clamped bicubic
interpolation and a BFECC correction. Diffusion is implicit, and since the
conductivity is constant it decouples under the same azimuthal DFT into one
tridiagonal solve per mode.

Both are unconditionally stable, so the time step is an accuracy parameter and
not a stability limit.

### What the clock means

Everything above is nondimensional, and deliberately so: under the Boussinesq
scaling $$\mathrm{Ra}$$ is the only parameter, and a single run stands for every
combination of gravity, viscosity, layer depth and expansivity that produces it.
The equations are written with diffusion of unit strength and lengths scaled by
whichever length is 1 in code units, so one unit of time is one thermal diffusion
time across that length, $$\ell^2/\kappa$$.

**Which length that is differs between the two geometries**, and since the clock
goes as its square, the difference is a factor of nearly five rather than a
detail. The annulus has $$r_o = 1$$, so $$\ell$$ is the outer radius; the box has
depth 1, so $$\ell$$ is the thickness of the convecting layer.

Nothing in the model fixes $$\ell$$ or $$\kappa$$, so putting years on the clock
is a *display* choice and the page states the reference it uses rather than
burying it: Earth's mantle, with $$\kappa = 10^{-6}\,\mathrm{m^2/s}$$. The
geometry already matches, twice over and consistently — the annulus' radius ratio
0.55 is the core–mantle boundary against the surface, 3480/6371, and the box's
depth is that same layer, $$6371 - 3480 = 2891$$ km. The readout names the symbol
it is using, so the two clocks cannot be mistaken for one.

The resulting figures are enormous, and that is the physics rather than a slipped
decimal point. $$R_o^2/\kappa \approx 1.3\times10^{12}$$ years: thermal diffusion
across the mantle is some two orders of magnitude slower than the age of the
Earth, which is precisely why the mantle convects instead of conducting. A run
here settles after order $$0.1$$ of that, so its dimensional clock passes a
hundred billion years — no contradiction, because at
$$\mathrm{Ra} = 2\times10^4$$ this model is thousands of times more viscous than
the mantle, whose $$\mathrm{Ra}$$ is $$10^7$$–$$10^9$$. It needs many diffusion
times to do what the Earth does in a fraction of one. The dimensional reading is
honest about the scaling; it is not a claim that the simulation is the Earth.

## Streamlines for free

The contours the *view* pane draws over the temperature field are level sets of
$$\psi$$, and because $$\vec{u} = \nabla\times(\psi\hat{\vec{z}})$$ is tangent to
those level sets, they *are* the streamlines — exactly, not approximately. There
is no particle tracing, no seeding heuristic and no second buffer: the fragment
shader evaluates the spline once per pixel and uses screen-space derivatives to
give the lines a constant width at any zoom. Their spacing follows
$$\max|\psi|$$, reduced on the GPU, so the density stays legible as
$$\mathrm{Ra}$$ is swept across decades.

The mesh overlay beside it is drawn the same way, as a distance field rather than
as geometry — both discretisations are uniform in $$(r, \varphi)$$, so a line
family is the distance to the nearest multiple of an element width. It offers the
two meshes separately because they are two different objects: $$\psi$$ lives in a
spline space of $$n_r - 3$$ by $$n_\varphi$$ elements, while $$T$$ is carried on a
finer grid. Either family fades out where it approaches one line per pixel —
past that it would alias into moiré, and a blank domain is a more honest picture
of an unresolvable mesh than a false texture is.

## Verification

Each check isolates one mechanism, coarse to fine.

| Check | Result |
|---|---|
| Velocity recovery from an interpolated $$\psi$$ | 3rd order, $$\lVert\nabla\cdot\vec{u}\rVert \approx 7\times10^{-15}$$ |
| Stokes solve against a manufactured solution | 4th order, the optimal $$h^{p+1}$$ |
| Free slip: $$\varepsilon_{r\varphi} \to 0$$ under refinement | 2nd order, $$= h^{p-1}$$ as $$\psi''$$ dictates |
| Diffusion relaxing to the analytic conduction profile | $$1.2\times10^{-6}$$, with $$\mathrm{Nu}\to1$$ |
| Inner against outer Nusselt number at steady state | agree to $$1.9\times10^{-6}$$ |
| GPU against the double-precision CPU reference | $$\max\lvert\Delta T\rvert = 1.1\times10^{-6}$$ over 25 steps |
| Matrix-free $$\mu$$ operator against the assembled solve | $$A(A^{-1}b) = b$$ to $$10^{-11}$$ |
| $$\mu(T)$$ as the contrast $$\to 1$$ | recovers the direct solve in one iteration |
| $$\mu(T,\dot\varepsilon)$$ at $$n = 1$$ | recovers $$\mu(T)$$ exactly, not approximately |
| Free slip with $$\mu(T)$$ | still 2nd order |
| Free slip with $$\mu(T,\dot\varepsilon)$$ | wall traction $$\approx 3\%$$ of the interior; see below |
| Onset of convection in the box | decays at $$\mathrm{Ra}=400$$, grows at 1200, either side of $$\mathrm{Ra}_c = 27\pi^4/4$$ |
| Blankenbach 1a, $$\mathrm{Ra} = 10^4$$ | $$\mathrm{Nu} = 4.8947$$ against 4.884409, 0.2% |

The two operator checks are the ones carrying weight for variable viscosity. The
matrix-free gather and the assembled per-mode solve share no code, so composing
them — apply the operator to the direct solver's answer and see whether the load
comes back — establishes that they are the same operator, which no downstream
comparison would. And with a contrast of exactly 1 the preconditioner *is* the
operator, so conjugate gradients must land on the direct solution in a single
step; asking for one iteration turns "the answers eventually agree" into a
statement that both tiers pose the same linear problem.

The free-slip test is the one worth dwelling on. A manufactured solution carries
*inhomogeneous* boundary data and will therefore happily converge with the wrong
condition imposed; only a run with the physical homogeneous conditions
discriminates, and only through the observed *rate* — an $$\omega = 0$$
implementation plateaus at the nonzero limit $$2u_\varphi/r$$ rather than
converging. Each of these is pinned as an automated regression test asserting
the observed order, not eyeballed from a table.

With the power law that rate is not available, and the honest thing is to say so.
The viscosity floor and ceiling put a kink in $$\mu$$ along a level set of the
strain rate — and the position of that level set is itself part of the answer, so
it moves with the mesh. The coefficient is then Lipschitz rather than smooth, the
solution loses the regularity that the second-order rate depends on, and the
measured boundary shear falls roughly first order and then stops falling. It is
the *rate* that is lost, not the condition: the wall traction settles at a few
per cent of the interior shear rate, where an $$\omega = 0$$ implementation would
leave it comparable. Three other explanations were measured and excluded — an
unconverged linear solve, an unconverged Picard iteration, and the width of the
stiff layer the regularisation creates at the wall — so the test asserts what
survives rather than a number that would have to be explained away.

Convergence of the Nusselt number itself is only first order, limited by the
operator splitting and the semi-Lagrangian interpolation. In the annulus at
$$\mathrm{Ra} = 10^4$$ it extrapolates to $$\mathrm{Nu} \approx 1.52$$, a figure
that is *self*-consistent — the two boundary fluxes agree to six digits, which is
a genuine global heat balance and not a fitted quantity, and the plot in the
corner of the simulation above is that agreement being reached, live — but that
has no published annulus benchmark to be checked against.

The box is where that gap closes, and it is the reason it is here. The last two
rows of the table are the linear and the finite-amplitude halves of the same
question. Linear theory for a free-slip layer of unit depth gives a growth rate
$$\sigma = \mathrm{Ra}\,k^2/q^2 - q$$ with $$q = \pi^2 + k^2$$, hence
$$\mathrm{Ra}_c(k) = q^3/k^2$$, minimised at $$k = \pi/\sqrt2$$ where
$$\mathrm{Ra}_c = 27\pi^4/4 = 657.5$$; a box of length $$2\sqrt2$$ is one
wavelength of that mode, so a single seeded roll sits at the minimum and its fate
is decided by $$\mathrm{Ra}$$ alone. Every metric term is implicated in that
number — a stray $$1/h$$ in the operator, a wrong Jacobian in the assembly, a
transverse wavenumber that still assumed a $$2\pi$$ period — and none of them
would be visible in a run that merely looked like convection.

Above onset, **Blankenbach et al. (1989) case 1a** — unit square, free slip,
isothermal, $$\mathrm{Ra} = 10^4$$, accepted $$\mathrm{Nu} = 4.884409$$ — is
reproduced to 0.2% on a $$33\times64$$ grid, using the periodic-box trick above:
length 2, one full wavelength, and the symmetry plane between the two cells is
the wall. So the scheme is now validated against a published number and not only
against itself.

## Performance and source

At 96 × 256 spline degrees of freedom for $$\psi$$ and a 193 × 256 grid for
$$T$$, a constant-viscosity step costs well under a millisecond on a discrete
GPU, so several steps and the render fit comfortably inside a 60 fps frame. With
$$\mu(T)$$ the direct solve becomes an iteration: about 2 ms of fixed cost and a
further 0.9 ms per conjugate-gradient iteration, so roughly 13 ms per step at the
default budget of twelve — one step per frame rather than several. Adding the
strain-rate dependence on top is free within measurement noise, which is the
whole point of taking the invariant from a pass the operator was performing
anyway. The start-up pause is the double-precision factorisation and shader
compilation, paid once.

The source, together with the design document recording why each choice above
was made and what was rejected, is at
[github.com/nate-sime/nate-sime.github.io](https://github.com/nate-sime/nate-sime.github.io/tree/master/apps/mantle).
Strain-rate dependence, $$\mu(T,\dot\varepsilon)$$, adds a non-linear outer
iteration around the solver described here, and is next.
