
## FEniCS project

The [FEniCS project](https://fenicsproject.org/) is a general toolbox
comprising a number of libraries for the computation of the finite element
discretisation of partial differential euqations. For example, consider the
Poisson problem. The the finite element formulation reads: find $$u \in V^h$$
such that

$$
a(u, v) = l(v) \quad \forall v \in V^h.
$$

where in the context of the Poisson problem we have

$$
\begin{aligned}
a(u, v) &= \int_\Omega \nabla u \cdot \nabla v \; \mathrm{d} x, \\
l(v) &= \int_\Omega f \; v \; \mathrm{d} x.
\end{aligned}
$$

By using FEniCS this formulation is represented in python exploiting
symbolic algebra as implemented in the unified form language
([UFL](https://github.com/FEniCS/ufl))


```python
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
l = ufl.inner(f, v) * ufl.dx
```

This symbolic representation is then translated to `C` code implementation of
the local element kernel by the FEniCSx form compiler
([FFCx](https://github.com/FEniCS/ffcx))

<details closed>
	<summary><i>Click to expand FFCx generated C code of the bilinear formulation</i></summary>
{% highlight C %}

{% include_relative assets/codeblock/poisson_tensor.c %}

{% endhighlight %}
</details>
<p/>


The generated code may then be used by
[DOLFINx](https://github.com/FEniCS/dolfinx) to assemble the global system
matrix and vector for solution by external linear algebra packages.
The following figure is generated from the computation of the FE
approximation of the solution of the Poisson problem as exhibited
in the DOLFINx [Poisson demo](https://github.com/FEniCS/dolfinx/blob/main/python/demo/demo_poisson.py).


<x3d width='$(window).height();' height='$(window).width();'>
    <scene>
        <Viewpoint 
        position="1.0 0.5 2.0"
        centerOfRotation="1.0 0.5 0.3"
        description="Poisson"></Viewpoint>
        <Inline nameSpaceName="Poisson" mapDEFToID="true" url="assets/plots/poisson.x3d" />
    </scene>
</x3d>

## Scalable FEM

Crucial to modern numerical simulation is computational *scaling*. This means
that the time and memory requried to compute the solution must have a linear
relationship with problem size. In the context of FE simulations, one must
ensure that every component of the computation scales, e.g.: file I/O, mesh
partitioning and distribution, degree of freedom map construction, linear
system assembly, and its solution by computational linear algebra. Every
component of the [DOLFINx](https://github.com/FEniCS/dolfinx) library in the
[FEniCS project](https://fenicsproject.org/) project is designed to scale
optimally or provide an interface to scalable third party libraries.

Consider, for example, that FE discretisations of second order elliptic
problems which yield $$n$$ unknowns (degrees of freedom) in the underlying
linear system. Solution by a direct method (e.g. Gaussian factorisation) has
complexities $$\mathcal{O}(n^{\frac{3}{2}})$$ and $$\mathcal{O}(n^2)$$
for 2D and 3D domains, respectively.

Although the direct solver may perform exceptionally well for a problem with $$n
\backsim 10^5$$, what if our scientific application demands the fidelity only
offered by a problem with $$n \geq 10^9$$? For these large problems, even if
the growth constant is small, the growth rates of $$\mathcal{O}(n^\frac{3}{2})$$
and $$\mathcal{O}(n^2)$$ cannot be considered merely a disadvantage. They ensure
that computing a solution of a large system with a direct method is impossible
in reasonable time.

The key is to implement appropriate iterative schemes combined with a suitably
constructed preconditioner. This technique offers (close to) optimal
$$\mathcal{O}(n)$$ complexity; however, at the cost of carefully designing the
preconditioner. With this scalable implementation we are able to compute large
scale 3D simulations of physical models in science and engineering. The following
videos show the temperatures and von Mises stresses computed in a turbocharger
over one cycle of operation where $$n > 3\times10^9$$.

<p align="center">
<iframe src="https://streamable.com/s/pstwv/ncphbn" frameborder="0" allowfullscreen></iframe>
<iframe src="https://streamable.com/s/b3adh/tiixzj" frameborder="0" allowfullscreen></iframe>
</p>


## Automated discontinuous Galerkin (DG) formulation

Standard conforming FE methods employ a basis which is $$C^0$$ continuous.
This implies that the FE solution space is a finite dimensional subspace of
the Sobolev space $$V^h \subset H^1$$ found in the weak formulation. DG
methods employ a basis which is discontinuous at the boundaries of the cells
in the FE mesh. In this setting the finite dimensional DG space
$$V^h_{\text{DG}} \subset L_2$$ and therefore $$V^h_{\text{DG}} \not\subset
V^h$$. This offers the key advantage that we may seek the FE solution in
a richer space of functions than standard conforming methods.

However the question remains, how do we tie together these discontinuous basis
functions? One such scheme is the interior penalty method which derives from
Nitsche's method. This yields a formulation which conserves numerical fluxes
through the facets between cells. These additional terms defined on the facets are
typically difficult to write out, even with tools such as
[UFL](https://github.com/FEniCS/ufl). This problem of complexity of implementation
is particularly pertinent in the cases of multiphysics problems.

The library [dolfin_dg](https://github.com/nate-sime/dolfin_dg) seeks to
address this verbosity for which DG methods are notorious. By providing the
flux functions found in the underlying PDE, the DG FE formulation is
automatically generated. This gives us tools to solve highly nonlinear problem
with naturally stabilised interior penalty methods.

![Image: Compressible Euler cylinder](img/fine.png) *Density field shock front
capture of a supersonic cylinder modelled by the compressible Euler equations*

![Image: Compressible Navier Stokes NACA0012](img/naca0012.png) *Density field of
a NACA0012 airfoil modelled by the compressible Navier Stokes system (Cf. [dg_naca0012_2d.py](https://github.com/nate-sime/dolfin_dg/blob/master/demo/dolfin/compressible_navier_stokes_naca0012/dg_naca0012_2d.py)).*

<!--  This also offers stabilisation
(affectionately called
'[upwinding](https://en.wikipedia.org/wiki/Upwind_scheme)') as an inherent
component of the formulation. -->

## Tracer methods & Mantle convection

Tracer methods are widespread in computational geodynamics, particularly for
modeling the advection of chemical data. However, they present certain
numerical challenges: the necessity for mass conservation of composition
fields and the need for the velocity field to be pointwise divergence free to
avoid gaps in tracer coverage. The issues are especially noticed in
simulations over many time steps.

The issue of conservation may be addressed by $$l_2$$ projection of the tracer
data to a FE function whilst constrained by the linear advection PDE. I.e.,
given $$N_p$$ particles, composition field $$\phi_h$$ and velocity field
$$\vec{u}$$, at time $$t$$ we solve

$$
\begin{gather}
\min_{\phi_h \in V^h_\text{DG}} \mathcal{J}(\phi_h) =
\sum^{N_p}_{p} \frac{1}{2} (\phi_h(\vec{x}_p(t), t) - \phi_p(t))^2, \\
\text{subject to:} \quad \frac{\partial \phi_h}{\partial t} + \vec{u} \cdot \nabla \phi_h = 0.
\end{gather}
$$

The linear advection equation presents us with the subsequent challenge of
mass conservation. Note that

$$
\vec{u} \cdot \nabla \phi_h = \nabla \cdot (\vec{u} \phi_h) - \phi_h \nabla \cdot \vec{u}.
$$

A velocity field which exactly solves the Stokes system is solenoidal
($$\nabla \cdot \vec{u}$$) and therefore will exactly conserve mass. However,
our computed approximations of the velocity, $$\vec{u}_h$$, may not. We
therefore must choose a numerical scheme which exactly satisfies the
conservation of mass pointwise such that $$\nabla \cdot \vec{u}_h(\vec{x}) = 0$$
for all $$\vec{x}$$ in the domain of interest.


![Image: Taylor Hood vs Solenoidal velocity piles](img/piles.png) *Snapshot of
a simulation of the Earth's mantle where tracers track compositional species.
The top and base of the domain correspond to the Earth's surface and the
core-mantle boundary, respectively. Left: spurious results encountered in
mantle convection where $$\nabla \cdot \vec{u}_h(\vec{x}) \neq 0$$. Right:
$$\nabla \cdot \vec{u}_h \equiv 0$$ showing aggregation of material into
"piles" at the core-mantle boundary.*

<img width="60%" height="65%" src="img/out_short325_15.gif" class="center">
*Mantle convection simulation. The top slice shows advection of compositional
tracers. The lower slice depicts the evolution of the mantle temperature where
lighter and darker colours represent hot and cold material, respectively. Note
the hot upwellings entraining chemical material to the Earth's surface, and
material sinking from the surface in cold downwellings to be trapped into
piles at the core-mantle boundary.*