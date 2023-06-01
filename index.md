---
layout: default
title: Home
---


# Research interests

Scientific computing, finite element methods, automatic code generation,
mathematical modelling, uncertainty quantification, parallel computing for
large scale problems, computational mathematics.


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
        <Inline nameSpaceName="Deer" mapDEFToID="true" url="assets/plots/poisson.x3d" />
    </scene>
</x3d>

## Scalable FEM

## Automated discontinuous Galerkin formulation

## Tracer methods

<!-- <div id='myDiv'>
<script src="assets/plots/divvel.js"></script>
</div> -->


## Mantle convection

<img width="50%" height="50%" src="img/out_short325_15.gif" class="center">


## Subduction zone dynamics


# CV

- 2018-February 2023, Postdoctoral Fellow, Carnegie Institution for Science.
- 2015-2018, Research Associate, The University of Cambridge.
- 2012-2015, Ph.D. Mathematics, The University of Nottingham.
- 2011-2012, P.G.Dip. Mathematics, The University of Nottingham.
- 2009-2011, Finance sector, Vancouver, B.C.
- 2005-2009, M.Sci. Theoretical Physics, The University of Nottingham.


# Software

- [FEniCS project](https://fenicsproject.org/) ([on GitHub](https://github.com/orgs/FEniCS/repositories))
- [dolfin_dg](https://github.com/nate-sime/dolfin_dg)
- [LEoPart](https://bitbucket.org/nate-sime/leopart)
- [GeoPart](https://bitbucket.org/nate-sime/geopart)
- [DOLFINx-MPC](https://github.com/jorgensd/dolfinx_mpc)
- [Subduction zone forearc modelling](https://bitbucket.org/nate-sime/subduction-zone-forea
rc-thermal-structure)
