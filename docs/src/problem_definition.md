# Problem Definition

## Problem format

Before attempting to solve numerically a problem described by partial differential equation(s), it is necessary 
to express the equation(s) in a form that can be understood and processed by the solver.

Let us first define 1D discretization grids, for ``n \in \mathbb{N}``:
```math
a =: x_1 < \cdots < x_n := b \quad \text{ and } \quad t_0 < \cdots < t_{end} \; .
```

The number of points ``n`` in the spatial discretization should be chosen based on the desired level 
of accuracy and the complexity of the solution.

We then need to consider the following system of quasilinear partial differential equations:
```math
c(x,t,u,u_x)u_t = x^{-m}(x^m f(x,t,u,u_x))_x + s(x,t,u,u_x)
```
where ``m`` designates the symmetry of the problem (``m=0``, ``m=1`` or ``m=2`` denoting cartesian, 
polar cylindrical and polar spherical coordinates respectively), while the functions 
``c(x,t,u,u_x)``, ``f(x,t,u,u_x)`` and ``s(x,t,u,u_x)`` represent the capacity, flux 
and source terms respectively.\
The capacity term ``c(x,t,u,u_x)`` must be expressed as a diagonal matrix (see 
[`SkeelBerzins.mass_matrix`](@ref)), this means that the relationship between the partial 
derivatives with respect to time is limited to being multiplied by a diagonal matrix 
``c(x,t,u,u_x)``.


Similarly, the boundary conditions associated to the PDE(s) must conform to a specific format, 
which is outlined as follows:
```math
p^i(x,t,u) + q^i(x,t)f^i(x,t,u,u_x) = 0 \quad \text{for } x=a \text{ or } x=b
```
with ``i=1,\cdots,``npde (number of PDEs).

!!! warning
    The reference [[1]](index.md) states that if ``m>0`` then we require 
    ``a \geq 0`` for the method to work properly.

The package is based on the spatial discretization method derived in reference 
[[1]](index.md), which results in a second-order accurate method in space.

Please refer to the section on [solvers](solvers.md) to understand how to implement the problem in the required 
format.
