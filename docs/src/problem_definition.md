# Problem Definition

## Problem format

Before attempting to solve a problem described by partial differential equation(s) (PDE(s)), it is necessary to express the equation(s)
in a form that can be understood and processed by the solver.

Let us first define 1D discretization grids, for ``n \in \mathbb{N}``:
```math
a =: x_1 < \cdots < x_n := b \quad \text{ and } \quad t_0 < \cdots < t_{end} \; .
```

The number of points ``n`` in the spatial discretization should be chosen based on the desired level of accuracy and the complexity of the solution.

We then need to consider the following system of quasilinear partial differential equations:
```math
c(x,t,u,u_x)u_t = x^{-m}(x^m f(x,t,u,u_x))_x + s(x,t,u,u_x)
```
where ``m`` denotes the symmetry of the problem (``m=0``, ``m=1`` or ``m=2`` denoting cartesian, polar cylindrical and polar spherical coordinates respectively), while the functions ``c(x,t,u,u_x)``, ``f(x,t,u,u_x)`` and ``s(x,t,u,u_x)`` represent the capacity term, the flux term and the source term respectively.\
The capacity term ``c(x,t,u,u_x)`` must be represented as a diagonal matrix (see [`SkeelBerzins.mass_matrix`](@ref)), this means that the relationship between the partial derivatives with respect to time is limited to being multiplied by a diagonal matrix ``c(x,t,u,u_x)``.


Similarly, the boundary conditions associated to the PDE(s) must conform to a specific format, which is outlined as follows:
```math
p^i(x,t,u) + q^i(x,t)f^i(x,t,u,u_x) = 0 \quad \text{for } x=a \text{ or } x=b
```
with ``i=1,\cdots,``npde (number of PDEs).

!!! warning
    The reference [1] states that if ``m>0`` then we require ``a \geq 0`` for the method to work properly.

The package utilizes the spatial discretization method outlined in reference [1], which results in a second-order accurate method in space.

Please refer to the example provided below and the section on solvers to understand how to implement the problem in the required format.

## Example

### Linear Diffusion for ``m=0`` with homogeneous Neumann boundary conditions

We rewrite a linear diffusion problem in cartesian coordinates, with homogeneous Neumann boundary conditions, to conform to the format outlined above.

We have for ``x \in [a,b]``, ``t \in [t_0,t_{end}]``
```math
u_t(x,t) = u_{xx}(x,t) \\
\frac{du}{dx}(a,t) = \frac{du}{dx}(b,t) = 0 \\
u(x,t_0) = u_0(x)
```

In the format required by the solver, it gives:
```math
m = 0 \quad ; \\
c(x,t,u,u_x) = 1 \quad ; \\
f(x,t,u,u_x) = u_x \quad ; \\
s(x,t,u,u_x) = 0 \quad ; \\
p(x,t,u) = 0 \quad \text{ for } x=a \text{ and } x=b \quad ; \\
q(x,t) = 1 \quad \text{ for } x=a \text{ and } x=b \; .
```
