# Solvers

Now that we have an understanding of the overall problem definition, we can explore the various 
solvers that are included in the package and investigate their specific input and output parameters 
in detail.

The solvers expect the PDE problem to be described in the following way.

## PDE formulation

In order to define the PDE(s), we have to follow the format introduced in the previous section on 
[Problem definition](problem_definition.md#problem-definition).
For the purpose of this explanation, we use the function `pdefunction(x,t,u,dudx)` to describe the 
PDE(s). The inputs of the function are self-explanatory. It will then return the capacity 
`c(x,t,u,dudx)`, the flux `f(x,t,u,dudx)` and the source `s(x,t,u,dux)` terms.

This function will be passed as an argument to the solver.

## Define the initial conditions

To define the initial condition(s), we introduce the function `icfunction(x)` (once again, arbitrary 
name). For problems that contain at least one parabolic equation, it will return the evaluation of 
the initial condition on the spatial mesh `xmesh` at the initial time ``t_0``.\
For stationary problems it will return the evaluation of the initial value on the spatial mesh 
`xmesh` used for the newton solver.

This function will be passed as an argument to the solver.

## Define the boundary conditions

To represent the boundary condition, we introduce the function `bdfunction(xl,ul,xr,ur,t)`. The 
input arguments are:
- `xl`: left boundary point of the problem.
- `ul`: estimate of the solution evaluated at the left boundary of the domain.
- `xr`: right boundary point of the problem.
- `ur`: estimate of the solution evaluated at the right boundary of the domain.
- `t`: evaluates the boundary conditions at time ``t \in [t_0,t_{end}]``.

The function will return the boundary condition terms introduced in the problem definition section, 
i.e. `p(x,t,u)` and `q(x,t)` for the left and right part of the spatial mesh.

!!! warning
    If ``m>0`` and the left boundary point of the domain is ``a=0``, then the solver ignores the 
    given boundary condition to enforce the symmetry condition ensuring second-order accuracy in 
    space near x=0.

This function will be passed as an argument to the solver.

## Obtaining Solutions

With the complete PDE formulation defined, we can now introduce the solver function `pdepe`. 
Look [`pdepe`](@ref).


### Solve with internal method: implicit Euler method

#### Parabolic equation(s)

The package contains an implementation of the implicit Euler method which can be used to solve 
parabolic equation(s). The method has a first order error with respect to time.

#### Elliptic Equation(s)

Only the internal method (implicit Euler method) can be used to solve stationary problems and is 
written as follows:

```math
M \frac{u^{k+1}-u^k}{\Delta t} = A(u^{k+1})
```
with ``M`` the mass matrix, ``\Delta t`` the time step used for the time discretization, ``A`` the 
(non)linear operator resulting from the space discretization, ``u^k`` and ``u^{k+1}`` the estimate 
solutions at time ``t_0 + k \Delta t`` and ``t_0 + (k+1) \Delta t`` respectively.

In Julia, positive infinity is defined as `Inf`. By setting ``\Delta t = `` Inf, it follows that 
``\frac{1}{\Delta t} = 0`` and thus we are left with the stationary problem which can solved by 
using the Newton solver (see [`SkeelBerzins.newton`](@ref)).

It results that the solution for the stationary problem can be obtained by running one iteration of 
the implicit Euler method.

### Solve with DifferentialEquations.jl

[SkeelBerzins.jl](https://github.com/gregoirepourtier/SkeelBerzins.jl) is also compatible with the 
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package.

It is possible to return the data from the problem in a [`SkeelBerzins.ProblemDefinition`](@ref) 
structure, then to define an ODEProblem and solve it using an ODE/DAE solver from 
DifferentialEquations.jl.

## Examples

### Solve with internal method: implicit Euler method
```
# Example of solving a linear diffusion PDE using the internal implicit Euler method.

using SkeelBerzins

# Define symmetry of the problem
m = 0

# Define PDE Formulation
function pdefunction(x,t,u,dudx)
    c = 1
    f = dudx
    s = 0

    c,f,s
end

# Define the initial condition
icfunction(x) = exp(-100*(x-0.25)^2)

# Define the boundary condtions
function bdfunction(xl,ul,xr,ur,t)
    pl = 0
    ql = 1
    pr = 0
    qr = 1

    pl,ql,pr,qr
end

# Define the spatial discretization
xmesh = collect(range(0,1,length=21))

# Define the time interval
tspan = (0,1)

# Solve
sol = pdepe(m,pdefunction,icfunction,bdfunction,xmesh,tspan)
```

### Solve with DifferentialEquations.jl
```
# Example of solving a linear diffusion PDE using the DifferentialEquations.jl package.

using SkeelBerzins, DifferentialEquations

# Define symmetry of the problem
m = 0

# Define PDE Formulation
function pdefunction(x,t,u,dudx)
    c = 1
    f = dudx
    s = 0

    c,f,s
end

# Define the initial condition
icfunction(x) = exp(-100*(x-0.25)^2)

# Define the boundary condtions
function bdfunction(xl,ul,xr,ur,t)
    pl = 0
    ql = 1
    pr = 0
    qr = 1

    pl,ql,pr,qr
end

# Define the spatial discretization
xmesh = collect(range(0,1,length=21))

# Define the time interval
tspan = (0,1)

# Define Keyword Arguments
params = SkeelBerzins.Params()
params.solver = :DiffEq

# Solve
problem_data = pdepe(m,pdefunction,icfunction,bdfunction,xmesh,tspan ; params=params)
problem_ode = DifferentialEquations.ODEProblem(problem_data)
sol_ode = DifferentialEquations.solve(problem,Rosenbrock23())
sol = reshape(sol_ode,problem_data)
```
