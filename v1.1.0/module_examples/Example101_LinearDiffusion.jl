#=

# Example 101: Linear Diffusion 1D

Solve the following linear diffusion equation
```math
u_t  = u_{xx}
```
for $x \in \Omega=(0,1)$ with homogeneous Neumann boundary conditions.

We take for our problem the following initial condition:
```math
u(x,0) = \exp(-100*(x-0.25)^2)
```
=#

module Example101_LinearDiffusion

using SkeelBerzins, DifferentialEquations

function pdefun(x, t, u, dudx)
    c = 1
    f = dudx
    s = 0

    return c, f, s
end

icfun(x) = exp(-100 * (x - 0.25)^2)

function bdfun(xl, ul, xr, ur, t)
    pl = 0
    ql = 1
    pr = 0
    qr = 1

    return pl, ql, pr, qr
end

function main()
    Nx = 21

    xmesh = LinRange(0, 1, Nx)
    tspan = (0, 1)

    m = 0

    ## Solve Problem using Rosenbrock-W solver from DifferentialEquations.jl
    pb = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; solver=:DiffEq)
    problem = DifferentialEquations.ODEProblem(pb)
    sol_diffEq = DifferentialEquations.solve(problem, Rosenbrock23())

    ## Solve Problem using the built-in implicit Euler method
    sol_euler = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan)

    (sum(sol_diffEq.u[end]), sum(sol_euler.u[end]))
end

using Test

function runtests()
    testval_diffEq = 3.7210048739504296

    testval_euler = 3.721004873950427

    sol_diffEq, sol_euler = main()

    @test sol_diffEq ≈ testval_diffEq
    @test sol_euler ≈ testval_euler
end

end
