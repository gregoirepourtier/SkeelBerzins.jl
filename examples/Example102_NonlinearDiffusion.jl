#=

# Example 102: Nonlinear Diffusion 1D

Solve the following nonlinear diffusion equation
```math
u_t  = (2uu_x)_{x}
```
for $x \in \Omega=(-1,1)$ with homogeneous Neumann boundary conditions.

We take for our problem the following initial condition (exact solution named Barenblatt solution):
```math
u(x,0.001) = \max\left(0,t^{-\alpha}\left(1-\frac{\alpha(m-1)x^2}{2mt^{2\alpha}}\right)^{\frac{1}{m-1}}\right)
```
for ``m=2`` and ``\alpha = \left(m+1\right)^{-1}``.
=#

module Example102_NonlinearDiffusion

using SkeelBerzins, DifferentialEquations

function main()
    Nx = 21

    L = 1
    T = 0.01

    x_mesh = collect(range(-1, L; length=Nx))
    tspan = (0.001, T)

    m = 0

    function barenblatt(x, t, p)
        tx = t^(-1.0 / (p + 1.0))
        xx = x * tx
        xx = xx * xx
        xx = 1 - xx * (p - 1) / (2.0 * p * (p + 1))
        if xx < 0.0
            xx = 0.0
        end
        return tx * xx^(1.0 / (p - 1.0))
    end

    function pdefun(x, t, u, dudx)
        c = 1
        f = 2 * u * dudx
        s = 0

        return c, f, s
    end

    function icfun(x)
        u0 = barenblatt(x, 0.001, 2)

        return u0
    end

    function bdfun(xl, ul, xr, ur, t)
        pl = 0
        ql = 1
        pr = 0
        qr = 1

        return pl, ql, pr, qr
    end

    params_diffEq = SkeelBerzins.Params(; solver=:DiffEq)

    pb = pdepe(m, pdefun, icfun, bdfun, x_mesh, tspan; params=params_diffEq)
    problem = DifferentialEquations.ODEProblem(pb)
    sol_diffEq = DifferentialEquations.solve(problem, Tsit5())

    params_euler = SkeelBerzins.Params(; tstep=1e-4)

    sol_euler = pdepe(m, pdefun, icfun, bdfun, x_mesh, tspan; params=params_euler)

    return (sum(sol_diffEq.u[end]), sum(sol_euler.u[end]))
end

using Test

function runtests()
    testval_diffEq = 46.66666666671536
    testval_euler = 46.66666666678757
    approx_diffEq, approx_euler = main()

    @test approx_diffEq ≈ testval_diffEq && approx_euler ≈ testval_euler
end

end
