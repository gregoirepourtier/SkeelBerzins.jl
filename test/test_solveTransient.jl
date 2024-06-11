module test_solveTransient

using SkeelBerzins, DifferentialEquations, DoubleFloats
using Test

# PDE Formulation
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

    xmesh = collect(range(0, 1; length=Nx))
    xmesh_DF = collect(Double64, range(0, 1; length=Nx))

    tspan = (0, 1)

    m = 0

    # Tests DiffEq
    params_DiffEq = SkeelBerzins.Params(; solver=:DiffEq)
    sol_DiffEq = solve_problem(m, pdefun, icfun, bdfun, xmesh, tspan, params_DiffEq)

    sol_DiffEq_DF = solve_problem(m, pdefun, icfun, bdfun, xmesh_DF, tspan, params_DiffEq; linsolve=SparspakFactorization())

    tstep = collect(0:1e-3:1)
    params_euler_vecTstep = SkeelBerzins.Params(; solver=:euler, tstep=tstep, hist=true)
    params_euler_fixTstep = SkeelBerzins.Params(; solver=:euler, tstep=1e-3)

    # Test built-in implicit Euler method
    sol_euler_vecTstep, hist1 = solve_problem(m, pdefun, icfun, bdfun, xmesh, tspan, params_euler_vecTstep)
    sol_euler_fixTstep = solve_problem(m, pdefun, icfun, bdfun, xmesh, tspan, params_euler_fixTstep)
    sol_euler_DF = pdepe(m, pdefun, icfun, bdfun, xmesh_DF, tspan; solver=:euler, tstep=1e-3, linsolve=SparspakFactorization())
    sol_euler, pb_euler, hist2 = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; solver=:euler, tstep=tstep, data=true, hist=true)

    if VERSION >= VersionNumber(1, 9, 0)
        params_DiffEq_banded = SkeelBerzins.Params(; solver=:DiffEq, sparsity=:banded)
        sol_DiffEq_banded = solve_problem(m, pdefun, icfun, bdfun, xmesh, tspan, params_DiffEq_banded; linsolve=LUFactorization())

        sol_euler_banded = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; solver=:euler, sparsity=:banded, tstep=1e-3,
                                 linsolve=LUFactorization())

        @test sum(sol_DiffEq.u[end]) ≈ 3.7210048739504296
        @test eltype(eltype(sol_DiffEq_DF.u)) == eltype(eltype(sol_euler_DF.u)) == Double64
        @test sum(sol_DiffEq_banded.u[end]) ≈ 3.72100487395043
        @test sum(sol_euler_vecTstep.u[end]) ≈ sum(sol_euler_fixTstep.u[end]) ≈ 3.721004873950427
        @test sum(sol_euler_banded.u[end]) ≈ 3.7210048739504265
        @test sum(hist1) ≈ sum(hist2)
    else
        @test sum(sol_DiffEq.u[end]) ≈ 3.7210048739504296
        @test eltype(eltype(sol_DiffEq_DF.u)) == eltype(eltype(sol_euler_DF.u)) == Double64
        @test sum(sol_euler_vecTstep.u[end]) ≈ sum(sol_euler_fixTstep.u[end]) ≈ 3.721004873950427
        @test sum(hist1) ≈ sum(hist2)
    end
end

function solve_problem(m, pdefun, icfun, bdfun, xmesh, tspan, params; linsolve=nothing)
    if params.solver == :DiffEq
        pb = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; params=params)
        problem = DifferentialEquations.ODEProblem(pb)
        sol = DifferentialEquations.solve(problem, Rosenbrock23(; linsolve=linsolve))
    else
        sol = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; params=params)
    end

    sol
end

function runtests()
    main()
end

end
