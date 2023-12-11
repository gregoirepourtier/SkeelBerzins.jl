module Example401_TwoScaleProblem

using SkeelBerzins, DifferentialEquations

## Equations defined on the Macro-Scale
function pdefun_macro(x, t, u, dudx)
    c = SVector(1, 1, 1)
    f = SVector(0.5, 0.1, 1) .* dudx
    s = SVector(-u[1] + u[2], u[1] - u[2], u[2] - u[3]) # ignored only on marked meshpoints --> given by coupling_macro function

    return c, f, s
end

icfun_macro(x) = [0, 0, 0]

function bdfun_macro(xl, ul, xr, ur, t)
    pl = [ul[1] - 1.0, 0, ul[3]]
    ql = [0, 1, 0]
    pr = [0, ur[2], 0]
    qr = [1, 0, 1]

    return pl, ql, pr, qr
end

## Equation defined on the Micro-Scale
function pdefun_micro(x, t, u, dudx)
    c = 1
    f = dudx
    s = 0

    return c, f, s
end

icfun_micro(x, ν) = 0.0

function bdfun_micro(xl, ul, xr, ur, t)
    pl = 0 # ignored symmetry condition
    ql = 1 # ignored symmetry condition
    pr = 0 # ignored --> given by coupling_micro function
    qr = 1

    return pl, ql, pr, qr
end

## Solve the two scale problem
function main_two_scale()

    Nx = 50
    xmesh = LinRange(0, 1, Nx)

    Nr = 15
    rmesh = LinRange(0, 1, Nr)

    tspan = (0, 2)

    m_x = 0
    m_r = 2

    coupling_macro(x, t, u, v) = SVector(-u[1] + v[end], u[1] - u[2], u[2] - u[3])
    coupling_micro(x, t, u, v) = -u[1] + v[end]

    params = SkeelBerzins.Params(; solver=:DiffEq, sparsity=:symb)
    markers = [xmesh[i] < 0.3 || xmesh[i] > 0.6 for i ∈ eachindex(xmesh)]

    pb = pdepe(m_x, pdefun_macro, icfun_macro, bdfun_macro, xmesh, tspan;
               params=params,
               mr=m_r,
               rmesh=rmesh,
               pdefun_micro=pdefun_micro,
               icfun_micro=icfun_micro,
               bdfun_micro=bdfun_micro,
               coupling_macro=coupling_macro,
               coupling_micro=coupling_micro,
               markers_micro=markers)
    problem = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(problem, Rosenbrock23())

    return sol, pb, markers, pb.npde
end

using Test

function runtests()
    sol, pb, markers_test, npde_macro = main_two_scale()
    tsol = reshape(sol, pb)

    @test sum(tsol[1].u[end]) ≈ 41.653053489954075
    @test sum(tsol[2].u[end]) ≈ 23.995904466062015
    @test sum(tsol[3].u[end]) ≈ 4.4143546273724255
    @test sum(tsol[end].u[end]) ≈ 9.6914384940458
end

end
