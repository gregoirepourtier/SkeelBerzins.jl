#=

# Example 108: Define Species on subregions
=#

module Subregions

using SkeelBerzins, ExtendableGrids, DifferentialEquations, GLMakie

function pdefun_dom(x, t, u, dudx, domain)
    eps = SVector(1, 1, 1)
    k = SVector(1, 1, 1)
    d1, d2, d3 = domain

    c = SVector(1, 1, 1)
    f = SVector(eps[1] * dudx[1], eps[2] * dudx[2], eps[3] * dudx[3])

    if d1[1] <= x <= d1[2]
        s = SVector(-k[1] * u[1] + 1.0e-4 * (3.0 - x), k[1] * u[1], 0)
    elseif d3[1] <= x <= d3[2]
        s = SVector(0, -k[3] * u[2], k[3] * u[2])
    else
        s = SVector(0, 0, 0)
    end

    c, f, s
end

icfun(x) = SVector(0, 0, 0)

function bdfun(xl, ul, xr, ur, t)
    pl = SVector(0, 0, 0)
    ql = SVector(1, 1, 1)
    pr = SVector(0, 0, ur[3])
    qr = SVector(1, 1, 0)

    pl, ql, pr, qr
end

function main(Nx; plots=false)
    xmesh = collect(range(0, 3; length=Nx))
    tspan = (0, 10)

    m = 0

    domain = ((0.0, 1.0), (1.0, 2.1), (1.9, 3.0))
    d1, d2, d3 = domain
    markers = hcat([d1[1] <= x <= d1[2] for x ∈ xmesh],
                   [d1[1] <= x <= d3[2] for x ∈ xmesh],
                   [d3[1] <= x <= d3[2] for x ∈ xmesh])

    pdefun(x, t, u, dudx) = pdefun_dom(x, t, u, dudx, domain)

    params = SkeelBerzins.Params(; solver=:DiffEq, sparsity=:symb)

    pb = pdepe(m, pdefun, icfun, bdfun, xmesh, tspan; params=params, markers_macro=markers)
    problem = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(problem, QNDF2())

    if plots
        grid = simplexgrid(xmesh)
        cellmask!(grid, [0.0], [1.0], 1)
        cellmask!(grid, [1.0], [2.1], 2)
        cellmask!(grid, [1.9], [3.0], 3)

        subgrid1 = subgrid(grid, [1])
        subgrid2 = subgrid(grid, [1, 2, 3])
        subgrid3 = subgrid(grid, [3])

        mesh1 = subgrid1.components[Coordinates][:]
        mesh2 = subgrid2.components[Coordinates][:]
        mesh3 = subgrid3.components[Coordinates][:]

        index = Observable(1)

        U1 = @lift(view(tres_sys.u[$index][1, :], subgrid1)[1:end])
        U2 = @lift(view(tres_sys.u[$index][2, :], subgrid2)[1:end])
        U3 = @lift(view(tres_sys.u[$index][3, :], subgrid3)[1:end])

        fig = Figure()
        s1 = Slider(fig[2, 1]; range=1:length(sol.t), startvalue=1)
        connect!(index, s1.value)

        title_obs = @lift("t = "*string(sol.t[$index]))

        ax = Axis(fig[1, 1]; title=title_obs, xlabel="x", ylabel="u(x,t)")

        lines!(ax, mesh1, U1)
        lines!(ax, mesh2, U2)
        lines!(ax, mesh3, U3)

        xlims!(-0.1, 3.1)
        ylims!(-5e-5, 1e-3)

        fig
    end

    xmesh, sol, pb
end

using Test

function runtests()
    xmesh, sol, pb = main(51)
    tres_sys = reshape(sol, pb)

    @test sum(tres_sys.u[end]) ≈ 0.02498057313386877
    @test length(tres_sys.t) ≈ 18
end

end
