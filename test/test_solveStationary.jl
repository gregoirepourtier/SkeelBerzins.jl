module test_solveStationary

using SkeelBerzins
using Test

function pdefun(x, t, u, dudx)
    c = 0
    f = dudx
    s = 1

    c, f, s
end

icfun(x) = 0.1

function bdfun(xl, ul, xr, ur, t)
    pl = ul - 0.1
    ql = 0
    pr = ur - 0.1
    qr = 0

    pl, ql, pr, qr
end

function main()
    Nx = 21

    xmesh = collect(range(0, 1; length=Nx))

    m = 0

    sol1 = pdepe(m, pdefun, icfun, bdfun, xmesh)
    sol2, hist2 = pdepe(m, pdefun, icfun, bdfun, xmesh; hist=true)
    sol3, pb3, hist3 = pdepe(m, pdefun, icfun, bdfun, xmesh; hist=true, data=true)
    sol4, pb4 = pdepe(m, pdefun, icfun, bdfun, xmesh; data=true)

    @test sum(sol1) ≈ sum(sol2) ≈ sum(sol3) ≈ sum(sol4) ≈ 3.7624999999999997
    @test sum(hist2) ≈ sum(hist3)
    @test typeof(pb3) == typeof(pb4)
end

function runtests()
    main()
end

end
