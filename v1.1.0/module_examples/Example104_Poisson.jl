#=

# Example 104: 1D Poisson equation 

Solve the following 1D Poisson equation:
```math
\begin{aligned}
-u_{xx} &= 1\\
u(0) &= 0.1\\
u(1) &= 0.1.
\end{aligned}
```
for $x \in \Omega=(0,1)$ with inhomogeneous Dirichlet boundary conditions using the implicit Euler method (internal method).
=#

module Example104_Poisson

using SkeelBerzins

function pdefun(x, t, u, dudx)
    c = 1
    f = dudx
    s = 1

    return c, f, s
end

function icfun(x)
    u0 = 0.1

    return u0
end

function bdfun(xl, ul, xr, ur, t)
    pl = ul - 0.1
    ql = 0
    pr = ur - 0.1
    qr = 0

    return pl, ql, pr, qr
end

function main()
    Nx = 21

    xmesh = collect(range(0, 1; length=Nx))

    m = 0

    sol = pdepe(m, pdefun, icfun, bdfun, xmesh)

    return sum(sol)
end

using Test

function runtests()
    testval = 3.7624999999999997
    @test main() ≈ testval
end

end
