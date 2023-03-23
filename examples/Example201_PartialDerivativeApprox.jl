#=


# Example 201: Interpolation of Partial Derivatives

Solve the following problem
```math
u_t = (Du_x)_x - (D\frac{η}{L})u_x
```
for $x \in \Omega=(0,1)$ with homogeneous Dirichlet boundary conditions for $x=0$
and $x=1$ using the DAE solvers of the DifferentialEquations.jl package.

We take for our problem the following initial condition:
```math
u(x,0) = (K\frac{L}{D})\frac{(1 - exp(-η*(1 - \frac{x}{L})))}{η}
```

Then the function `pdeval` interpolates in space the obtained solution and its partial derivative.
=#

module Example201_PartialDerivativeApprox

using SkeelBerzins
using LinearAlgebra

function main()

    n = 50

    L = 1
    D = 0.1
    eta = 10
    K = 1
    Ip = 1

    function pdefun(x,t,u,dudx)
        c = 1
        f = D*dudx
        s = -(D*eta/L)*dudx

        c,f,s
    end

    icfun(x) = (K*L/D)*(1 - exp(-eta*(1 - x/L)))/eta

    function bcfun(xl,ul,xr,ur,t)
        pl = ul
        ql = 0
        pr = ur
        qr = 0

        pl,ql,pr,qr
    end

    m = 0
    xmesh = collect(range(0,1,length=n))
    tspan = (0,1)

    params=SkeelBerzins.Params()
    params.solver = :DiffEq

    pb = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan ; params=params)
    problem = DifferentialEquations.ODEProblem(pb)
    sol = DifferentialEquations.solve(problem,Rosenbrock23(),saveat=1/(n-1))

    function analytical(t)     
        It = 0
        for n ∈ 1:40
            m = (n*pi)^2 + 0.25*eta^2
            It = It + ((n*pi)^2 / m)* exp(-(D/L^2)*m*t)
        end
        It = 2*Ip*((1 - exp(-eta))/eta)*It
    end

    
    tmesh = collect(range(0,1,length=length(sol.t)))
    
    seriesI = analytical.(tmesh)

    u1    = zeros(length(sol.t))
    dudx1 = zeros(length(sol.t))
    u2    = zeros(length(sol.t))
    dudx2 = zeros(length(sol.t))
    for t ∈ 1:n
        u1[t], dudx1[t] = pdeval(pb.m,pb.xmesh,sol.u[t],0,pb)
        u2[t], dudx2[t] = sol(0,sol.t[t],pb)
    end

    dudx1 .= (Ip*D/K) .* dudx1
    dudx2 .= (Ip*D/K) .* dudx2

    err1 = norm(seriesI[2:n] - dudx1[2:n], Inf)
    err2 = norm(seriesI[2:n] - dudx2[2:n], Inf)

    return err1, err2
end

function test()
    testval = 0.07193510317047391
    err1,err2 = main()
    err1 ≈ err2 ≈ testval
end

end
