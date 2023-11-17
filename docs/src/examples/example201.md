# Example 201: Interpolation of Partial Derivatives


Solve the following problem
```math
u_t = (Du_x)_x - (D\frac{η}{L})u_x
```
for ``x \in \Omega=(0,1)`` with homogeneous Dirichlet boundary conditions for ``x=0`` and ``x=1`` using the DAE solvers of the DifferentialEquations.jl package.

We take for our problem the following initial condition:
```math
u(x,0) = \frac{K*L*(1 - exp(-η*(1 - \frac{x}{L})))}{D*η}
```

Then the function `pdeval` interpolates in space the obtained solution and its partial derivative.


```
module Example201_PartialDerivativeApprox

using SkeelBerzins, DifferentialEquations
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

    params=SkeelBerzins.Params(solver=:DiffEq)

    pb         = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan ; params=params)
    problem    = DifferentialEquations.ODEProblem(pb)
    sol_diffEq = DifferentialEquations.solve(problem,Rosenbrock23(),saveat=1/(n-1))

    sol_euler, pb_data_euler = pdepe(m,pdefun,icfun,bcfun,xmesh,tspan ; data=true, tstep=1/(n-1))

    @assert sol_euler.t == sol_diffEq.t "error: different time steps"

    function analytical(t)
        It = 0
        for n ∈ 1:40
            m = (n*pi)^2 + 0.25*eta^2
            It = It + ((n*pi)^2 / m)* exp(-(D/L^2)*m*t)
        end
        It = 2*Ip*((1 - exp(-eta))/eta)*It
    end


    tmesh = collect(range(0,1,length=length(sol_diffEq.t)))

    seriesI = analytical.(tmesh)

    u1,dudx1    = (zeros(length(sol_diffEq.t)),zeros(length(sol_diffEq.t)))
    u2,dudx2    = (copy(u1),copy(dudx1))
    u3,dudx3    = (copy(u1),copy(dudx1))
    u4,dudx4    = (copy(u1),copy(dudx1))
    for t ∈ 1:n
        u1[t], dudx1[t] = pdeval(pb.m,pb.xmesh,sol_diffEq.u[t],0,pb)
        u2[t], dudx2[t] = sol_diffEq(0,sol_diffEq.t[t],pb)
        u3[t], dudx3[t] = pdeval(pb_data_euler.m,pb_data_euler.xmesh,sol_euler.u[t],0,pb_data_euler)
        u4[t], dudx4[t] = sol_euler(0,sol_euler.t[t],pb_data_euler)
    end

    dudx1 .= (Ip*D/K) .* dudx1
    dudx2 .= (Ip*D/K) .* dudx2
    dudx3 .= (Ip*D/K) .* dudx3
    dudx4 .= (Ip*D/K) .* dudx4

    err1 = norm(seriesI[2:n] - dudx1[2:n], Inf)
    err2 = norm(seriesI[2:n] - dudx2[2:n], Inf)
    err3 = norm(seriesI[2:n] - dudx3[2:n], Inf)
    err4 = norm(seriesI[2:n] - dudx4[2:n], Inf)

    return err1, err2, err3, err4
end

function test()
    testval_diffEq = 0.07193510317047391
    testval_euler  = 0.6621164947313215
    err1,err2,err3,err4 = main()
    err1 ≈ err2 ≈ testval_diffEq && err3 ≈ err4 ≈ testval_euler
end

end
```
