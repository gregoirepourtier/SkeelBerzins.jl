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
    problem = DifferentialEquations.ODEProblem(pb,Rosenbrock23())
    sol = DifferentialEquations.solve(problem,saveat=1/(n-1))

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

    u    = zeros(length(sol.t))
    dudx = zeros(length(sol.t))
    for t ∈ 1:n
        u[t], dudx[t] = pdeval(pb.m,pb.xmesh,sol.u[t],0,pb)
    end

    dudx .= (Ip*D/K) .* dudx

    err = norm(seriesI[2:n] - dudx[2:n], Inf)

    return err
end


function test()
    testval = 0.0733179231620893
    main() ≈ testval
end

end
```
