#=


# Example 111: Linear Diffusion Problem in Spherical Coordinates

Solve the following problem
```math
u_t = \frac{1}{x^2}(x^2 u_x)_x \\
```
for $x \in \Omega=(0,1)$ with the imposed symmetry condition in $x=0$ (since use of spherical coordinates)
and Dirichlet condition in $x=1$ using the DAE solvers of the DifferentialEquations.jl package.

We take for our problem the following initial condition:
```math
u(x,0) = x^2
```
=#

module Example111_LinearDiffusionSpherical_DiffEq

using SkeelBerzins, DifferentialEquations
using LinearAlgebra


function main()

    Nx = 21

	L = 1
	T = 0.8

	x_mesh = collect(range(0, L, length=Nx))
	tspan  = (0.0, T)

	m = 2

    function pdefun(x,t,u,dudx)
        c = 1
        f = dudx
        s = 0
        
        return c,f,s
    end

    function icfun(x)
        u0 = x^2
        
        return u0
    end

    function bdfun(xl,ul,xr,ur,t)
        pl = 0 # ignored by solver since m=1
        ql = 0 # ignored by solver since m=1
        pr = ur-(1+6*t)
        qr = 0
    
        return pl,ql,pr,qr
    end

    params = SkeelBerzins.Params()
	params.solver = :DiffEq

    pb      = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan ; params=params)
	problem = DifferentialEquations.ODEProblem(pb)
	sol     = DifferentialEquations.solve(problem,Rosenbrock23())

    exact(x,t) = x^2 + 6*t

    return norm(sol.u[end] - exact.(x_mesh,T)) < 1.0e-14
end


function test()
	main()
end

end
