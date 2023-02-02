#=


# Example 109: Reaction-Diffusion System

Solve the system of PDEs
```math
\partial_t u_1 = 0.5 \partial^2_x u_1 - u_1 + u_2 = 0 \\
\partial_t u_2 = 0.1 \partial^2_x u_2 + u_1 - u_2 = 0 \\
u_1(0,t) = 1 \\
\partial_x u_1(1,t) = 0 \\
\partial_x u_2(0,t) = 0 \\
u_2(0,t) = 0
````
for $x \in \Omega=(0,10)$ using the DAE solvers of the DifferentialEquations.jl package.

We take for our problem the following initial condition:
```math
u_1(x,0) = 0 \\
u_2(x,0) = 0
```
=#

module Example109_SystemReactionDiffusion_DiffEq

using SkeelBerzins


function main()

	N_x = 21
		
	L = 1
	T = 10

	x_mesh = collect(range(0,L,length=N_x))
	tspan  = (0, T)

	m = 0

	function pdefun_test(x,t,u,dudx)
		c = SVector(1,1)
		f = SVector(0.5,0.1) .* dudx
		y = u[1] - u[2]
		s = SVector(-y, y)
		
		return c,f,s
	end

	function icfun_test(x)
		u0 = SVector(0.0, 0.0)
		
		return u0
	end


	function bdfun_test(xl,ul,xr,ur,t)
		pl = SVector(ul[1]-1.0, 0)
    	ql = SVector(0, 1)
		pr = SVector(0, ur[2])
		qr = SVector(1, 0)

		return pl,ql,pr,qr
	end

	params = SkeelBerzins.Params()
	params.solver = :DiffEq

	pb = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh,tspan ; params=params)
	problem = DifferentialEquations.ODEProblem(pb)
	sol = DifferentialEquations.solve(problem,Rosenbrock23())

	return sum(sol.u[end])
end

function test()
    testval=29.035923566365785
    main() ≈ testval
end

end
