# Example 106: System of Reaction-Diffusion equations

Solve the following system of PDEs:
```math
\partial_t u_1 = 0.5 \partial^2_x u_1 - u_1 + u_2 = 0 \\
\partial_t u_2 = 0.1 \partial^2_x u_2 + u_1 - u_2 = 0 \\
u_1(0,t) = 1 \\
\partial_x u_1(1,t) = 0 \\
\partial_x u_2(0,t) = 0 \\
u_2(0,t) = 0
```
for ``x \in \Omega=(0,10)``.

We take for our problem the following initial condition:
```math
u_1(x,0) = 0 \\
u_2(x,0) = 0
```

```
module Example106_SystemReactionDiffusion

using SkeelBerzins, DifferentialEquations


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

	params_diffEq = SkeelBerzins.Params()
	params_diffEq.solver = :DiffEq

	pb = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh,tspan ; params=params_diffEq)
	problem = DifferentialEquations.ODEProblem(pb)
	sol_diffEq = DifferentialEquations.solve(problem,Rosenbrock23())

	params_euler = SkeelBerzins.Params()
	params_euler.tstep = 1e-2

	sol_euler = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh,tspan ; params=params_euler)

	return (sum(sol_diffEq.u[end]),sum(sol_euler.u[end]))
end

function test()
    testval_diffEq = 29.035923566365785
	testval_euler  = 29.034702247833415

    approx_diffEq, approx_euler = main()

    approx_diffEq ≈ testval_diffEq && approx_euler ≈ testval_euler
end

end
```
