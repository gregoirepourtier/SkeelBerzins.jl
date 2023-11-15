# Example 101: Linear Diffusion equation

Solve the following linear diffusion equation:
```math
u_t  = u_{xx}
```
for ``x \in \Omega=(0,1)`` with homogeneous Neumann boundary conditions.

We take for our problem the following initial condition:
```math
u(x,0) = exp(-100*(x-0.25)^2)
```

```
module Example101_LinearDiffusion

using SkeelBerzins, DifferentialEquations


function main()

	N_x = 21
		
	L = 1
	T = 1

	x_mesh = collect(range(0,L,length=N_x))
	tspan  = (0, T)

	m = 0
	
	fpeak(x)=exp(-100*(x-0.25)^2)

	function pdefun(x,t,u,dudx)
		c = 1
		f = dudx 
		s = 0
		
		return c,f,s
	end

	function icfun(x)
		u0 = fpeak(x)
		
		return u0
	end


	function bdfun(xl,ul,xr,ur,t)
		pl = 0
		ql = 1
		pr = 0
		qr = 1

		return pl,ql,pr,qr
	end

	params = SkeelBerzins.Params(solver=:DiffEq)

	pb         = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan ; params=params)
	problem    = DifferentialEquations.ODEProblem(pb)
	sol_diffEq = DifferentialEquations.solve(problem,Tsit5())

	sol_euler = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan)

	return (sum(sol_diffEq.u[end]), sum(sol_euler.u[end]))
end

function test()
	testval_diffEq = 3.720992375010752
	testval_euler  = 3.721004873950427
	approx_diffEq, approx_euler = main()

    approx_diffEq ≈ testval_diffEq && approx_euler ≈ testval_euler
end


end
```
