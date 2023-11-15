#=


# Example 101: Linear Diffusion 1D

Solve the linear diffusion equation
```math
u_t  = u_{xx}
```
for $x \in \Omega=(0,1)$ with homogeneous Neumann boundary conditions.

We take for our problem the following initial condition:
```math
u(x,0) = exp(-100*(x-0.25)^2)
```
=#

module Example101_LinearDiffusion

using SkeelBerzins, DifferentialEquations, DoubleFloats

function solve_problem(m,pdefun,icfun,bdfun,xmesh,tspan,params)
	if params.solver == :DiffEq
		pb         = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; params=params)
		problem    = DifferentialEquations.ODEProblem(pb)
		sol = DifferentialEquations.solve(problem,Tsit5()) # Rosenbrock23() for stiff Problems
	else
		sol = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; params=params)
	end

	sol
end


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

	params_DiffEq           = SkeelBerzins.Params(solver=:DiffEq)
	params_DiffEq_banded    = SkeelBerzins.Params(solver=:DiffEq,sparsity=:banded)
	params_DiffEq_SparseCSC = SkeelBerzins.Params(solver=:DiffEq,sparsity=:SparseArrays)

	sol_DiffEq           = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_DiffEq)
	sol_DiffEq_banded    = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_DiffEq_banded)
	sol_DiffEq_SparseCSC = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_DiffEq_SparseCSC)

	tstep = collect(0:0.1:1)
	params_euler_vecTstep    = SkeelBerzins.Params(solver=:euler,tstep=tstep,hist=true)
	params_euler_fixTstep    = SkeelBerzins.Params(solver=:euler,tstep=1e-3)

	sol_euler_vecTstep = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_euler_vecTstep)
	sol_euler_fixTstep = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_euler_fixTstep)

	return (sum(sol_DiffEq.u[end]), sum(sol_euler_fixTstep.u[end]))
end

function test()
    testval_diffEq = 3.720992375010752
	testval_euler  = 3.721004873950427
	approx_diffEq, approx_euler = main()

    approx_diffEq ≈ testval_diffEq && approx_euler ≈ testval_euler
end


end
