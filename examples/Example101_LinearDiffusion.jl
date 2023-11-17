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

function solve_problem(m,pdefun,icfun,bdfun,xmesh,tspan,params ; linsolve=nothing)
	if params.solver == :DiffEq
		pb         = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; params=params)
		problem    = DifferentialEquations.ODEProblem(pb)
		sol = DifferentialEquations.solve(problem,Rosenbrock23(linsolve=linsolve))
	else
		sol = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; params=params)
	end

	sol
end


function main()

	N_x = 21

	x_mesh             = collect(range(0,1,length=N_x))
	x_mesh_doubleFloat = collect(Double64,range(0,1,length=N_x))

	tspan  = (0, 1)

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

	params_DiffEq             = SkeelBerzins.Params(solver=:DiffEq,sparsity=:sparseArrays)

	sol_DiffEq             = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_DiffEq)
	sol_DiffEq_doubleFloat = solve_problem(m,pdefun,icfun,bdfun,x_mesh_doubleFloat,tspan,params_DiffEq ; linsolve=SparspakFactorization())

	tstep = collect(0:1e-3:1)
	params_euler_vecTstep    = SkeelBerzins.Params(solver=:euler,tstep=tstep,hist=true)
	params_euler_fixTstep    = SkeelBerzins.Params(solver=:euler,tstep=1e-3)

	sol_euler_vecTstep, hist = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_euler_vecTstep)
	sol_euler_fixTstep       = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_euler_fixTstep)
	sol_euler_doubleFloat    = pdepe(m,pdefun,icfun,bdfun,x_mesh_doubleFloat,tspan ; solver=:euler, tstep=1e-3, linsolve=SparspakFactorization())

    if VERSION >= VersionNumber(1,9,0)
        params_DiffEq_banded = SkeelBerzins.Params(solver=:DiffEq,sparsity=:banded)
        sol_DiffEq_banded    = solve_problem(m,pdefun,icfun,bdfun,x_mesh,tspan,params_DiffEq_banded ; linsolve=LUFactorization())

        sol_euler_banded         = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan ; solver=:euler, sparsity=:banded, tstep=1e-3, linsolve=LUFactorization())

        return (sum(sol_DiffEq.u[end]),sum(sol_DiffEq_banded.u[end]),eltype(eltype(sol_DiffEq_doubleFloat.u)), 
			    sum(sol_euler_vecTstep.u[end]),sum(sol_euler_fixTstep.u[end]),sum(sol_euler_banded.u[end]),eltype(eltype(sol_euler_doubleFloat.u)))
    else
        return (sum(sol_DiffEq.u[end]),eltype(eltype(sol_DiffEq_doubleFloat.u)), 
			    sum(sol_euler_vecTstep.u[end]),sum(sol_euler_fixTstep.u[end]),eltype(eltype(sol_euler_doubleFloat.u)))
    end
end

function test()
    testval_diffEq        = 3.7210048739504296
    testval_diffEq_banded = 3.702806314278916

    testval_euler         = 3.721004873950427
    testval_euler_banded  = 3.7210048106612303

    if VERSION >= VersionNumber(1,9,0)
        approx_diffEq, approx_diffEq_banded, sol_diffEq_doubleFloat, approx_euler_vec, approx_euler, approx_euler_banded, sol_euler_doubleFloat = main()

        all_tests = approx_diffEq ≈ testval_diffEq && 
        sol_diffEq_doubleFloat == Double64 == sol_euler_doubleFloat &&
        approx_diffEq_banded ≈ testval_diffEq_banded && 
        approx_euler ≈ testval_euler ≈ approx_euler_vec &&
        approx_euler_banded ≈ testval_euler_banded
    else
        approx_diffEq, sol_diffEq_doubleFloat, approx_euler_vec, approx_euler, sol_euler_doubleFloat = main()

        all_tests = approx_diffEq ≈ testval_diffEq && 
        sol_diffEq_doubleFloat == Double64 == sol_euler_doubleFloat &&
        approx_euler ≈ testval_euler ≈ approx_euler_vec
    end

    all_tests
end


end
