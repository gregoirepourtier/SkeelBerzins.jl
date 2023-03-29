#= 


# Example 107: 1D Poisson equation 

Solve the 1D Poisson equation
```math
-u_{xx} = 1 \\
u(0)=0.1 \\
u(1)=0.1 \\
````
for $x \in \Omega=(0,1)$ with inhomogeneous Dirichlet boundary conditions using the implicit Euler method (internal method).
=#

module Example107_Poisson

using SkeelBerzins


function main()

	N_x = 21
		
	L = 1
	T = 1

	x_mesh = collect(range(0,L,length=N_x))
	tspan  = (0, T)

	m = 0

	function pdefun_test(x,t,u,dudx)
		c = 1
		f = dudx 
		s = 1
		
		return c,f,s
	end

	function icfun_test(x)
		u0 = 0.1
		
		return u0
	end


	function bdfun_test(xl,ul,xr,ur,t)
		pl = ul - 0.1
		ql = 0
		pr = ur - 0.1
		qr = 0

		return pl,ql,pr,qr
	end

	params = SkeelBerzins.Params()
	params.tstep = Inf
	sol  = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh,tspan ; params=params)
	sol2 = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh)
	
	return sum(sol),sum(sol2)
end

function test()
    testval=3.7624999999999997
	sol1,sol2 = main()
    sol1 ≈ sol2 ≈ testval
end

end
