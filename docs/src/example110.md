# Example 110: System of Reaction-Diffusion equations

Solve the following system of PDEs:
```math
\partial_t u_1 = 0.5 \partial^2_x u_1 - u_1 + u_2 = 0 \\
\partial_t u_2 = 0.1 \partial^2_x u_2 + u_1 - u_2 = 0 \\
u_1(0,t) = 1 \\
\partial_x u_1(1,t) = 0 \\
\partial_x u_2(0,t) = 0 \\
u_2(0,t) = 0
```
for ``x \in \Omega=(0,10)`` using the implicit Euler method (internal method).

We take for our problem the following initial condition:
```math
u_1(x,0) = 0 \\
u_2(x,0) = 0
```

```
module Example110_SystemReactionDiffusion

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
	params.tstep = 1e-2

	sol = pdepe(m,pdefun_test,icfun_test,bdfun_test,x_mesh,tspan ; params=params)

	return sum(sol.u[end])
end

function test()
    testval=29.034702247833415
    main() ≈ testval
end

end
```
