# Example 106: Linear Diffusion equation cylindrical coordinates

Solve the following problem:
```math
u_t = \frac{1}{x}(xu_x)_x
```
for ``x \in \Omega=(0,1)`` with the imposed symmetry condition in ``x=0`` (since use of cylindrical coordinates)
and Dirichlet condition in ``x=1`` using the implicit Euler method (internal method).

We initialize our problem with the exact solution (Bessel function and its first zero):
```math
u(x,0) = J_0(nx)
```
where ``n = 2.404825557695773``.


```
module Example106_LinearDiffusionCylindrical

using SkeelBerzins
using SpecialFunctions

function main()

    Nx = 21
    
    L = 1
    T = 1
    
    x_mesh = collect(range(0, L, length=Nx))
    tspan  = (0, T)
    
    m = 1

    function pdefun(x,t,u,dudx)
        c = 1
        f = dudx
        s = 0
        
        return c,f,s
    end

    function icfun(x)
        n = 2.404825557695773
        u0 = besselj(0,n*x)
        
        return u0
    end

    function bdfun(xl,ul,xr,ur,t)
        n  = 2.404825557695773
        pl = 0 # ignored by solver since m=1
        ql = 0 # ignored by solver since m=1
        pr = ur-besselj(0,n)*exp(-n^2*t)
        qr = 0

        return pl,ql,pr,qr
    end

    sol = pdepe(m,pdefun,icfun,bdfun,x_mesh,tspan)

    return sum(sol.u[end])
end


function test()
    testval = 0.04020183138531086
    main() ≈ testval
end

end
```
