# Example 104: Poisson equation

Solve the following 1D Poisson equation:
```math
-u_{xx} = 1 \\
u(0)=0.1 \\
u(1)=0.1 \\
```
for ``x \in \Omega=(0,1)`` with inhomogeneous Dirichlet boundary conditions using the implicit Euler method (internal method).

```
module Example104_Poisson

using SkeelBerzins


function main()

    N_x = 21
        
    L = 1
    T = 1

    x_mesh = collect(range(0,L,length=N_x))

    m = 0

    function pdefun(x,t,u,dudx)
        c = 1
        f = dudx 
        s = 1
        
        return c,f,s
    end

    function icfun(x)
        u0 = 0.1
        
        return u0
    end


    function bdfun(xl,ul,xr,ur,t)
        pl = ul - 0.1
        ql = 0
        pr = ur - 0.1
        qr = 0

        return pl,ql,pr,qr
    end

    sol = pdepe(m,pdefun,icfun,bdfun,x_mesh)

    return sum(sol)
end

function test()
    testval=3.7624999999999997
    main() ≈ testval
end

end
```
