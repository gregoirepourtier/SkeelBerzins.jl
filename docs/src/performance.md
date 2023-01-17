# Achieve Performance

## Using StaticArrays.jl

One way to improve the performance of the package when solving system of PDEs, is to use the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to define the PDE formulation.

This package allows you to use arrays that are stored on the stack rather than on the heap, which eliminates the need for memory allocation when defining and evaluating the problem.
In contrast, when using arrays to define the problem, it creates heap-allocations for every interval where the function is evaluated.

## Example

```
using SkeelBerzins, StaticArrays

# Define symmetry of the problem
m = 0

# Define PDE Formulation
function pdefunction(x,t,u,dudx)
    c = SVector(1,1)
    f = SVector(0.5,0.1) .* dudx
    y = u[1] - u[2]
    s = SVector(-y, y)

    c,f,s
end

# Define the initial condition
function icfunction(x)
    u0 = SVector(0.0, 0.0)
    
    u0
end

# Define the boundary condtions
function bdfunction(xl,ul,xr,ur,t)
    pl = SVector(ul[1]-1.0, 0)
    ql = SVector(0, 1)
    pr = SVector(0, ur[2])
    qr = SVector(1, 0)

    pl,ql,pr,qr
end

# Define the spatial discretization
xmesh = collect(range(0,1,length=21))

# Define the time interval
tspan = (0,10)

# Solve
sol = pdepe(m,pdefunction,icfunction,bdfunction,xmesh,tspan ; solver=:euler,tstep=1e-2)
```
