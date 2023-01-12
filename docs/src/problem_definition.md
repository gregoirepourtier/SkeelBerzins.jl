# Problem Definition

Before solving the problem defined by partial differential equation(s), we define how to write it down.
Indeed the solver requires that the PDE(s) must respect a particular format.

We need to consider the following system of quasilinear partial differential equations:
```math
c\left(x,t,u,u_x\right)u_t = x^{-m}\left(x^m f(x,t,u,u_x)\right)_x + r\left(x,t,u,u_x\right)
```

with the following boundary conditions:
```math
d(x,t,u) + n(x,t)f(x,t,u,u_x) = 0
```


Warning: In the case where ``m>0``, it follows that ``a \geq 0``.


The spatial discretization from the package follows the one described in [1] resulting in second-order accurate method in space.
