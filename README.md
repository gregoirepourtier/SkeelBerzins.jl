# SkeelBerzins.jl

[![CI](https://github.com/gregoirepourtier/SkeelBerzins.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/gregoirepourtier/SkeelBerzins.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/)

Solver for one-dimensional parabolic and elliptic nonlinear partial differential equations.

This package is based on the spatial discretization for the method of lines introduced by Skeel and Berzins in [1]. For transient problems, the time discretization is performed either by the implicit Euler method (internal method) or by using an ODE/DAE solver from the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package.

# Prerequisite

The Julia version must be greater than or equal to 1.6. More information about downloading Julia can be found [here](https://julialang.org/downloads/).

# Installation

To install the package, use the following commands:

1\. From the Julian prompt:

```jldoctest
julia> using Pkg
julia> Pkg.add("SkeelBerzins")
```

2\. From the Pkg mode 

```jldoctest
julia> ]
(@v1.6) pkg> add SkeelBerzins
```

## Reference
[1] Skeel, Robert D. and Berzins, Martin, "A Method for the Spatial Discretization of Parabolic Equations in One Space Variable", _SIAM Journal on Scientific and Statistical Computing_, Vol. 11, 1990, pp.1–32.
