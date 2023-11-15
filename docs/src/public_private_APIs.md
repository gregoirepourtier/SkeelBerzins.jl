# Public and private APIs

We present here an index of all the methods that are present in the package. These are separated into two classes:
1. Public API
2. Private API

## Public API

### Solvers

```@docs
pdepe
```

### Parameters for the package
```@docs
Params
```

### Compatibility with DifferentialEquations.jl
To define an `ODEProblem` and solve it, please refer to the package extension  
[SkeelBerzinsDiffEq.jl](https://github.com/gregoirepourtier/SkeelBerzins.jl/blob/main/ext/SkeelBerzinsDiffEq.jl) 
and the following section [Solve with DifferentialEquations.jl](@ref).

### Interpolation of the obtained solution
```@docs
pdeval
```

## Private API

These methods should only be considered for developers or people trying to understand the inner 
workings of the package.

### Spatial Discretization
```@docs
SkeelBerzins.assemble!
SkeelBerzins.interpolation
SkeelBerzins.interpolation!
```

### Newton solvers
```@docs
SkeelBerzins.newton
SkeelBerzins.newton_stat
```

### Problem definition
```@docs
SkeelBerzins.ProblemDefinition
SkeelBerzins.problem_init
SkeelBerzins.get_sparsity_pattern
SkeelBerzins.get_quad_points_weights
```

### Implicit Euler method
```@docs
SkeelBerzins.implicitEuler!
SkeelBerzins.implicitEuler_stat!
SkeelBerzins.mass_matrix
```

### Post-processing
```@docs
SkeelBerzins.interpolate_sol_time
SkeelBerzins.interpolate_sol_space
```
