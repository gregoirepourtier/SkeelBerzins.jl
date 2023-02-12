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

```@docs
DifferentialEquations.ODEProblem
Base.reshape
```

### Interpolation of the obtained solution
```@docs
pdeval
```

## Private API

These methods should only be considered for developers or people trying to understand the inner workings of the package.

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
```

### Compatibility with DifferentialEquations.jl

The `ODEFunction` constructor is implicitely defined in the `ODEProblem` and so doesn't need to be considered by the user to solve the problem with the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) package.
```@docs
DifferentialEquations.ODEFunction
```