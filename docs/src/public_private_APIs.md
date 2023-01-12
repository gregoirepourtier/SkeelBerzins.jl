# Public and private APIs

We present here an index of all the methods that are present in the package. These are separated into two classes:
1. Public API
2. Private API

## Public API

### Internal implicit Euler method

```@docs
pdepe
```

### Compatibility with DifferentialEquations.jl

```@docs
SciMLBase.ODEProblem
Base.reshape
```


## Private API

These methods should only be considered for developers or people trying to understand the inner workings of the package.

### Spatial Discretization

```@docs
SkeelBerzins.assemble!
```

### Newton solvers
```@docs
SkeelBerzins.newton
SkeelBerzins.newton_stat
```

### Problem definition


### Compatibility with DifferentialEquations.jl

The `ODEFunction` constructor is implicitely defined in the `ODEProblem` and so doesn't need to be considered by the user to solve the problem with the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) package.
```@docs
SciMLBase.ODEFunction
```