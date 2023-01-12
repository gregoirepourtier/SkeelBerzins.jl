# Achieve Performance

## Using StaticArrays.jl
One way to improve the performance of the package when solving system of PDEs, is to use the package [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) to define the PDE formulation.

This package allows you to use arrays that are stored on the stack rather than on the heap, which eliminates the need for memory allocation when defining and evaluating the problem.
Indeed normally when using arrays to define the problem, it will create allocations for every interval where we evaluate the function.

## Example
