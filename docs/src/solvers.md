# Solvers

Now that we formulated the problem in the expected format, we can introduce the solvers enabling to solve
one-dimensional parabolic and elliptic nonlinear partial differential equations.


## Internal method: implicit Euler method


### Parabolic equations


### Elliptic Equations

Only the internal method will work to solve time independant problems. The implicit Euler method is writen in the following way:


It results that we get the solution from the stationary problem by running one iteration from the implicit Euler method.


## Solve with DifferentialEquations.jl

SkeelBerzins.jl is also compatible with the DifferentialEquations.jl package. It is possible to return the data from the problem in a `ProblemDefinition` structure, then to define an ODEFunction, an ODEProblem and solve it using any ODE/DAE solvers from DifferentialEquations.jl .
