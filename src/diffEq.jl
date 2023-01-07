# Enables the DifferentialEquations.jl package to be compatible with 
# the spatial discretization from this package (solve the system of ODEs/DAEs 
# generated in the assemble! function)


"""
    ODEFunction(problem)

Generate an [ODEFunction](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction) 
from the spatial discretization derived in the [`assemble!`](@ref) function.
It is expressed as a [mass matrix ODE](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dae_example/)
and defines the problem with respect to the sparsity pattern.

Inputs arguments:
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
"""
function SciMLBase.ODEFunction(pb::ProblemDefinition)
    DifferentialEquations.ODEFunction(assemble! ; jac_prototype=pb.jac, mass_matrix=mass_matrix(pb))
end


"""
    ODEProblem(problem,callback=DifferentialEquations.CallbackSet())

Generate an [ODEProblem](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem)
from the [`ODEFunction`](@ref).

Inputs arguments:
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
- `callback`: (optional) see [callback](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).
"""
function SciMLBase.ODEProblem(pb::ProblemDefinition,callback=DifferentialEquations.CallbackSet())
    odefunction=DifferentialEquations.ODEFunction(pb)
    DifferentialEquations.ODEProblem(odefunction,pb.inival,pb.tspan,pb,callback)
end


"""
    reshape(sol_ODE, problem)

Reshape the solution `sol_ODE` obtained by the DifferentialEquations.jl package to obtain
the solution at time t as a 2D-Array of size (npde,Nx). 

Indeed since in the spatial discretization, we flattened the problem, the solution has a similar size. 
So by reshaping the solution, we get a solution organised by unknows.
"""
function Base.reshape(sol::AbstractDiffEqArray, pb::ProblemDefinition)
    RecursiveArrayTools.DiffEqArray([reshape(sol.u[i],(pb.npde,pb.Nx)) for i=1:length(sol.u)] ,sol.t)
end
