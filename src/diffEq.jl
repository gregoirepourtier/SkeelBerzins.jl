# Enables the DifferentialEquations.jl package to be compatible with 
# the spatial discretization from this package (solve the system of ODEs/DAEs 
# generated in the assemble! function)


"""
    ODEFunction(problem)

Generate an [ODEFunction](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEFunction) 
from the spatial discretization derived in the [`assemble!`](@ref) function.
It is expressed as a [mass matrix ODE](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dae_example/) if the mass matrix is different
from the identity matrix or as a simple [system of ODEs](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/advanced_ode_example/#stiff) otherwise and defines the problem with respect to the sparsity pattern.

Input argument:
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
"""
function DifferentialEquations.ODEFunction(pb::ProblemDefinition)
    massMatrix, flag_DAE = mass_matrix(pb)
    if flag_DAE
        DifferentialEquations.ODEFunction(assemble! ; jac_prototype=pb.jac, mass_matrix=massMatrix)
    else
        DifferentialEquations.ODEFunction(assemble! ; jac_prototype=pb.jac)
    end
end


"""
    ODEProblem(problem,callback=DifferentialEquations.CallbackSet())

Generate an [ODEProblem](https://docs.sciml.ai/DiffEqDocs/stable/types/ode_types/#SciMLBase.ODEProblem)
from the [`ODEFunction`](@ref) which then can be solved by using the [DifferentialEquations.solve()](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#CommonSolve.solve-Tuple{SciMLBase.AbstractDEProblem,%20Vararg{Any}}) method.

Input arguments:
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
- `callback`: (optional) see [callback](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/).
"""
function DifferentialEquations.ODEProblem(pb::ProblemDefinition)
    odefunction=DifferentialEquations.ODEFunction(pb)
    DifferentialEquations.ODEProblem(odefunction,pb.inival,pb.tspan,pb)
end


"""
    reshape(sol_ODE, problem)

Reshape the solution `sol_ODE` obtained by the DifferentialEquations.jl package to obtain
the solution at time t as a 2D-Array of size (npde,Nx). 

Indeed since in the spatial discretization, we flattened the problem, the solution has a similar size. 
So by reshaping the solution, we get a solution organised by unknows.
"""
function Base.reshape(sol::AbstractDiffEqArray, pb::ProblemDefinition)
    if pb.Nr === nothing
        RecursiveArrayTools.DiffEqArray([reshape(sol.u[i],(pb.npde,pb.Nx)) for i=1:length(sol.u)] ,sol.t)
    else
        solutions = RecursiveArrayTools.DiffEqArray[]
        for j=1:pb.npde
            push!(solutions, RecursiveArrayTools.DiffEqArray([sol.u[i][j:pb.npde+pb.Nr:end] for i=1:length(sol.u)],sol.t))
        end
        for j=pb.npde+1:pb.Nr+pb.npde:pb.Nx*(pb.npde+pb.Nr)
            push!(solutions, RecursiveArrayTools.DiffEqArray([sol.u[i][j:j+pb.Nr-1] for i=1:length(sol.u)],sol.t))
        end
        solutions
    end
end
