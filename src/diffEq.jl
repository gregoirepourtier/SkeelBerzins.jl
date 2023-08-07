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
        total = pb.Nx*pb.npde + pb.Nx_marked*pb.Nr
        index_macro = zeros(Bool,(pb.npde, total))
        index_micro = zeros(Bool,(pb.Nx_marked, total))

        for j=1:pb.npde
            k    = j
            for i ∈ 1:pb.Nx
                index_macro[j,k] = true
                
                if pb.markers[i]
                    k += pb.npde + pb.Nr
                else
                    k += pb.npde
                end
            end
            push!(solutions, RecursiveArrayTools.DiffEqArray([sol.u[i][index_macro[j,:]] for i=1:length(sol.u)],sol.t))
        end

        cpt_marker = 1
        k          = pb.npde+1
        for idx ∈ 1:pb.Nx
            if pb.markers[idx]
                index_micro[cpt_marker, k:k+pb.Nr-1] .= true

                push!(solutions, RecursiveArrayTools.DiffEqArray([sol.u[i][index_micro[cpt_marker,:]] for i=1:length(sol.u)],sol.t))

                cpt_marker += 1
                k += pb.npde + pb.Nr
            else
                k += pb.npde
            end
        end

        solutions
    end
end
