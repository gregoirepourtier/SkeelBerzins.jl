# Enables the DifferentialEquations.jl package to solve
# the system of ODEs/DAEs generated in the assemble! function


function DifferentialEquations.ODEFunction(pb::ProblemDefinition)
    DifferentialEquations.ODEFunction(assemble! ; jac_prototype=pb.jac, mass_matrix=mass_matrix(pb))
end


function DifferentialEquations.ODEProblem(pb::ProblemDefinition,callback=DifferentialEquations.CallbackSet())
    odefunction=DifferentialEquations.ODEFunction(pb)
    DifferentialEquations.ODEProblem(odefunction,pb.inival,pb.tspan,pb,callback)
end


# To get the solution in time t as a 2D Array of size npde*Nx
function Base.reshape(sol::AbstractDiffEqArray, pb::ProblemDefinition)
    RecursiveArrayTools.DiffEqArray([reshape(sol.u[i],(pb.npde,pb.Nx)) for i=1:length(sol.u)] ,sol.t)
end
