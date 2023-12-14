module SkeelBerzins

using BandedMatrices
using LinearAlgebra
using SparseArrays
using SparseDiffTools

using RecursiveArrayTools
using DocStringExtensions
using Reexport
using Symbolics

@reexport using LinearSolve
@reexport using StaticArrays

include("pdepe.jl")
export pdepe
export solve_two_scale

include("utils.jl")
export Params
export pdeval

include("problem_definition.jl")

include("bd_assembler.jl")
include("local_assembler.jl")
include("assembler_one_scale.jl")
include("assembler_two_scale.jl")

include("two_scale.jl")

if !isdefined(Base, :get_extension)
    include("../ext/SkeelBerzinsDiffEq.jl")
end

end # module
