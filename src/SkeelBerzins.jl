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

include("utils.jl")
export Params
export pdeval

include("bd_assembler.jl")
include("local_assembler.jl")
include("assembler.jl")

include("two_scale.jl")

if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" begin
            include("../ext/SkeelBerzinsDiffEq.jl")
        end
    end
end

end # module
