module SkeelBerzins

using BandedMatrices
using LinearAlgebra
using SparseArrays
using SparseDiffTools

using RecursiveArrayTools
using DocStringExtensions
using Reexport
using Symbolics

@reexport using DifferentialEquations
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

include("diffEq.jl")

include("two_scale.jl")

end # module
