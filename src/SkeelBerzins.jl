module SkeelBerzins

using BandedMatrices
using LinearAlgebra
using SparseArrays
using SparseDiffTools

using RecursiveArrayTools
using DocStringExtensions
using Reexport

# @reexport using DifferentialEquations
@reexport using LinearSolve
@reexport using StaticArrays

include("pdepe.jl")
export pdepe

include("utils.jl")
export Params
export pdeval

include("assembler.jl")

# include("diffEq.jl")


end # module
