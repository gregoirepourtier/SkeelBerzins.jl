module SkeelBerzins

using BandedMatrices
using LinearAlgebra
using SparseArrays
using SparseDiffTools
using ExtendableSparse

using RecursiveArrayTools
using DocStringExtensions
using Parameters
using Reexport

@reexport using DifferentialEquations
@reexport using LinearSolve
@reexport using StaticArrays

# see https://github.com/JuliaDiff/SparseDiffTools.jl#note-about-sparse-differentiation-of-gpuarrays-bandedmatrices-and-blockbandedmatrices
using ArrayInterfaceBandedMatrices

include("pdepe.jl")
export pdepe

include("utils.jl")
export Params

include("assembler.jl")

include("diffEq.jl")


end # module
