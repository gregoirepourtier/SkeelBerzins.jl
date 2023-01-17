module SkeelBerzins

using BandedMatrices
using LinearAlgebra
using SparseArrays
using SparseDiffTools
using StaticArrays
using ExtendableSparse

using RecursiveArrayTools
using DocStringExtensions

import SciMLBase
import DifferentialEquations

# see https://github.com/JuliaDiff/SparseDiffTools.jl#note-about-sparse-differentiation-of-gpuarrays-bandedmatrices-and-blockbandedmatrices
using ArrayInterfaceBandedMatrices

# give the user the possibility to use a specific linear solver in Newton's method
using LinearSolve
using KLU

# using Symbolics

include("pdepe.jl")
export pdepe

include("utils.jl")

include("assembler.jl")

include("diffEq.jl")


end # module
