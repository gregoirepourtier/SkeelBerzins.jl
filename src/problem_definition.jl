## ProblemDefinition Structs for different type of Problems

"""
$(TYPEDEF)

Abstract type for ProblemDefinition.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Structure storing the problem definition for one-scale problems.

$(TYPEDFIELDS)
"""
mutable struct ProblemDefinitionOneScale{T1, T2, T3, Tv <: AbstractVector, Ti <: Integer, Tm <: Number, elTv <: Number,
                                         pdeFunction <: Function,
                                         icFunction <: Function,
                                         bdFunction <: Function} <: AbstractProblemDefinition

    """
    Number of unknowns
    """
    npde::Ti

    """
    Number of discretization points
    """
    Nx::Ti

    """
    Grid of the problem
    """
    xmesh::Tv

    """
    Time interval
    """
    tspan::Tuple{Tm, Tm}

    """
    Flag to know if the problem is singular or not
    """
    singular::Bool

    """
    Symmetry of the problem
    """
    m::Ti

    """
    Jacobi matrix
    """
    jac::Union{SparseMatrixCSC{elTv, Ti}, BandedMatrix{elTv, Matrix{elTv}, Base.OneTo{Ti}}}

    """
    Evaluation of the initial condition
    """
    inival::Vector{elTv}

    """
    Interpolation points from the paper
    """
    ξ::Vector{elTv}
    ζ::Vector{elTv}

    """
    Function defining the coefficients of the PDE
    """
    pdefunction::pdeFunction

    """
    Function defining the initial condition
    """
    icfunction::icFunction

    """
    Function defining the boundary conditions
    """
    bdfunction::bdFunction

    markers_macro::Union{Vector{Bool}, Matrix{Bool}}

    function ProblemDefinitionOneScale{T1, T2, T3, Tv, Ti, Tm, elTv, pdeFunction, icFunction, bdFunction}() where {T1, T2, T3, Tv,
                                                                                                                   Ti, Tm,
                                                                                                                   elTv,
                                                                                                                   pdeFunction,
                                                                                                                   icFunction,
                                                                                                                   bdFunction}
        new()
    end
end

"""
$(TYPEDEF)

Structure storing the problem definition.

$(TYPEDFIELDS)
"""
mutable struct ProblemDefinitionTwoScale{T1, T2, T3, T4, T5, T6, Tv <: AbstractVector, Ti <: Integer, Tm <: Number,
                                         elTv <: Number, pdeFunction_macro <: Function, pdeFunction_micro <: Function,
                                         bdFunction_macro <: Function, bdFunction_micro <: Function, Coupling_macro <: Function,
                                         Coupling_micro <: Function} <: AbstractProblemDefinition
    """
    Number of unknowns
    """
    npde_macro::Ti
    npde_micro::Ti

    """
    Number of discretization points
    """
    Nx::Ti
    Nr::Ti

    """
    Number of discretization points where a "micro" PDE is defined
    """
    Nx_marked::Ti

    """
    Grid of the problem
    """
    xmesh::Tv
    rmesh::Tv

    """
    Time interval
    """
    tspan::Tuple{Tm, Tm}

    """
    Flag to know if the problem is singular or not
    """
    singular_macro::Bool
    singular_micro::Bool

    """
    Symmetry of the problem
    """
    mx::Ti
    mr::Ti

    """
    Jacobi matrix
    """
    jac::Union{SparseMatrixCSC{elTv, Ti}, BandedMatrix{elTv, Matrix{elTv}, Base.OneTo{Ti}}}

    """
    Evaluation of the initial condition
    """
    inival::Vector{elTv}

    """
    Interpolation points from the paper
    """
    ξ_macro::Vector{elTv}
    ζ_macro::Vector{elTv}
    ξ_micro::Vector{elTv}
    ζ_micro::Vector{elTv}

    """
    Function defining the coefficients of the PDE
    """
    pdefunction_macro::pdeFunction_macro
    pdefunction_micro::pdeFunction_micro

    """
    Function defining the boundary conditions
    """
    bdfunction_macro::bdFunction_macro
    bdfunction_micro::bdFunction_micro

    coupling_macro::Coupling_macro
    coupling_micro::Coupling_micro

    markers_macro::Union{Vector{Bool}, Matrix{Bool}}
    markers_micro::Vector{Bool}

    function ProblemDefinitionTwoScale{T1, T2, T3, T4, T5, T6, Tv, Ti, Tm, elTv, pdeFunction_macro, pdeFunction_micro,
                                       bdFunction_macro, bdFunction_micro, Coupling_macro, Coupling_micro}() where {T1, T2, T3,
                                                                                                                    T4, T5,
                                                                                                                    T6, Tv, Ti,
                                                                                                                    Tm, elTv,
                                                                                                                    pdeFunction_macro,
                                                                                                                    pdeFunction_micro,
                                                                                                                    bdFunction_macro,
                                                                                                                    bdFunction_micro,
                                                                                                                    Coupling_macro,
                                                                                                                    Coupling_micro
                                                                                                                    }
        new()
    end
end
