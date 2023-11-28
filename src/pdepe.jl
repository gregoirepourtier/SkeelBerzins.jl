# Solver for one-dimensional ellipitic and parabolic nonlinear partial differential equations
# API methods built upon MATLAB's pdepe API

"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh, tspan ; params=SkeelBerzins.Params(), kwargs...)

Solve 1D elliptic and parabolic partial differential equation(s) using the spatial discretization
method described in [1]. Note that to use this method, one of the PDEs must be parabolic.
The time discretization is either done by the implicit Euler method (internal method) or by using a
ODE/DAE solver from the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
package. For more information on how to define the different inputs to solve a problem, look at the
following sections: [Problem Definition](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/problem_definition/)
and [Solvers](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/solvers/).

Input arguments:

  - `m`: scalar referring to the symmetry of the problem. It can either take the value `m=0`, `m=1` or
    `m=2` representing cartesian, cylindrical or spherical coordinates respectively.
  - `pdefunction`: Function. Defines the formulation of the PDE(s) using capacity, flux and source
    terms.
  - `icfunction`: Function. Defines the initial condition(s) of the problem to solve (if
    `tstep!=Inf` - initial condition(s) from the ODE/DAE problem, else if
    `tstep=Inf` - initial value(s) used for the newton solver).
  - `bdfunction`: Function. Defines the boundary conditions of the problem.
  - `xmesh`: 1D array representing the spatial mesh on which the user intends to obtain the solution.
  - `tspan`: tuple ``(t_0, t_{end})`` representing the time interval of the problem.

Keyword argument:

  - `params`: defines a [`SkeelBerzins.Params`](@ref) struct containing the keyword arguments from the
    solvers.
  - `kwargs`: instead of using the [`SkeelBerzins.Params`](@ref) struct, the user can pass directly
    fields from this particular struct to the solver.

Returns a [`RecursiveArrayTools.DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray)
or a [`SkeelBerzins.ProblemDefinition`](@ref) struct, depending on the selected solver (either :euler or :DiffEq).
The obtained solution includes a linear interpolation method in time and can be used to evaluate the
solution at any time step within the interval ``(t_0,t_{end})`` (accessible using `sol(t)`).
A spatial interpolation similar as the [`pdeval`](@ref) function is available on the solution object
using the command `sol(x_eval,t,pb)`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh, tspan; params=SkeelBerzins.Params(), kwargs...) where {T1, T2, T3}

    params = (solver=haskey(kwargs, :solver) ? kwargs[:solver] : params.solver,
              tstep=haskey(kwargs, :tstep) ? kwargs[:tstep] : params.tstep,
              hist=haskey(kwargs, :hist) ? kwargs[:hist] : params.hist,
              sparsity=haskey(kwargs, :sparsity) ? kwargs[:sparsity] : params.sparsity,
              linsolve=haskey(kwargs, :linsolve) ? kwargs[:linsolve] : params.linsolve,
              maxit=haskey(kwargs, :maxit) ? kwargs[:maxit] : params.maxit,
              tol=haskey(kwargs, :tol) ? kwargs[:tol] : params.tol,
              data=haskey(kwargs, :data) ? kwargs[:data] : params.data)

    # Check if the paramater m is valid
    @assert m == 0||m == 1 || m == 2 "Parameter m invalid"
    # Check conformity of the mesh with respect to the given symmetry
    @assert m == 0||m > 0 && xmesh[1] ≥ 0 "Non-conforming mesh"
    # Check validity of time step
    @assert params.tstep≠Inf "Adjust time step or try the alternative method for solving elliptic problems"

    # Initialize Problem
    Nx, npde, inival, elTv, Ti, pb = problem_init(m, xmesh, tspan, pdefun, icfun, bdfun, params)

    # Solve time dependent problem via Implicit Euler method
    if params.solver == :euler
        colors = matrix_colors(pb.jac)::Vector{Ti}

        # Preallocations for Newton's method
        rhs = zeros(elTv, npde * Nx)
        cache = ForwardColorJacCache(implicitEuler!, rhs; dx=rhs, colorvec=colors, sparsity=pb.jac)

        # To store the history of the Newton solver
        if params.hist
            storage = Vector{elTv}[]
        end

        # Build the mass matrix
        mass_mat_diag = mass_matrix(pb)[1]
        mass_mat_vec = reshape(mass_mat_diag.diag, (npde, Nx))

        # Process Time steps
        flag_tstep = params.tstep isa Real

        timeSteps, tstep = (flag_tstep ? (collect(tspan[1]:(params.tstep):tspan[2]), params.tstep) :
                            (params.tstep, diff(params.tstep)))::Tuple{Vector{Float64}, Union{Float64, Vector{Float64}}}

        un = copy(inival)
        results = pb.npde == 1 ? Vector{elTv}[un] : Matrix{elTv}[reshape(inival, (npde, Nx))]

        for i ∈ eachindex(timeSteps[1:(end - 1)])
            time_val = timeSteps[i + 1]
            tstep_val = flag_tstep ? tstep : tstep[i]

            lmul!(mass_mat_diag, un)

            unP1 = newton(un, tstep_val, time_val, pb, mass_mat_vec, cache, rhs; tol=params.tol, maxit=params.maxit,
                          hist_flag=params.hist, linSol=params.linsolve)

            if params.hist
                pb.npde == 1 ? push!(results, unP1[1]) :
                push!(results, reshape(unP1[1], (npde, Nx)))
                push!(storage, unP1[2])
                un = unP1[1]
            else
                pb.npde == 1 ? push!(results, unP1) : push!(results, reshape(unP1, (npde, Nx)))
                un = unP1
            end
        end

        sol = RecursiveArrayTools.DiffEqArray(results, timeSteps)

        if params.hist && params.data
            return sol, pb, storage
        elseif params.hist
            return sol, storage
        elseif params.data
            return sol, pb
        end

        return sol

        # Use the DifferentialEquations.jl package to solve the system of differential equations 
    elseif params.solver == :DiffEq

        # returns data from the problem to define an ODEProblem
        return pb
    end
end

"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh ; params=SkeelBerzins.Params(tstep=Inf), kwargs...)

Solve 1D elliptic PDE(s) using the spatial discretization method described in [1] - [`pdepe`](@ref)
variant to solve stationary problems. Performs one step of the implicit Euler method.

Input arguments:

  - `m`: scalar referring to the symmetry of the problem. It can either take the value `m=0`, `m=1` or
    `m=2` representing cartesian, cylindrical or spherical coordinates respectively.
  - `pdefunction`: Function. Defines the formulation of the PDE(s) using capacity, flux and source
    terms (capacity term should be set to 0).
  - `icfunction`: Function. It defines the initial value(s) used for the Newton solver.
  - `bdfunction`: Function. Defines the boundary conditions of the problem.
  - `xmesh`: 1D array representing the spatial mesh on which the user intends to obtain the solution.

Keyword argument:

  - `params`: defines a [`SkeelBerzins.Params`](@ref) structure containing the keyword arguments from
    the solvers.
  - `kwargs`: instead of using the [`SkeelBerzins.Params`](@ref) struct, the user can pass directly
    fields from this particular struct to the solver.

Returns a 1D Array with the solution available at the points defined by the spatial discretization
`xmesh`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh; params=SkeelBerzins.Params(; tstep=Inf), kwargs...) where {T1, T2, T3}

    params = (solver=haskey(kwargs, :solver) ? kwargs[:solver] : params.solver,
              tstep=haskey(kwargs, :tstep) ? kwargs[:tstep] : Inf,
              hist=haskey(kwargs, :hist) ? kwargs[:hist] : params.hist,
              sparsity=haskey(kwargs, :sparsity) ? kwargs[:sparsity] : params.sparsity,
              linsolve=haskey(kwargs, :linsolve) ? kwargs[:linsolve] : params.linsolve,
              maxit=haskey(kwargs, :maxit) ? kwargs[:maxit] : params.maxit,
              tol=haskey(kwargs, :tol) ? kwargs[:tol] : params.tol,
              data=haskey(kwargs, :data) ? kwargs[:data] : params.data)

    @assert params.solver==:euler "Elliptic problems can only be solved using the Euler method"
    @assert params.tstep==Inf "Time step should be set to Inf to obtain the stationary solution"

    # Initialize Problem
    Nx, npde, inival, elTv, Ti, pb = problem_init(m, xmesh, (0, 1), pdefun, icfun, bdfun, params)

    colors = matrix_colors(pb.jac)::Vector{Ti}

    # Preallocations for Newton's method
    rhs = zeros(elTv, npde * Nx)
    cache = ForwardColorJacCache(implicitEuler!, rhs; dx=rhs, colorvec=colors, sparsity=pb.jac)

    unP1 = newton_stat(inival, params.tstep, 1, pb, cache, rhs; tol=params.tol, maxit=params.maxit, hist_flag=params.hist,
                       linSol=params.linsolve)

    if params.hist && params.data
        return unP1[1], pb, unP1[2]
    elseif params.hist
        return unP1
    elseif params.data
        return unP1, pb
    end

    return unP1
end
