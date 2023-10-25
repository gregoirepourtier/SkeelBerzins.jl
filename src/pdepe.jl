# Solver for one-dimensional ellipitic and parabolic nonlinear partial differential equations
# API methods built upon MATLAB's pdepe API


"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh, tspan ; params=nothing)

Solve 1D elliptic and/or parabolic partial differential equation(s) using the spatial discretization 
method described in [1].
The time discretization is either done by the implicit Euler method (internal method) or by using a 
ODE/DAE solver from the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) 
package. For more information on how to define the different inputs to solve a problem, look at the 
following sections: [Problem Definition](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/
problem_definition/) and [Solvers](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/solvers/).

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

Returns a [`RecursiveArrayTools.DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/
array_types/#RecursiveArrayTools.DiffEqArray), a [`SkeelBerzins.ProblemDefinition`](@ref) struct
or a 1D Array, depending on the selected solver. Moreover, if the solution is obtained from a time 
dependent problem, a linear interpolation method can be use to evaluate the solution at any time step 
within the interval ``(t_0,t_{end})`` (accessible using `sol(t)`). An interpolation similar as the 
[`pdeval`](@ref) function is available on the solution object using the command `sol(x_eval,t,pb)`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh, tspan ; params=SkeelBerzins.Params()) where {T1,T2,T3}

    # Check if the paramater m is valid
    @assert m==0 || m==1 || m==2 "Parameter m invalid"
    # Check conformity of the mesh with respect to the given symmetry
    @assert m==0 || m>0 && xmesh[1]≥0 "Non-conforming mesh"

    # Size of the space discretization
    Nx = length(xmesh)

    # Regular case: m=0 or a>0 (Galerkin method)
    # Singular case: m≥1 and a=0 (Petrov-Galerkin method)
    singular = m≥1 && xmesh[1]==0

    α     = @view xmesh[1:end-1]
    β     = @view xmesh[2:end]
    γ = (α .+ β) ./ 2

    ξ, ζ = get_quad_points_weights(m, α, β, γ, singular)

    Tv   = typeof(xmesh)
    elTv = eltype(xmesh)
    Tm   = eltype(tspan)

    # Number of unknows in the PDE problem
    npde = length(icfun(xmesh[1]))

    # Reshape inival as a 1D array for compatibility with the solvers from DifferentialEquations.jl
    inival = npde==1 ? icfun.(xmesh) : vec(reduce(hcat,icfun.(xmesh))) 

    Ti = eltype(npde)

    if params.sparsity == :sparseArrays
        Tjac = SparseMatrixCSC{elTv, Ti}
    else
        Tjac = BandedMatrix{elTv, Matrix{elTv}, Base.OneTo{Ti}}
    end

    # Choosing how to initialize the jacobian with sparsity pattern
    jac = get_sparsity_pattern(Tjac, Nx, npde, elTv)

    pb = ProblemDefinition{npde, Tv, Ti, Tm, Tjac, T1, T2, T3}(
        npde,
        Nx,
        xmesh,
        tspan,
        singular,
        m,
        jac,
        inival,
        ξ,
        ζ,
        pdefun,
        icfun,
        bdfun,
        similar(xmesh,npde),
        similar(xmesh,npde),
    )

    # Solve via implicit Euler method or an ODE/DAE solver of DifferentialEquations.jl
    if params.solver == :euler # implicit Euler method

        colors = matrix_colors(pb.jac)::Vector{Ti}

        # Preallocations for Newton's method
        rhs   = zeros(elTv,npde*Nx)
        cache = ForwardColorJacCache(implicitEuler!, rhs ; dx=rhs, colorvec=colors, sparsity=pb.jac)
        
        # Solve time independent problem (stationary case); set tstep to Inf -> 1 step Implicit Euler
        if params.tstep==Inf 

            # To store the history of the Newton solver
            if params.hist
                storage = elTv[]
            end
            
            un = inival

            if params.hist
                unP1, history = newton_stat(un, 
                                            params.tstep, 
                                            1, 
                                            pb, 
                                            cache, 
                                            rhs ; 
                                            tol=params.tol, 
                                            maxit=params.maxit, 
                                            hist_flag=params.hist, 
                                            linSol=params.linSolver)
            else
                unP1 = newton_stat(un, 
                                   params.tstep, 
                                   1, 
                                   pb, 
                                   cache, 
                                   rhs ; 
                                   tol=params.tol, 
                                   maxit=params.maxit, 
                                   hist_flag=params.hist, 
                                   linSol=params.linSolver)
            end

            if params.hist && params.data
                return unP1, pb, history
            elseif params.hist
                return unP1, history
            elseif params.data
                return unP1, pb
            end

            return unP1
        
        # Solve time dependent problems (instationary case) via Implicit Euler method
        else 

            # To store the history of the Newton solver
            if params.hist
                storage = Vector{elTv}[]
            end

            # Build the mass matrix
            mass_mat_diag, flag_DAE = mass_matrix(pb)
            mass_mat_vec  = reshape(mass_mat_diag.diag,(npde,Nx))

            flag_tstep = params.tstep isa Real

            if flag_tstep
                timeSteps = collect(tspan[1]:params.tstep:tspan[2])
            else
                timeSteps = params.tstep
                params.tstep = diff(params.tstep)
            end
            
            un = copy(inival)
            if pb.npde == 1
                results   = Vector{elTv}[]
                push!(results,un)
            else
                results   = Matrix{elTv}[]
                push!(results,reshape(inival,(npde,Nx)))
            end


            for i ∈ eachindex(timeSteps[2:end])

                t = timeSteps[i+1]

                lmul!(mass_mat_diag,un)

                if params.hist && flag_tstep
                    unP1, history = newton(un, 
                                           params.tstep, 
                                           t, 
                                           pb, 
                                           mass_mat_vec, 
                                           cache, 
                                           rhs ; 
                                           tol=params.tol, 
                                           maxit=params.maxit, 
                                           hist_flag=params.hist, 
                                           linSol=params.linSolver)
                elseif flag_tstep
                    unP1          = newton(un, 
                                           params.tstep, 
                                           t, 
                                           pb, 
                                           mass_mat_vec, 
                                           cache, 
                                           rhs ; 
                                           tol=params.tol, 
                                           maxit=params.maxit, 
                                           hist_flag=params.hist, 
                                           linSol=params.linSolver)
                elseif params.hist && !flag_tstep
                    unP1, history = newton(un, 
                                           params.tstep[i+1], 
                                           t, 
                                           pb, 
                                           mass_mat_vec, 
                                           cache, 
                                           rhs ; 
                                           tol=params.tol, 
                                           maxit=params.maxit, 
                                           hist_flag=params.hist, 
                                           linSol=params.linSolver)
                elseif !flag_tstep
                    unP1          = newton(un, 
                                           params.tstep[i+1], 
                                           t, 
                                           pb, 
                                           mass_mat_vec, 
                                           cache, 
                                           rhs ; 
                                           tol=params.tol, 
                                           maxit=params.maxit, 
                                           hist_flag=params.hist, 
                                           linSol=params.linSolver)
                end

                if pb.npde == 1
                    push!(results,unP1)
                else
                    push!(results,reshape(unP1,(npde,Nx)))
                end
                
                un = unP1
                
                if params.hist
                    push!(storage,history)
                end
            end

            sol = RecursiveArrayTools.DiffEqArray(results,timeSteps)

            if params.hist && params.data
                return sol, pb, storage
            elseif params.hist
                return sol, storage
            elseif params.data
                return sol, pb
            end

            return sol
        end
    
    # Use the DifferentialEquations.jl package to solve the system of differential equations 
    elseif params.solver == :DiffEq 

        # returns data from the problem to define an ODEProblem
        return pb
    end

end



"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh ; params=nothing)

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

Returns a 1D Array with the solution available at the points defined by the spatial discretization 
`xmesh`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh ; params=nothing) where {T1,T2,T3}

    tspan = (0.0,1.0) # define an arbitrary time interval

    if params === nothing
        params = Params()
        params.tstep = Inf
    else
        @assert params.tstep == Inf "Time step should be set to Inf to obtain the stationary solution"
    end

    sol = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; params=params)

    return sol
end
