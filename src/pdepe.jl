# Solver for one-dimensional ellipitic and parabolic nonlinear partial differential equations
# API methods


"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh, tspan ; params=nothing)

Solve 1D elliptic and/or parabolic partial differential equation(s) using the spatial discretization method described in [1].
The time discretization is either done by the implicit Euler method (internal method) or by using a ODE/DAE solver from the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package.
For more information on how to define the different inputs to solve a problem, look at the following sections: [Problem Definition](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/problem_definition/) and [Solvers](https://gregoirepourtier.github.io/SkeelBerzins.jl/dev/solvers/).

Input arguments:
- `m`: scalar refering to the symmetry of the problem. It can either take the value `m=0`, `m=1` or `m=2` representing 
       cartesian, cylindrical or spherical coordinates respectively.
- `pdefunction`: Function. Defines the PDE(s) formulation that incorporates capacity, flux and source terms.
- `icfunction`: Function. Defines the initial condition of the system to solve (if `tstep ```\\neq``` Inf` - initial condition from the ODE/DAE problem, 
                else if `tstep = Inf` - initial value used for the newton solver).
- `bdfunction`: Function. Defines the boundary conditions of the problem.
- `xmesh`: one-dimensional array representing the mesh on which the user wants to obtain the solution.
- `tspan`: tuple ``(t_0, t_{end})`` representing the time interval of the problem.

Keyword argument:
- `params`: defines a [`SkeelBerzins.Params`](@ref) structure containing the keyword arguments from the solvers.

Returns a [`RecursiveArrayTools.DiffEqArray`](https://docs.sciml.ai/RecursiveArrayTools/stable/array_types/#RecursiveArrayTools.DiffEqArray), a [`SkeelBerzins.ProblemDefinition`](@ref) structure
or a 1D Array, depending on the chosen solver.
Moreover, if the solution is obtained from a time dependent problem, a linear interpolation method can be use to evaluate the solution 
at any time step within the interval ``(t_0,t_{end})`` (accessible using `sol(t)`). An interpolation similar as the [`pdeval`](@ref) function is available on the solution object using the command `sol(x_eval,t,pb)`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh, tspan ; params=SkeelBerzins.Params(), mr=nothing, rmesh=nothing, pdefun_micro::T4   = nothing,
                                                                                                                            icfun_micro::T5    = nothing,
                                                                                                                            bdfun_micro::T6    = nothing,
                                                                                                                            coupling_macro::T7 = nothing,
                                                                                                                            coupling_micro::T8 = nothing,
                                                                                                                            markers            = nothing) where {T1,T2,T3,T4,T5,T6,T7,T8}

    Nx, singular, α, β, γ, npde = init_problem(m, xmesh, icfun)

    Tv = eltype(xmesh)
    Tm = eltype(tspan)
    Ti = eltype(npde)

    inival = npde==1 ? icfun.(xmesh) : vec(reduce(hcat,icfun.(xmesh))) # Reshape inival as a 1D array for compatibility with the solvers from DifferentialEquations.jl

    # display(inival)

    Nr           = nothing
    inival_micro = nothing
    if mr !== nothing
        Nr, singular_micro, α_micro, β_micro, γ_micro, npde_micro = init_problem(mr, rmesh, icfun_micro)
        inival_micro = icfun_micro.(rmesh)
        xmesh_marked = xmesh[markers]
        nx_marked = length(xmesh_marked)
    end

    inival = Nr === nothing ? inival : init_inival(inival, inival_micro, Nx, Nr, npde, markers, nx_marked, Tv)

    pb = ProblemDefinition{npde, Tv, Ti, Tm, T1, T4, T2, T5, T3, T6, T7, T8}()

    pb.npde          = npde
    pb.Nx            = Nx
    pb.xmesh         = xmesh
    pb.tspan         = tspan
    pb.singular      = singular
    pb.m             = m
    pb.inival        = inival
    
    pb.pdefunction   = pdefun
    pb.icfunction    = icfun
    pb.bdfunction    = bdfun

    pb.interpolant   = zeros(Tv,npde)
    pb.d_interpolant = zeros(Tv,npde)

    pb.Nr             = Nr
    if mr !== nothing
        pb.npde_micro     = npde_micro # only considered for npde_micro=1
        pb.rmesh          = rmesh
        pb.singular_micro = singular_micro
        pb.mr             = mr

        pb.pdefunction_micro = pdefun_micro
        pb.icfunction_micro  = icfun_micro
        pb.bdfunction_micro  = bdfun_micro

        pb.coupling_macro = coupling_macro
        pb.coupling_micro = coupling_micro

        pb.ξ_micro, pb.ζ_micro = init_quadrature_pts(mr, α_micro, β_micro, γ_micro, singular_micro)

        pb.markers      = markers
        pb.Nx_marked    = nx_marked
        pb.xmesh_marked = xmesh_marked
    end

    pb.ξ, pb.ζ = init_quadrature_pts(m, α, β, γ, singular)


    # Choosing how to initialize the jacobian with sparsity pattern (to review)
    if params.sparsity == :sparseArrays
        row    = Int64[]
	    column = Int64[]
	    vals   = Tv[]

        for i ∈ 1:npde
            for j ∈ 1:2*npde
                push!(row,i)
                push!(column,j)
                push!(vals,one(Tv))
            end
        end
        for i ∈ (Nx-1)*npde+1:Nx*npde
            for j ∈ (Nx-2)*npde+1:Nx*npde
                push!(row,i)
                push!(column,j)
                push!(vals,one(Tv))
            end
        end
        for i ∈ npde+1:npde:(Nx-1)*npde
            for k ∈ i:i+npde-1
                for j ∈ i-npde:i+2*npde-1
                    push!(row,k)
                    push!(column,j)
                    push!(vals,one(Tv))
                end
            end
        end
        pb.jac = sparse(row,column,vals)
    elseif params.sparsity == :banded
        pb.jac = BandedMatrix{Tv}(Ones(Nx*npde,Nx*npde),(2*npde-1,2*npde-1)) # Not working for general numeric datatypes
    elseif params.sparsity == :symb
        du0    = copy(inival)
	    pb.jac = Tv.(Symbolics.jacobian_sparsity((du,u)->assemble!(du,u,pb,0.0),du0,inival))
    end
    
    # tmp = 14
    # display(pb.jac[tmp:tmp+13,tmp:tmp+13])

    # Solve via implicit Euler method or an ODE/DAE solver of DifferentialEquations.jl
    if params.solver == :euler # implicit Euler method

        colors = matrix_colors(pb.jac)

        # Preallocations for Newton's method
        rhs   = zeros(Tv,npde*Nx)
        cache = ForwardColorJacCache(implicitEuler!, rhs, dx=rhs, colorvec=colors, sparsity=pb.jac)
        
        if params.tstep==Inf # Solve time independent problems (stationary case) if user set tstep to Inf -> 1 step Implicit Euler

            # To store the history of the Newton solver
            if params.hist
                storage = Tv[]
            end
            
            un = inival

            if params.hist
                unP1, history = newton_stat(un, params.tstep, 1, pb, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
            else
                unP1 = newton_stat(un, params.tstep, 1, pb, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
            end

            if params.hist && params.data
                return unP1, pb, history
            elseif params.hist
                return unP1, history
            elseif params.data
                return unP1, pb
            end

            return unP1
        
        else # Solve time dependent problems (instationary case) via Implicit Euler method

            # To store the history of the Newton solver
            if params.hist
                storage = Vector{Tv}[]
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
                results   = Vector{Tv}[]
                push!(results,un)
            else
                results   = Matrix{Tv}[]
                push!(results,reshape(inival,(npde,Nx)))
            end


            for i ∈ eachindex(timeSteps[2:end])

                t = timeSteps[i+1]

                lmul!(mass_mat_diag,un)

                if params.hist && flag_tstep
                    unP1, history = newton(un, params.tstep, t, pb, mass_mat_vec, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
                elseif flag_tstep
                    unP1          = newton(un, params.tstep, t, pb, mass_mat_vec, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
                elseif params.hist && !flag_tstep
                    unP1, history = newton(un, params.tstep[i+1], t, pb, mass_mat_vec, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
                elseif !flag_tstep
                    unP1          = newton(un, params.tstep[i+1], t, pb, mass_mat_vec, cache, rhs ; tol=params.tol, maxit=params.maxit, hist_flag=params.hist, linSol=params.linSolver)
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

    elseif params.solver == :DiffEq # Use the DifferentialEquations.jl package to solve the system of differential equations 

        # returns data from the problem to define an ODEProblem
        return pb
    end

end



"""
    pdepe(m, pdefunction, icfunction, bdfunction, xmesh ; params=nothing)

Solve 1D elliptic PDE(s) using the spatial discretization method described in [1] - [`pdepe`](@ref) variant to solve stationary problems.
Performs one step of the implicit Euler method.
For more information, look at link implicit Euler...

Input arguments:
- `m`: scalar refering to the symmetry of the problem. It can either take the value `m=0`, `m=1` or `m=2` representing 
       cartesian, cylindrical or spherical coordinates respectively.
- `pdefunction`: Function. Defines the PDE(s) formulation which includes the capacity, flux and source terms (capacity term should be set to 0).
- `icfunction`: Function. It defines the initial value used for the Newton solver.
- `bdfunction`: Function. Defines the boundary conditions of the problem.
- `xmesh`: one-dimensional array representing the mesh on which the user wants to obtain the solution.

Keyword argument:
- `params`: defines a [`SkeelBerzins.Params`](@ref) structure containing the keyword arguments from the solvers.


Returns a 1D Array with the solution at the points from the spatial discretization `xmesh`.
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh ; params=nothing, mr=nothing, rmesh=nothing, pdefun_micro::T4 = nothing,
                                                                                                       icfun_micro::T5  = nothing,
                                                                                                       bdfun_micro::T6  = nothing,
                                                                                                       coupling::T7     = nothing) where {T1,T2,T3,T4,T5,T6,T7}

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
