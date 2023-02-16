# Tools for the package


"""
    interpolation(xl, ul, xr, ur, quadrature_point, problem)

Interpolate u and ``\\frac{du}{dx}`` between two discretization points at some specific quadrature point.

Input arguments:
- `xl`: left boundary of the current interval.
- `ul`: solution evaluated at the left boundary of the current interval.
- `xr`: right boundary of the current interval.
- `ur`: solution evaluated at the right boundary of the current interval.
- `quadrature_point`: quadrature point chosen according to the method described in [1].
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
"""
function interpolation(xl, ul, xr, ur, qd_point, pb)

    if pb.singular
        h   = xr^2 - xl^2
        tmp = qd_point^2 - xl^2

        u    = ul*(1 - (tmp/h)) + ur*(tmp/h)
        dudx = (2*qd_point*(ur-ul))/h

    elseif pb.m==0
        h = xr-xl

        u    = (ul*(xr - qd_point) + ur*(qd_point - xl))/h
        dudx = (ur-ul)/h

    elseif pb.m==1
        tmp      = log(xr/xl)
        test_fct = log(qd_point/xl)/tmp

        u    = ul*(1 - test_fct) + ur*test_fct
        dudx = (ur-ul)/(qd_point*tmp)

    elseif pb.m==2
        h        = xr-xl
        test_fct = (xr*(qd_point-xl)) / (qd_point*h)

        u    = ul*(1 - test_fct) + ur*test_fct
        dudx = (ur-ul)*((xr*xl)/(qd_point*qd_point*h))

    end

    return u, dudx
end


"""
    interpolation!(interp, d_interp, xl, ul, xr, ur, quadrature_point, problem)

Mutating version of the [`interpolation`](@ref) function.
"""
function interpolation!(interp, dinterp, xl, ul, xr, ur, qd_point, pb)

    if pb.singular
        h   = xr^2 - xl^2
        tmp = qd_point^2 - xl^2

        for i ∈ 1:pb.npde
            interp[i]  = ul[i]*(1 - (tmp/h)) + ur[i]*(tmp/h)
            dinterp[i] = (2*qd_point*(ur[i]-ul[i]))/h
        end

    elseif pb.m==0
        h = xr-xl

        for i ∈ 1:pb.npde
            interp[i]  = (ul[i]*(xr - qd_point) + ur[i]*(qd_point - xl))/h
            dinterp[i] = (ur[i]-ul[i])/h
        end

    elseif pb.m==1
        tmp      = log(xr/xl)
        test_fct = log(qd_point/xl)/tmp

        for i ∈ 1:pb.npde
            interp[i]  = ul[i]*(1 - test_fct) + ur[i]*test_fct
            dinterp[i] = (ur[i]-ul[i])/(qd_point*tmp)
        end

    elseif pb.m==2
        h        = xr-xl
        test_fct = (xr*(qd_point-xl)) / (qd_point*h)

        for i ∈ 1:pb.npde
            interp[i]  = ul[i]*(1 - test_fct) + ur[i]*test_fct
            dinterp[i] = (ur[i]-ul[i])*((xr*xl)/(qd_point*qd_point*h))
        end

    end

end


"""
$(TYPEDEF)

Mutable structure storing the problem definition.

$(TYPEDFIELDS)
"""
mutable struct ProblemDefinition{ T, Tv<:Number, Ti<:Integer, Tm<:Number, pdeFunction<:Function, 
                                                                          icFunction<:Function, 
                                                                          bdFunction<:Function }
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
    xmesh::Vector{Tv}

    """
    Time interval
    """
    tspan::Tuple{Tm,Tm}

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
    jac::Union{SparseMatrixCSC{Tv, Ti}, ExtendableSparseMatrix{Tv, Ti}, BandedMatrix{Tv}}

    """
    Evaluation of the initial condition
    """
    inival::Vector{Tv}

    """
    Interpolation points from the paper
    """
    ξ::Vector{Tv}
    ζ::Vector{Tv}
    
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

    """
    Preallocated vectors for interpolation in assemble! function when solving system of PDEs
    """
    interpolant::Vector{Tv}
    d_interpolant::Vector{Tv}
    
    ProblemDefinition{T,Tv,Ti,Tm,pdeFunction,icFunction,bdFunction}() where {T,Tv,Ti,Tm,pdeFunction,icFunction,bdFunction} = new()
end


"""
$(TYPEDEF)

Mutable structure containing all the keyword arguments for the solver [`pdepe`](@ref).

$(TYPEDFIELDS)
"""
@with_kw mutable struct Params

    """
    Choice of the time discretization either use `:euler` for internal implicit Euler method or `:DiffEq` for the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package.
    """
    solver::Symbol = :euler

    """
    Defines a time step (either pass a `Float64` or a `Vector`) when using the implicit Euler method. 
    When set to `tstep=Inf`, it solves the stationary version of the problem.
    """
    tstep::Union{Float64,Vector{Float64}} = 1e-3

    """
    Flag, returns with the solution, a list of 1d-array with the history from the newton solver.
    """
    hist::Bool = false

    """
    Choice of the type of matrix (`:sparseArrays`, `:exSparse`, `:banded`) use to store the jacobian.
    """
    sparsity::Symbol = :sparseArrays

    """
    Choice of the solver for the LSE in the newton method, see [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/).
    """
    linSolver::Union{LinearSolve.SciMLLinearSolveAlgorithm, Nothing} = nothing

    """
    Maximum number of iterations for the Newton solver.
    """
    maxit::Int = 100

    """
    Tolerance used for Newton method.
    Returns solution if ``||\\; u_{i+1} - u_{i} \\;||_2 <`` `tol`.
    """
    tol::Float64 = 1e-10

    """
    Returns the data of the PDE problem
    """
    data::Bool = false

end


"""
    implicitEuler!(y,u,problem,tau,mass_matrix,timeStep)

Assemble the system for the implicit Euler method.
"""
function implicitEuler!(y,u,pb,tau,mass,timeStep)
    assemble!(y, u, pb, timeStep)

    y = reshape(y,(pb.npde,pb.Nx))
    u = reshape(u,(pb.npde,pb.Nx))

    for i ∈ 1:pb.Nx
        for j ∈ 1:pb.npde
            y[j,i] *= -1
            y[j,i] += (1/tau)*mass[j,i]*u[j,i]
        end
    end
end


"""
    implicitEuler_stat!(y,u,problem,tau,timeStep)

Assemble the system for the implicit Euler method (variant method for stationary problems).
"""
function implicitEuler_stat!(y,u,pb,tau,timeStep)
    assemble!(y, u, pb, timeStep)

    y = reshape(y,(pb.npde,pb.Nx))
    u = reshape(u,(pb.npde,pb.Nx))

    for i ∈ 1:pb.Nx
        for j ∈ 1:pb.npde
            y[j,i] *= -1
            y[j,i] += (1/tau)*u[j,i]
        end
    end
end


"""
    newton(b, tau, timeStep, problem, mass_matrix, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false, linSol=nothing)

Newton method solving nonlinear system of equations. The Jacobi matrix used for the iteration rule is computed with
the help of the [`SparseDiffTools.jl`](https://github.com/JuliaDiff/SparseDiffTools.jl) package.

Input arguments:
- `b`: right-hand side of the system to solve.
- `tau`: constant time step used for the time discretization.
- `timeStep`: current time step of tspan.
- `problem`: Structure of type [`SkeelBerzins.ProblemDefinition`](@ref).
- `mass_matrix`: mass matrix of the problem, see [`SkeelBerzins.mass_matrix`][@ref].
- `cache`: `SparseDiffTools.ForwardColorCache`. To avoid allocating the cache in each iteration of the newton solver when computing the jacobian.
- `rhs`: preallocated vector to avoid creating allocations.

Keyword arguments:
- `tol`: tolerance or stoppping criteria (by default to `1.0e-10`).
- `maxit`: maximum number of iterations (by default to `100`).
- `hist_flag`: flag to save the history and returns it (by default to `false`).
- `linSol`: choice of the solver for the LSE, see [`LinearSolve.jl`](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/) (by default `nothing`).

Returns the solution of the nonlinear system of equations and if `hist_flag=true`, the history of the solver.
"""
function newton(b, tau, timeStep, pb, mass, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false, linSol=nothing)

    if hist_flag
	    history = Float64[]
    end

    unP1 = copy(b)

    for _ ∈ 1:maxit
        forwarddiff_color_jacobian!(pb.jac, (y,u)->implicitEuler!(y,u,pb,tau,mass,timeStep), unP1, cache)
        value!(rhs, cache)
        rhs .=  rhs .- (1 ./tau).*b

        # Solving the LSE using the LinearSolve.jl package
        prob = LinearProblem(pb.jac,rhs)
        sol1 = solve(prob,linSol)
        h = sol1.u

        unP1 .= unP1 .- h

        nm = norm(h)

        if hist_flag
		    push!(history,nm)
        end

        if nm<tol && hist_flag
            return unP1,history
        elseif nm<tol
            return unP1
        end

    end

    throw("convergence failed")
end


"""
    newton_stat(b, tau, timeStep, problem, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false, linSol=nothing)

Newton method solving nonlinear system of equations (variant of [`newton`](@ref) for stationary problems).
"""
function newton_stat(b, tau, timeStep, pb, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false, linSol=nothing)

    if hist_flag
	    history = Float64[]
    end

    unP1 = copy(b)

    for _ ∈ 1:maxit
        forwarddiff_color_jacobian!(pb.jac, (y,u)->implicitEuler_stat!(y,u,pb,tau,timeStep), unP1, cache)
        value!(rhs, cache)
        rhs .= rhs .- (1 ./tau).*b

        # Solving the LSE using the LinearSolve.jl package
        prob = LinearProblem(pb.jac,rhs)
        sol1 = solve(prob,linSol)
        h = sol1.u

        unP1 .= unP1 .- h

        nm = norm(h)

        if hist_flag
		    push!(history,nm)
        end

        if nm<tol && hist_flag
            return unP1,history
        elseif nm<tol
            return unP1
        end

    end

    throw("convergence failed")
end


"""
    pdeval(m, xmesh, u_approx, x_eval, pb)

Function that interpolates with respect to the space component the solution obtained by the solver function [`pdepe`](@ref).

Input arguments:
- `m`: symmetry of the problem (m=0,1,2 for cartesian, cylindrical or spherical).
- `xmesh`: space discretization.
- `u_approx`: approximate solution obtained by the solver `pdepe`.
- `x_eval`: point or vector of points where to interpolate the approximate solution.
- `pb`: structure defining the problem definition, see [`SkeelBerzins.ProblemDefinition`](@ref).

Returns a tuple ``(u,dudx)`` corresponding to the solution and its partial derivative with respect to the space component
evaluated in `x_eval`.
"""
function pdeval(m, xmesh, u_approx, x_eval, pb)

    if m != pb.m
        error("The symmetry argument is different from the one used to solve the PDE problem.")
    end

    u_interp_list  = []
    du_interp_list = []
    for pt_x in x_eval
        if isapprox(pt_x,xmesh[1],atol=1.0e-10*abs(xmesh[2]-xmesh[1]))
            u_interp,dudx_interp = interpolation(xmesh[1], u_approx[1], xmesh[2], u_approx[2], pt_x, pb)

            push!(u_interp_list,u_interp)
            push!(du_interp_list,dudx_interp)
        else
            idx = searchsortedfirst(xmesh,pt_x)

            if idx==1 || idx > length(xmesh)
                error("Evaluation point oustide of the spatial discretization.")
            end

            if pt_x == xmesh[idx-1]
                push!(u_interp_list,u_approx[idx-1])
                u_interp, dudx_interp = interpolation(xmesh[idx-1], ul, xmesh[idx], ur, pt_x, pb)
                push!(du_interp_list,dudx_interp)
            else
                ul = u_approx[idx-1]
                ur = u_approx[idx]
        
                u_interp, dudx_interp = interpolation(xmesh[idx-1], ul, xmesh[idx], ur, pt_x, pb)
        
                push!(u_interp_list,u_interp)
                push!(du_interp_list,dudx_interp)
            end
        end
    end

    type_tmp = typeof(x_eval)
    if type_tmp <: AbstractVector
        return u_interp_list, du_interp_list
    else
        return u_interp_list[1], du_interp_list[1]
    end
end


"""
    interpolate_sol_space(u_approx, x_eval, pb, times ; val=true, deriv=true)

Similar as [`pdeval`](@ref), i.e. interpolate the solution and its derivative with respect to the space component.

Input arguments:
- `u_approx`: approximate solution obtained by the solver `pdepe`.
- `x_eval`: point or vector of points where to interpolate the approximate solution.
- `times`: point or vector of points which denotes the time step of the solution.
- `pb`: structure defining the problem definition, see [`SkeelBerzins.ProblemDefinition`](@ref).

Keyword arguments:
- `val`: flag to return the evaluation of the solution at the `x_eval` position.
- `deriv`: flag to return the derivative with respect to the space component  of the solution at the `x_eval` position.

Returns a tuple ``(u,dudx)`` corresponding to the solution and its partial derivative with respect to the space component
evaluated in `x_eval` at time `times`.
"""
function interpolate_sol_space(u_approx, x_eval, times, pb ; val=true, deriv=true)

    u_interp_list = []
    du_interp_list = []

    u_interp_list_time  = []
    du_interp_list_time = []

    type_tmp1 = typeof(x_eval)

    cpt_tmp = 1

    for t in times

        sol = interpolate_sol_time(u_approx, t)

        for pt_x in x_eval
            if isapprox(pt_x,pb.xmesh[1],atol=1.0e-10*abs(pb.xmesh[2]-pb.xmesh[1]))
                u_interp,dudx_interp = interpolation(pb.xmesh[1], sol[1], pb.xmesh[2], sol[2], pt_x, pb)
                if val && !deriv
                    push!(u_interp_list,u_interp)
                elseif deriv && !val
                    push!(du_interp_list,dudx_interp)
                else
                    push!(u_interp_list,u_interp)
                    push!(du_interp_list,dudx_interp)
                end
            else
                idx = searchsortedfirst(pb.xmesh,pt_x)
                if idx==1 || idx > length(pb.xmesh)
                    error("Evaluation point oustide of the spatial discretization.")
                end

                ul = sol[idx-1]
                ur = sol[idx]

                if pt_x == pb.xmesh[idx-1]
                    u_interp, dudx_interp = interpolation(pb.xmesh[idx-1], ul, pb.xmesh[idx], ur, pt_x, pb)
                    if val && !deriv
                        push!(u_interp_list,u_approx[idx-1])
                    elseif deriv && !val
                        push!(du_interp_list,dudx_interp)
                    else
                        push!(u_interp_list,u_approx[idx-1])
                        push!(du_interp_list,dudx_interp)
                    end
                else
                    u_interp, dudx_interp = interpolation(pb.xmesh[idx-1], ul, pb.xmesh[idx], ur, pt_x, pb)
            
                    if val && !deriv
                        push!(u_interp_list,u_interp)
                    elseif deriv && !val
                        push!(du_interp_list,dudx_interp)
                    else
                        push!(u_interp_list,u_interp)
                        push!(du_interp_list,dudx_interp)
                    end
                end
            end
        end
        if type_tmp1 <: AbstractVector # To review for type_tmp1 and type_tmp2 <: AbstractVector?
            if val && !deriv
                push!(u_interp_list_time,u_interp_list)
            elseif deriv && !val
                push!(du_interp_list_time,du_interp_list)
            else
                push!(u_interp_list_time,u_interp_list)
                push!(du_interp_list_time,du_interp_list)
            end

        else
            if val && !deriv
                push!(u_interp_list_time,u_interp_list[cpt_tmp])
            elseif deriv && !val
                push!(du_interp_list_time,du_interp_list[cpt_tmp])
            else
                push!(u_interp_list_time,u_interp_list[cpt_tmp])
                push!(du_interp_list_time,du_interp_list[cpt_tmp])
            end
        end
        cpt_tmp += 1
    end

    type_tmp2 = typeof(times)
    if type_tmp2 <: AbstractVector
        if val && !deriv
            return u_interp_list_time
        elseif deriv && !val
            return du_interp_list_time
        else
            return u_interp_list_time, du_interp_list_time
        end
    else
        if val && !deriv
            return u_interp_list_time[1]
        elseif deriv && !val
            return du_interp_list_time[1]
        else
            return u_interp_list_time[1], du_interp_list_time[1]
        end
    end
end

(sol::AbstractDiffEqArray)(x_eval, t, pb ; val=true, deriv=true) = interpolate_sol_space(sol, x_eval, t, pb ; val=val, deriv=deriv)


"""
    interpolate_sol_time(u_approx,t)

Linear interpolatation of the solution with respect to the time component.

Input arguments:
- `u_approx`: approximate solution obtained by the solver `pdepe`.
- `t`: time ``t \\in [t_0,t_{end}]``.
"""
function interpolate_sol_time(u_approx, t)
    if isapprox(t, u_approx.t[1], atol=1.0e-10*abs(u_approx.t[2]-u_approx.t[1]))
        return u_approx[1]
    end

    idx = searchsortedfirst(u_approx.t, t)

    if idx==1 || idx > length(u_approx)
        error("The chosen time is not within the interval.")
    end

    if t==u_approx.t[idx-1]
        return u_approx[idx-1]
    else
        new_sol_interp = similar(u_approx[idx])
        dt = u_approx.t[idx]-u_approx.t[idx-1]
        x1 = (u_approx.t[idx]-t)/dt
        x0 = (t-u_approx.t[idx-1])/dt
        new_sol_interp .= x1.*u_approx[idx-1] .+ x0.*u_approx[idx]
    end
end

(sol::AbstractDiffEqArray)(t) = interpolate_sol_time(sol,t)
