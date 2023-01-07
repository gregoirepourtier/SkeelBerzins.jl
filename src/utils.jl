# Tools for the package


"""
    interpolation(xl, ul, xr, ur, quadrature_point, problem)

Interpolate u and ``\frac{du}{dx}`` between two discretization points at some specific quadrature point.

Input arguments:
- 
- 
- 
- 
-
- 
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

Mutating version of the [interpolation function](@ref).
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

Mutable structure storing the problem definition
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
    Preallocated vectors for interpolation in assemble fct°
    """
    interpolant::Vector{Tv}
    d_interpolant::Vector{Tv}
    
    ProblemDefinition{T,Tv,Ti,Tm,pdeFunction,icFunction,bdFunction}() where {T,Tv,Ti,Tm,pdeFunction,icFunction,bdFunction} = new()
end


"""
    implicitEuler!(y,u,problem,tau,mass_matrix,timeStep)

Assemble the system for the implicit Euler method
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
$(SIGNATURES)

Assemble the system for the implicit Euler method for stationary problem
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
$(SIGNATURES)

Newton method solving nonlinear system of equations
"""
function newton(b, tau, timeStep, pb, mass, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false)

    if hist_flag
	    history = Float64[]
    end

    unP1 = copy(b)

    for _ ∈ 1:maxit
        forwarddiff_color_jacobian!(pb.jac, (y,u)->implicitEuler!(y,u,pb,tau,mass,timeStep), unP1, cache)
        value!(rhs, cache)
        rhs .=  rhs .- (1 ./tau).*b

        ## Solving the LSE

        # BackSlash operator
        h = pb.jac\rhs

        # LinearSolve package
        # prob = LinearProblem(pb.jac,rhs)
        # sol1 = solve(prob)
        # h = sol1.u

        # KLU Linear Solver
        # factor = klu(pb.jac)
        # h = factor\vec(rhs)

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
$(SIGNATURES)

Newton method solving nonlinear system of equations (stationary problems)
"""
function newton_stat(b, tau, timeStep, pb, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=false)

    if hist_flag
	    history = Float64[]
    end

    unP1 = copy(b)

    for _ ∈ 1:maxit
        forwarddiff_color_jacobian!(pb.jac, (y,u)->implicitEuler_stat!(y,u,pb,tau,timeStep), unP1, cache)
        value!(rhs, cache)
        rhs .= rhs .- (1 ./tau).*b

        # Solving the LSE

        # BackSlash operator    
        h = pb.jac\rhs

        # LinearSolve package
        # prob = LinearProblem(pb.jac,rhs) # add method?
        # sol1 = solve(prob)
        # h = sol1.u

        # KLU Linear Solver
        # factor = klu(pb.jac)
        # h = factor\rhs

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
