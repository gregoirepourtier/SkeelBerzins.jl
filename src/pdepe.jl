# Solver for one-dimensional ellipitic and parabolic nonlinear partial differential equations
# API method


"""
$(SIGNATURES)

Solve 1D elliptic and parabolic PDE(s) using the spatial discretization method described in [1].
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh, tspan ; solver=:euler, tstep=1e-3, hist=false, sparsity=:sparseArray) where {T1,T2,T3}

    # Check if the paramater m is valid
    @assert m==0 || m==1 || m==2 "Parameter m invalid"
    # Check conformity of the mesh with respect to the given symmetry
    @assert m==0 || m>0 && xmesh[1]≥0 "Non-conforming mesh"

    # Size of the space discretization
    Nx = length(xmesh)

    # Regular case: m=0 or a>0 (Galerkin method) ; Singular case: m≥1 and a=0 (Petrov-Galerkin method)
    singular = m≥1 && xmesh[1]==0

    α     = @view xmesh[1:end-1]
    β     = @view xmesh[2:end]
    gamma = (α .+ β) ./ 2

    # Number of unknows in the PDE problem
    npde = length(icfun.(xmesh[1]))

    inival_tmp = icfun.(xmesh)
    inival     = zeros(npde,Nx)
    for i ∈ 1:Nx
        inival[:,i] .= inival_tmp[i]
    end

    # Reshape inival as a one-dimensional array to make it compatible with the solvers from DifferentialEquations.jl
    inival = vec(inival)

    Tv = eltype(xmesh)
    Ti = eltype(npde)
    Tm = eltype(tspan)

    pb = ProblemDefinition{npde, Tv, Ti, Tm, T1, T2, T3}()

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

    pb.interpolant   = zeros(npde)
    pb.d_interpolant = zeros(npde)

    # Cartesian Coordinates
    if m==0 # (Always Regular case)

        # Quadrature point ξ and weight ζ for m=0
        ξ = gamma
        ζ = gamma
    
    # Cylindrical Polar Coordinates
    elseif m==1
        
        # Quadrature point ξ and weight ζ for m=1
        if singular
            ξ = (2/3) .* (α .+ β  .- ((α .* β) ./ (α .+ β)))
        else # Regular case
            ξ = (β .- α) ./ log.(β ./ α)
        end
        ζ = (ξ .* gamma ).^(0.5)
    
    # Spherical Polar Coordinates
    else m==2

        # Quadrature point ξ and weight ζ for m=2
        if singular
            ξ = (2/3) .* (α .+ β  .- ((α .* β) ./ (α .+ β)))
        else # Regular case
            ξ = (α .* β .* log.(β ./ α)) ./ (β .- α)
        end
        ζ = (α .* β .* gamma).^(1/3)
    end

    pb.ξ = ξ
    pb.ζ = ζ


    # Choosing how to initialize the jacobian with sparsity pattern (to review)
    if sparsity == :sparseArray
        pb.jac = spzeros(Nx*npde,Nx*npde)

        for i ∈ 1:npde
            for j ∈ 1:2*npde
                pb.jac[i,j] = 1
            end
        end

        for i ∈ (Nx-1)*npde+1:Nx*npde
            for j ∈ (Nx-2)*npde+1:Nx*npde
                pb.jac[i,j] = 1
            end
        end

        for i ∈ npde+1:npde:(Nx-1)*npde
            for j ∈ i-npde:i+2*npde-1
                pb.jac[i:i+npde-1,j] .= 1
            end
        end
    elseif sparsity == :exSparse # issue with forwarddiff_color_jacobian!
        pb.jac = ExtendableSparseMatrix(Nx*npde,Nx*npde)

        for i ∈ 1:npde
            for j ∈ 1:2*npde
                pb.jac[i,j] = 1
            end
        end

        for i ∈ (Nx-1)*npde+1:Nx*npde
            for j ∈ (Nx-2)*npde+1:Nx*npde
                pb.jac[i,j] = 1
            end
        end

        for i ∈ npde+1:npde:(Nx-1)*npde
            for j ∈ i-npde:i+2*npde-1
                pb.jac[i:i+npde-1,j] .= 1
            end
        end
        flush!(pb.jac)
        pb.jac = pb.jac.cscmatrix

    # elseif sparsity == :symbolics # Sparsity detection for the assemble! fct° using the pkg Symbolics.jl
    #     du0 = copy(inival)
    #     pb.jac = float.(Symbolics.jacobian_sparsity((du,u)->assemble!(du,u,pb,tspan[1]),du0,inival))
    
    elseif sparsity == :banded
        pb.jac = BandedMatrix{Float64}(Ones(Nx*npde,Nx*npde),(2*npde-1,2*npde-1))
    end
    

    # Solve via implicit Euler method or an ODE/DAE solver of DifferentialEquations.jl
    if solver == :euler # implicit Euler method

        colors = matrix_colors(pb.jac)

        # Preallocations for Newton's method
        rhs   = zeros(npde*Nx)
        cache = ForwardColorJacCache(implicitEuler!, rhs, dx=rhs, colorvec=colors, sparsity=pb.jac)

        # To store the history of the Newton solver
        if hist
            storage = []
        end
        
        if tstep==Inf # Solve time independent problems (stationary case) if user set tstep to Inf -> 1 step Implicit Euler
            
            un = inival

            if hist
                unP1, history = newton_stat(un, tstep, 1, pb, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=hist)
            else
                unP1 = newton_stat(un, tstep, 1, pb, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=hist)
            end

            if hist
                return unP1, history
            end

            return unP1
        
        else # Solve time dependent problems (instationary case) via Implicit Euler method

            # Build the mass matrix
            mass_mat_diag = mass_matrix(pb)
            mass_mat_vec  = reshape(mass_mat_diag.diag,(npde,Nx))

            timeSteps = collect(tspan[1]:tstep:tspan[2])
            results   = []
                                    
            un = copy(inival)
            push!(results,reshape(inival,(npde,Nx)))


            for t ∈ timeSteps[2:end]

                lmul!(mass_mat_diag,un)

                if hist
                    unP1, history = newton(un, tstep, t, pb, mass_mat_vec, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=hist)
                else
                    unP1          = newton(un, tstep, t, pb, mass_mat_vec, cache, rhs ; tol=1.0e-10, maxit=100, hist_flag=hist)
                end

                push!(results,reshape(unP1,(npde,Nx)))
                un .= unP1
                
                if hist
                    push!(storage,history)
                end
            end

            if hist
                return results, storage
            end

            return RecursiveArrayTools.DiffEqArray(results,timeSteps)
        end

    elseif solver == :DiffEq # Use the DifferentialEquations.jl package to solve the system of differential equations 

        # returns data from the problem to define an ODEProblem
        return pb
    end

end


# Another way to solve stationary problems 
"""
$(SIGNATURES)

Solve 1D elliptic PDE(s) using the spatial discretization method from Skeel and Berzins
"""
function pdepe(m, pdefun::T1, icfun::T2, bdfun::T3, xmesh ; hist=false, sparsity=:sparseArray) where {T1,T2,T3}

    tspan = (0.0,1.0) # define an arbitrary time interval
    sol = pdepe(m,pdefun,icfun,bdfun,xmesh,tspan ; solver=:euler, tstep=Inf, hist=hist, sparsity=sparsity)

    return sol
end
