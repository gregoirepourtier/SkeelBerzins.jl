# Discretization of the problem according to the method introduced in [1].
# Assemble the system of differential equations (operator and mass matrix).


"""
    assemble!(du, u, problem, t)

Performs space discretization following the difference equations described in [1].

Assemble the right-hand side ``f`` to generate an ODE/DAE problem:

```math
(M) \\frac{du}{dt} = f(u,problem,t)
```
where the input `problem` is defined as a [`SkeelBerzins.ProblemDefinition`](@ref) structure.

This function is specified in a way that it is compatible with the DifferentialEquations.jl package.
"""
function assemble!(du, u, pb::ProblemDefinition{npde}, t) where {npde}

    type_tmp = eltype(u[1])

    du = reshape(du,(pb.npde,pb.Nx))
    u  = reshape(u,(pb.npde,pb.Nx))

    # Evaluate the boundary conditions of the problem and interpolate u and du/dx for the first 
    # interval of the discretization
    if pb.npde==1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1,1], pb.xmesh[end], u[1,end], t)

        interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1,1], pb.xmesh[2], u[1,2], pb.ξ[1], pb)
        cl, fl, sl  = pb.pdefunction(pb.ξ[1], t, interpolant, d_interpolant)
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[:,1], pb.xmesh[end], u[:,end], t)

        if type_tmp <: AbstractFloat
            interpolant   = pb.interpolant
            d_interpolant = pb.d_interpolant
        else
            interpolant   = Array{type_tmp,1}(undef,pb.npde)
            d_interpolant = Array{type_tmp,1}(undef,pb.npde)
        end
        @views interpolation!(interpolant, d_interpolant, pb.xmesh[1], u[:,1], pb.xmesh[2], u[:,2], pb.ξ[1], pb)
        @views cl, fl, sl  = pb.pdefunction(pb.xmesh[1], t, SVector{npde}(interpolant), SVector{npde}(d_interpolant))
    end

    # Left boundary of the domain
    frac = (pb.ζ[1]^(pb.m+1) - pb.xmesh[1]^(pb.m+1))/(pb.m+1)

    if pb.singular # ignores the given boundary condition to enforce the symmetry condition
        for i ∈ 1:pb.npde
            if cl[i] ≠ 0
                du[i,1] = ((pb.m+1)*fl[i]/pb.ξ[1] + sl[i]) / cl[i]
            else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[i,1] = ((pb.m+1)*fl[i]/pb.ξ[1] + sl[i])
            end
        end
    else # Regular Case
        for i ∈ 1:pb.npde
            if ql[i] ≠ 0 && cl[i] ≠ 0
                du[i,1] = (pl[i] + ql[i]/pb.xmesh[1]^(pb.m) * ((pb.ξ[1]^pb.m)*fl[i] + frac*sl[i])) / 
                          (ql[i]/(pb.xmesh[1]^(pb.m))*frac*cl[i])
            elseif ql[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[i,1] = (pl[i] + ql[i]/pb.xmesh[1]^(pb.m) * ((pb.ξ[1]^pb.m)*fl[i] + frac*sl[i]))
            else # Dirichlet boundary conditions
                du[i,1] = pl[i]
            end
        end
    end

    # Interior meshpoints of the domain
    for i ∈ 2:pb.Nx-1

        if pb.npde==1
            interpolant, d_interpolant = interpolation(pb.xmesh[i], u[1,i], pb.xmesh[i+1], u[1,i+1], pb.ξ[i], pb)
            cr, fr, sr = pb.pdefunction(pb.ξ[i], t, interpolant, d_interpolant)
        else
            @views interpolation!(interpolant, d_interpolant, pb.xmesh[i], u[:,i], pb.xmesh[i+1], u[:,i+1], pb.ξ[i], pb)
            @views cr, fr, sr = pb.pdefunction(pb.xmesh[i], t, SVector{npde}(interpolant), SVector{npde}(d_interpolant))
        end

        frac1 = (pb.ζ[i]^(pb.m+1) - pb.xmesh[i]^(pb.m+1))/(pb.m+1)
        frac2 = (pb.xmesh[i]^(pb.m+1) - pb.ζ[i-1]^(pb.m+1))/(pb.m+1)
        
        if pb.singular
            for j ∈ 1:pb.npde
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[j,i] = (pb.ζ[i]^(pb.m+1)/pb.ξ[i] *fr[j] - pb.ζ[i-1]^(pb.m+1)/pb.ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j]) / 
                              (frac1*cr[j] + frac2*cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j,i] = (pb.ζ[i]^(pb.m+1)/pb.ξ[i] *fr[j] - pb.ζ[i-1]^(pb.m+1)/pb.ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            end
        else # Regular Case
            for j ∈ 1:pb.npde
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[j,i] = (pb.ξ[i]^(pb.m)*fr[j] - pb.ξ[i-1]^(pb.m)*fl[j] + frac1*sr[j] + frac2*sl[j]) / 
                              (frac1*cr[j] + frac2*cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j,i] = (pb.ξ[i]^(pb.m)*fr[j] - pb.ξ[i-1]^(pb.m)*fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            end
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    frac = (pb.xmesh[end]^(pb.m+1) - pb.ζ[end]^(pb.m+1))/(pb.m+1)

    if pb.singular
        for i ∈ 1:pb.npde
            if qr[i] ≠ 0 && cl[i] ≠ 0
                du[i,end] = (pr[i] + qr[i]/pb.xmesh[end]^(pb.m) * (pb.ζ[end]^(pb.m+1)/pb.ξ[end] *fl[i] - frac*sl[i])) / 
                            (-qr[i]/pb.xmesh[end]^(pb.m) * frac*cl[i])
            elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[i,end] = (pr[i] + qr[i]/pb.xmesh[end]^(pb.m) * (pb.ζ[end]^(pb.m+1)/pb.ξ[end] *fl[i] - frac*sl[i]))
            else
                du[i,end] = pr[i]
            end
        end
    else # Regular Case
        for i ∈ 1:pb.npde
            if qr[i] ≠ 0 && cl[i] ≠ 0
                du[i,end] = (pr[i] + qr[i]/pb.xmesh[end]^(pb.m) * (pb.ξ[end]^pb.m *fl[i] - frac*sl[i])) / (-qr[i]/pb.xmesh[end]^(pb.m) * frac*cl[i])
            elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[i,end] = (pr[i] + qr[i]/pb.xmesh[end]^(pb.m) * (pb.ξ[end]^pb.m *fl[i] - frac*sl[i]))
            else
                du[i,end] = pr[i]
            end
        end
    end

    return du
end


"""
    mass_matrix(problem)

Assemble the diagonal mass matrix M of the system of differential equations when solving a problem 
with at least one parabolic PDE. The coefficients from M either take the value 0 or 1 since it is 
scaled in the difference equations in the right-hand side.

The entries of the matrix are set to 0 when the corresponding equation of the system is elliptic
or the boundary condition is pure Dirichlet leading to solve a Differential-Algebraic system of Equations.
In the case where the mass matrix is identity, we solve a system of ODEs.
"""
function mass_matrix(pb::ProblemDefinition{npde}) where {npde}

    # Initialize the mass matrix M
    M = ones(pb.npde, pb.Nx)
    flag_DAE = false

    inival_tmp = reshape(pb.inival,(pb.npde,pb.Nx))

    if pb.npde==1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival_tmp[1,1], pb.xmesh[end], inival_tmp[1,end], pb.tspan[1])

        interpolant, d_interpolant = interpolation(pb.xmesh[1], inival_tmp[1,1], pb.xmesh[2], inival_tmp[1,2], pb.ξ[1], pb)
        c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival_tmp[:,1], pb.xmesh[end], inival_tmp[:,end], pb.tspan[1])

        @views interpolation!(pb.interpolant, pb.d_interpolant, pb.xmesh[1], inival_tmp[:,1], pb.xmesh[2], inival_tmp[:,2], pb.ξ[1], pb)
        @views c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], SVector{npde}(pb.interpolant), SVector{npde}(pb.d_interpolant))
    end

    for i ∈ 1:pb.npde
        if c[i] == 0 # elliptic equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
            M[i,2:end-1] .= 0
            flag_DAE = true
        end

        if ql[i] == 0 && !pb.singular # For the left boundary, Dirichlet BC(s) leads to a DAE (ignoring singular case)
            M[i,1] = 0
            flag_DAE = true
        end

        if qr[i] == 0  # For the right boundary, Dirichlet BC(s) leads to a DAE
            M[i,end] = 0
            flag_DAE = true
        end
    end

    return Diagonal(vec(M)),flag_DAE
end
