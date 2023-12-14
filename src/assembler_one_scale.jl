# Discretization of the problem according to the method introduced in [1].
# Assemble the system of differential equations (operator and mass matrix) for one-scale problems.

"""
    assemble!(du, u, problem, t)

Performs space discretization for one-scale problems following the difference equations described in [1].

Assemble the right-hand side ``f`` to generate an ODE/DAE problem:

```math
(M) \\frac{du}{dt} = f(u,problem,t)
```

where the input `problem` is defined as a [`SkeelBerzins.ProblemDefinitionOneScale`](@ref) structure.

This function is specified in a way that it is compatible with the DifferentialEquations.jl package.
"""
function assemble!(du, u, pb::ProblemDefinitionOneScale{m, npde, singular}, t) where {m, npde, singular}

    # Evaluate the boundary conditions of the problem and interpolate u and du/dx for the first interval of the discretization
    if npde == 1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1], pb.xmesh[end], u[end], t)
        interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1], pb.xmesh[2], u[2], pb.ξ[1], pb)
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1:(pb.npde)], pb.xmesh[end], u[(end - pb.npde + 1):end], t)
        @views interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1:(pb.npde)], pb.xmesh[2],
                                                          u[(pb.npde + 1):(2 * pb.npde)], pb.ξ[1], Val(m),
                                                          Val(singular), Val(npde))
    end
    cl, fl, sl = pb.pdefunction(pb.ξ[1], t, interpolant, d_interpolant)

    # Left boundary of the domain
    for i ∈ 1:npde
        @views assemble_left_bd!(du, u, i, m, pb.xmesh[1], pb.ξ[1], pb.ζ[1], cl[i], fl[i], sl[i], pl[i], ql[i], Val(singular),
                                 Val(pb.markers_macro[1, i]))
    end

    # Interior meshpoints of the domain
    for i ∈ 2:(pb.Nx - 1)
        interpolant, d_interpolant = pb.npde == 1 ?
                                     interpolation(pb.xmesh[i], u[1 + (i - 1) * pb.npde], pb.xmesh[i + 1], u[1 + i * pb.npde],
                                                   pb.ξ[i], pb) :
                                     interpolation(pb.xmesh[i], view(u, (1 + (i - 1) * pb.npde):(i * pb.npde)), pb.xmesh[i + 1],
                                                   view(u, (1 + i * pb.npde):((i + 1) * pb.npde)), pb.ξ[i],
                                                   Val(m), Val(singular), Val(npde))

        cr, fr, sr = pb.pdefunction(pb.ξ[i], t, interpolant, d_interpolant)

        for j ∈ 1:npde
            idx_u = j + (i - 1) * npde
            if !pb.markers_macro[i - 1, j] && pb.markers_macro[i, j]
                @views assemble_left_bd!(du, u, idx_u, pb.m, pb.xmesh[i], pb.ξ[i], pb.ζ[i], cr[j], fr[j], sr[j], pl[j], ql[j],
                                         Val(pb.singular), Val(true))
            elseif pb.markers_macro[i, j] && pb.markers_macro[i + 1, j]
                @views assemble_local!(du, i, idx_u, m, pb.xmesh, pb.ξ, pb.ζ, cl[j], fl[j], sl[j], cr[j], fr[j], sr[j], pl[j],
                                       ql[j], pr[j], qr[j], Val(singular))
            elseif pb.markers_macro[i, j] && !pb.markers_macro[i + 1, j]
                @views assemble_right_bd!(du, u, idx_u, pb.m, pb.xmesh[i], pb.ξ[i], pb.ζ[i], cl[j], fl[j], sl[j], pr[j], qr[j],
                                          Val(pb.singular), Val(true))
            else
                du[idx_u] = u[idx_u]
            end
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    for i ∈ 1:npde
        idx_u = i + (pb.Nx - 1) * npde
        @views assemble_right_bd!(du, u, idx_u, m, pb.xmesh[end], pb.ξ[end], pb.ζ[end], cl[i], fl[i], sl[i], pr[i], qr[i],
                                  Val(pb.singular), Val(pb.markers_macro[pb.Nx, i]))
    end

    nothing
end

"""
    mass_matrix(problem)

Assemble the diagonal mass matrix M of the system of differential equations when solving a one-scale
problem with at least one parabolic PDE. The coefficients from M either take the value 0 or 1 since
it is scaled in the difference equations in the right-hand side.

The entries of the matrix are set to 0 when the corresponding equation of the system is elliptic
or the boundary condition is pure Dirichlet leading to solve a system of DAEs.
In the case where the mass matrix is identity, we solve a system of ODEs.
"""
function mass_matrix(pb::ProblemDefinitionOneScale{m, npde, singular}) where {m, npde, singular}

    inival = pb.inival

    # Initialize the mass matrix M
    M = ones(npde, pb.Nx)
    flag_DAE = false

    if npde == 1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1], pb.xmesh[end], inival[end], pb.tspan[1])
        interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1], pb.xmesh[2], inival[2], pb.ξ[1], pb)
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1:npde], pb.xmesh[end],
                                              inival[(end - npde + 1):end], pb.tspan[1])
        @views interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1:npde], pb.xmesh[2],
                                                          inival[(npde + 1):(2 * npde)], pb.ξ[1],
                                                          Val(m), Val(singular), Val(npde))
    end
    c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)

    for i ∈ 1:npde
        if c[i] == 0 # elliptic equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
            M[i, 2:(end - 1)] .= 0
            flag_DAE = true
        end

        if ql[i] == 0 && !pb.singular # For the left boundary, Dirichlet BC(s) leads to a DAE (ignoring singular case)
            M[i, 1] = 0
            flag_DAE = true
        end

        if qr[i] == 0  # For the right boundary, Dirichlet BC(s) leads to a DAE
            M[i, end] = 0
            flag_DAE = true
        end
    end

    for j ∈ 1:(pb.Nx)
        for i ∈ 1:npde
            if !pb.markers_macro[j, i]
                M[i, j] = 0
                flag_DAE = true
            end
        end
    end

    return Diagonal(vec(M)), flag_DAE
end
