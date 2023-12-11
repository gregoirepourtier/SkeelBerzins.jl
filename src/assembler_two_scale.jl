# Discretization of the problem according to the method introduced in [1].
# Assemble the system of differential equations (operator and mass matrix) for two-scale problems.

"""
    assemble_two_scale!(du, u, problem, t)

Performs space discretization for two-scale problems following the difference equations described in [1].

Assemble the right-hand side ``f`` to generate an ODE/DAE problem:

```math
(M) \\frac{du}{dt} = f(u,problem,t)
```

where the input `problem` is defined as a [`SkeelBerzins.ProblemDefinition`](@ref) structure.

This function is specified in a way that it is compatible with the DifferentialEquations.jl package.
"""
function assemble_two_scale!(du, u, pb::ProblemDefinition{m, npde, singular}, t) where {m, npde, singular}

    # Evaluate the boundary conditions of the problem and interpolate u and du/dx for the first interval of the discretization
    if npde == 1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1], pb.xmesh[end], u[end - pb.Nr], t)
        interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1], pb.xmesh[2], u[2 + pb.Nr], pb.ξ[1], pb)
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1:(pb.npde)], pb.xmesh[end],
                                              u[(end - pb.Nr - pb.npde + 1):(end - pb.Nr)], t)
        @views interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1:(pb.npde)], pb.xmesh[2],
                                                          u[(pb.npde + 1 + pb.Nr):(2 * pb.npde + pb.Nr)], pb.ξ[1], Val(m),
                                                          Val(singular), Val(npde))
    end
    cl, fl, sl = pb.pdefunction(pb.ξ[1], t, interpolant, d_interpolant)

    if pb.Nr ≠ 0 && pb.markers_micro[1]
        @views sl = pb.coupling_macro(pb.ξ[1], t, interpolant, u[(pb.npde + 1):(pb.npde + pb.Nr)]) # or pb.xmesh[1]
    end

    # Left boundary of the domain
    @views assemble_left_bd!(du, u, 1, 1, pb, cl, fl, sl, pl, ql)

    cpt_marker = 0
    if pb.Nr ≠ 0 && pb.markers_micro[1]
        idx_u = 1
        idx_uP1 = pb.npde + pb.Nr + 1

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pb.xmesh[1])

        cpt_marker += 1
    end

    # Interior meshpoints of the domain
    for i ∈ 2:(pb.Nx - 1)
        idx_u = 1 + (i - 1) * pb.npde + cpt_marker * pb.Nr
        if pb.Nr != 0
            idx_uP1 = pb.markers_micro[i] ? idx_u + pb.Nr + pb.npde : idx_u + pb.npde
        else
            idx_uP1 = idx_u + pb.npde
        end

        interpolant, d_interpolant = npde == 1 ?
                                     interpolation(pb.xmesh[i], u[idx_u], pb.xmesh[i + 1], u[idx_uP1], pb.ξ[i], pb) :
                                     interpolation(pb.xmesh[i], view(u, idx_u:(idx_u + pb.npde - 1)), pb.xmesh[i + 1],
                                                   view(u, idx_uP1:(idx_uP1 + pb.npde - 1)), pb.ξ[i],
                                                   Val(m), Val(singular), Val(npde))

        cr, fr, sr = pb.pdefunction(pb.ξ[i], t, interpolant, d_interpolant)
        if pb.Nr ≠ 0 && pb.markers_micro[i]
            @views sr = pb.coupling_macro(pb.ξ[i], t, interpolant, u[(idx_u + 1):(idx_uP1 - 1)])
        end

        @views assemble_local!(du, u, i, idx_u, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr)

        if pb.Nr ≠ 0 && pb.markers_micro[i]
            two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pb.xmesh[i])
            cpt_marker += 1
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    idx_last = 1 + (pb.Nx - 1) * pb.npde + cpt_marker * pb.Nr
    @views assemble_right_bd!(du, u, pb.Nx, idx_last, pb, cl, fl, sl, pr, qr)

    if pb.Nr ≠ 0 && pb.markers_micro[end]
        idx_u = idx_last
        idx_uP1 = idx_last + pb.Nr + pb.npde # doesn't actually exist

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pb.xmesh[end])
    end

    nothing
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
function mass_matrix_two_scale(pb::ProblemDefinition{m, npde, singular}) where {m, npde, singular}

    inival = pb.inival

    if pb.Nr ≠ 0
        n_total = pb.npde * pb.Nx + pb.Nx_marked * pb.Nr

        # Initialize the mass matrix M
        M = ones(n_total)
        flag_DAE = false

        if pb.npde == 1
            pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1], pb.xmesh[end], inival[end - pb.Nr], pb.tspan[1])
            interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1], pb.xmesh[2], inival[2 + pb.Nr], pb.ξ[1], pb)
        else
            @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1:(pb.npde)], pb.xmesh[end],
                                                  inival[(end - pb.Nr - pb.npde + 1):(end - pb.Nr)], pb.tspan[1])
            @views interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1:(pb.npde)], pb.xmesh[2],
                                                              inival[(pb.npde + 1 + pb.Nr):(2 * pb.npde + pb.Nr)], pb.ξ[1],
                                                              Val(m), Val(singular), Val(npde))
        end
        c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)

        # Assume here that npde_micro=1 and c_micro≠0
        pl_micro, ql_micro, pr_micro, qr_micro = pb.bdfunction_micro(pb.rmesh[1], inival[pb.npde + 1], pb.rmesh[end],
                                                                     inival[pb.npde + pb.Nr], pb.tspan[1])

        # interpolant_micro, d_interpolant_micro = interpolation(pb.rmesh[1], inival[pb.npde+1], pb.rmesh[2], inival[pb.npde+2], pb.ξ[1], pb.singular, pb.m)
        # c_micro, f_micro, s_micro = pb.pdefunction_micro(pb.ξ[1], pb.tspan[1], interpolant_micro, d_interpolant_micro)

        for j ∈ 1:(pb.npde)
            if ql[j] == 0 && !pb.singular # For the left boundary, Dirichlet BC(s) leads to a DAE (ignoring singular case)
                M[j] = 0
                flag_DAE = true
            end
        end

        if pb.markers_micro[1]
            if ql_micro == 0 && !pb.singular_micro
                M[pb.npde + 1] = 0
                flag_DAE = true
            end

            if qr_micro == 0
                M[pb.npde + pb.Nr] = 0
                flag_DAE = true
            end
        end

        i = pb.npde + pb.Nr + 1
        for idx_xmesh ∈ 2:(pb.Nx - 1)
            for j ∈ 1:(pb.npde)
                if c[j] == 0 # elliptic equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    M[i] = 0
                    flag_DAE = true
                end
                i += 1
            end

            if pb.markers_micro[idx_xmesh]
                if ql_micro == 0 && !pb.singular_micro
                    M[i] = 0
                    flag_DAE = true
                end

                i += pb.Nr - 1 # assume here that c_micro ≠ 0

                if qr_micro == 0
                    M[i] = 0
                    flag_DAE = true
                end
                i += 1
            end
        end

        for j ∈ 1:(pb.npde)
            if qr[j] == 0  # For the right boundary, Dirichlet BC(s) leads to a DAE
                M[i] = 0
                flag_DAE = true
            end
            i += 1
        end

        if pb.markers_micro[end]
            if ql_micro == 0 && !pb.singular_micro
                M[i] = 0
                flag_DAE = true
                i += pb.Nr - 1
            end

            if qr_micro == 0
                M[i] = 0
                flag_DAE = true
            end
        end

        return Diagonal(M), flag_DAE
    else
        # Initialize the mass matrix M
        M = ones(pb.npde, pb.Nx)
        flag_DAE = false

        if npde == 1
            pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1], pb.xmesh[end], inival[end], pb.tspan[1])
            interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1], pb.xmesh[2], inival[2], pb.ξ[1], pb)
        else
            @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1:(pb.npde)], pb.xmesh[end],
                                                  inival[(end - pb.npde + 1):end], pb.tspan[1])
            @views interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1:(pb.npde)], pb.xmesh[2],
                                                              inival[(pb.npde + 1):(2 * pb.npde)], pb.ξ[1],
                                                              Val(m), Val(singular), Val(npde))
        end
        c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)

        for i ∈ 1:(pb.npde)
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

            for j ∈ 1:(pb.Nx)
                if !pb.markers_macro[j, i]
                    M[i, j] = 0
                    flag_DAE = true
                end
            end
        end

        return Diagonal(vec(M)), flag_DAE
    end
end
