# Discretization of the problem according to the method introduced in [1].
# Assemble the system of differential equations (operator and mass matrix) for two-scale problems.

"""
    assemble!(du, u, problem, t)

Performs space discretization for two-scale problems following the difference equations described in [1].

Assemble the right-hand side ``f`` to generate an ODE/DAE problem:

```math
(M) \\frac{du}{dt} = f(u,problem,t)
```

where the input `problem` is defined as a [`SkeelBerzins.ProblemDefinitionTwoScale`](@ref) structure.

This function is specified in a way that it is compatible with the DifferentialEquations.jl package.
"""
function assemble!(du, u, pb::ProblemDefinitionTwoScale{mx, npde_macro, singular_macro, mr, npde_micro, singular_micro},
                   t) where {mx, npde_macro, singular_macro, mr, npde_micro, singular_micro}

    Nx = pb.Nx
    xmesh = pb.xmesh
    ξ_macro = pb.ξ_macro
    ζ_macro = pb.ζ_macro
    markers_macro = pb.markers_macro

    # Evaluate the boundary conditions of the problem and interpolate u and du/dx for the first interval of the discretization
    if npde_macro == 1
        pl, ql, pr, qr = pb.bdfunction_macro(xmesh[1], u[1], xmesh[end], u[end - pb.Nr], t)
        interpolant, d_interpolant = interpolation(xmesh[1], u[1], xmesh[2], u[2 + pb.Nr], ξ_macro[1], mx, singular_macro)
    else
        @views pl, ql, pr, qr = pb.bdfunction_macro(xmesh[1], u[1:npde_macro], xmesh[end],
                                                    u[(end - pb.Nr - npde_macro + 1):(end - pb.Nr)], t)
        @views interpolant, d_interpolant = interpolation(xmesh[1], u[1:npde_macro], xmesh[2],
                                                          u[(npde_macro + 1 + pb.Nr):(2 * npde_macro + pb.Nr)],
                                                          ξ_macro[1], Val(mx), Val(singular_macro), Val(npde_macro))
    end
    cl, fl, sl = pb.pdefunction_macro(ξ_macro[1], t, interpolant, d_interpolant)

    if pb.markers_micro[1]
        @views sl = pb.coupling_macro(ξ_macro[1], t, interpolant, u[(npde_macro + 1):(npde_macro + pb.Nr)]) # or pb.xmesh[1]
    end

    # Left boundary of the domain
    for i ∈ 1:npde_macro
        @views assemble_left_bd!(du, u, i, mx, xmesh[1], ξ_macro[1], ζ_macro[1], cl[i], fl[i], sl[i], pl[i], ql[i],
                                 Val(singular_macro),
                                 Val(markers_macro[1, i]))
    end

    cpt_marker = 0
    if pb.markers_micro[1]
        idx_u = 1
        idx_uP1 = npde_macro + pb.Nr + 1

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, xmesh[1])

        cpt_marker += 1
    end

    # Interior meshpoints of the domain
    for i ∈ 2:(Nx - 1)
        idx_u = 1 + (i - 1) * npde_macro + cpt_marker * pb.Nr
        idx_uP1 = pb.markers_micro[i] ? idx_u + pb.Nr + npde_macro : idx_u + npde_macro

        interpolant, d_interpolant = npde_macro == 1 ?
                                     interpolation(xmesh[i], u[idx_u], xmesh[i + 1], u[idx_uP1], ξ_macro[i], mx, singular_macro) :
                                     interpolation(xmesh[i], view(u, idx_u:(idx_u + npde_macro - 1)), xmesh[i + 1],
                                                   view(u, idx_uP1:(idx_uP1 + npde_macro - 1)), ξ_macro[i],
                                                   Val(mx), Val(singular_macro), Val(npde_macro))

        cr, fr, sr = pb.pdefunction_macro(ξ_macro[i], t, interpolant, d_interpolant)
        if pb.markers_micro[i]
            @views sr = pb.coupling_macro(ξ_macro[i], t, interpolant, u[(idx_u + 1):(idx_uP1 - 1)])
        end

        for j ∈ 1:npde_macro
            idx = idx_u + j - 1
            if !markers_macro[i - 1, j] && markers_macro[i, j]
                @views assemble_left_bd!(du, u, idx, mx, xmesh[i], ξ_macro[i], ζ_macro[i], cr[j], fr[j], sr[j], pl[j], ql[j],
                                         Val(singular_macro), Val(true))
            elseif markers_macro[i, j] && markers_macro[i + 1, j]
                @views assemble_local!(du, i, idx, mx, xmesh, ξ_macro, ζ_macro, cl[j], fl[j], sl[j], cr[j], fr[j], sr[j], pl[j],
                                       ql[j], pr[j], qr[j], Val(singular_macro))
            elseif markers_macro[i, j] && !markers_macro[i + 1, j]
                @views assemble_right_bd!(du, u, idx, mx, xmesh[i], ξ_macro[i], ζ_macro[i], cl[j], fl[j], sl[j], pr[j], qr[j],
                                          Val(singular_macro), Val(true))
            else
                du[idx] = u[idx]
            end
        end

        if pb.markers_micro[i]
            two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, xmesh[i])
            cpt_marker += 1
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    idx_last = 1 + (Nx - 1) * npde_macro + cpt_marker * pb.Nr
    for i ∈ 1:npde_macro
        idx = idx_last + i - 1
        @views assemble_right_bd!(du, u, idx, mx, xmesh[end], ξ_macro[end], ζ_macro[end], cl[i], fl[i], sl[i], pr[i], qr[i],
                                  Val(singular_macro), Val(markers_macro[Nx, i]))
    end

    if pb.markers_micro[end]
        idx_u = idx_last
        idx_uP1 = idx_last + pb.Nr + npde_macro # doesn't actually exist

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, xmesh[end])
    end

    nothing
end

"""
    mass_matrix(problem)

Assemble the diagonal mass matrix M of the system of differential equations when solving a two-scale
problem with at least one parabolic PDE. The coefficients from M either take the value 0 or 1 since it is
scaled in the difference equations in the right-hand side.

The entries of the matrix are set to 0 when the corresponding equation of the system is elliptic
or the boundary condition is pure Dirichlet leading to solve a Differential-Algebraic system of Equations.
In the case where the mass matrix is identity, we solve a system of ODEs.
"""
function mass_matrix(pb::ProblemDefinitionTwoScale{mx, npde_macro, singular_macro, mr, npde_micro, singular_micro}) where {mx,
                                                                                                                           npde_macro,
                                                                                                                           singular_macro,
                                                                                                                           mr,
                                                                                                                           npde_micro,
                                                                                                                           singular_micro
                                                                                                                           }
    Nx = pb.Nx
    xmesh = pb.xmesh
    ξ_macro = pb.ξ_macro
    markers_macro = pb.markers_macro

    inival = pb.inival

    n_total = npde_macro * Nx + pb.Nx_marked * pb.Nr

    # Initialize the mass matrix M
    M = ones(n_total)
    flag_DAE = false

    if npde_macro == 1
        pl, ql, pr, qr = pb.bdfunction_macro(xmesh[1], inival[1], xmesh[end], inival[end - pb.Nr], pb.tspan[1])
        interpolant, d_interpolant = interpolation(xmesh[1], inival[1], xmesh[2], inival[2 + pb.Nr], ξ_macro[1], mx,
                                                   singular_macro)
    else
        @views pl, ql, pr, qr = pb.bdfunction_macro(xmesh[1], inival[1:npde_macro], xmesh[end],
                                                    inival[(end - pb.Nr - npde_macro + 1):(end - pb.Nr)], pb.tspan[1])
        @views interpolant, d_interpolant = interpolation(xmesh[1], inival[1:npde_macro], xmesh[2],
                                                          inival[(npde_macro + 1 + pb.Nr):(2 * npde_macro + pb.Nr)], ξ_macro[1],
                                                          Val(mx), Val(singular_macro), Val(npde_macro))
    end
    c, f, s = pb.pdefunction_macro(ξ_macro[1], pb.tspan[1], interpolant, d_interpolant)

    # Assume here that npde_micro=1 and c_micro≠0
    pl_micro, ql_micro, pr_micro, qr_micro = pb.bdfunction_micro(pb.rmesh[1], inival[npde_macro + 1], pb.rmesh[end],
                                                                 inival[npde_macro + pb.Nr], pb.tspan[1])

    # interpolant_micro, d_interpolant_micro = interpolation(pb.rmesh[1], inival[pb.npde+1], pb.rmesh[2], inival[pb.npde+2], pb.ξ[1], pb.singular, pb.m)
    # c_micro, f_micro, s_micro = pb.pdefunction_micro(pb.ξ[1], pb.tspan[1], interpolant_micro, d_interpolant_micro)

    for j ∈ 1:npde_macro
        if ql[j] == 0 && !pb.singular_macro # For the left boundary, Dirichlet BC(s) leads to a DAE (ignoring singular case)
            M[j] = 0
            flag_DAE = true
        end
    end

    if pb.markers_micro[1]
        if ql_micro == 0 && !pb.singular_micro
            M[npde_macro + 1] = 0
            flag_DAE = true
        end

        if qr_micro == 0
            M[npde_macro + pb.Nr] = 0
            flag_DAE = true
        end
    end

    i = npde_macro + pb.Nr + 1
    for idx_xmesh ∈ 2:(Nx - 1)
        for j ∈ 1:npde_macro
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

    for j ∈ 1:npde_macro
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
end
