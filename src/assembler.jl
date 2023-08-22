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

    if pb.Nr === nothing
        Nr_tmp = 0
    else
        Nr_tmp = pb.Nr
    end


    # Evaluate the boundary conditions of the problem and interpolate u and du/dx for the first interval of the discretization
    if pb.npde==1
        pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1], pb.xmesh[end], u[end - Nr_tmp], t)

        interpolant, d_interpolant = interpolation(pb.xmesh[1], u[1], pb.xmesh[2], u[2 + Nr_tmp], pb.ξ[1], pb.singular, pb.m)
        cl, fl, sl  = pb.pdefunction(pb.ξ[1], t, interpolant, d_interpolant)
        if pb.Nr !== nothing && pb.markers[1]
            @views sl = pb.coupling_macro(pb.xmesh[1], t, interpolant, u[pb.npde+1 : pb.npde + Nr_tmp])
        end
    else
        @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], u[1 : pb.npde], pb.xmesh[end], u[end - Nr_tmp - pb.npde+1 : end - Nr_tmp], t)

        if type_tmp <: AbstractFloat
            interpolant   = pb.interpolant
            d_interpolant = pb.d_interpolant
        else
            interpolant   = Array{type_tmp,1}(undef,pb.npde)
            d_interpolant = Array{type_tmp,1}(undef,pb.npde)
        end
        @views interpolation!(interpolant, d_interpolant, pb.xmesh[1], u[1 : pb.npde], pb.xmesh[2], u[pb.npde+1 + Nr_tmp : 2*pb.npde + Nr_tmp], pb.ξ[1], pb.singular, pb.m, pb.npde)
        @views cl, fl, sl  = pb.pdefunction(pb.xmesh[1], t, SVector{npde}(interpolant), SVector{npde}(d_interpolant))
        if pb.Nr !== nothing && pb.markers[1]
            @views sl = pb.coupling_macro(pb.xmesh[1], t, SVector{npde}(interpolant), u[pb.npde+1 : pb.npde + Nr_tmp])
        end
    end

    type_check_c = eltype(cl[1])

    # Left boundary of the domain
    frac = (pb.ζ[1]^(pb.m+1) - pb.xmesh[1]^(pb.m+1))/(pb.m+1)

    @views assemble_left_bd!(du,u,1,1,pb,cl,fl,sl,pl,ql,frac,type_check_c)

    cpt_marker = 0

    if (pb.Nr !== nothing) && pb.markers[1]
        idx_u   = 1
        idx_uP1 = pb.npde + Nr_tmp + 1

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, 1)

        cpt_marker += 1 
    end
    
    # Interior meshpoints of the domain
    for i ∈ 2:pb.Nx-1

        idx_u   = 1 + (i-1)*pb.npde + cpt_marker*Nr_tmp
        if pb.Nr !== nothing
            idx_uP1 = pb.markers[i] ? idx_u + Nr_tmp + pb.npde : idx_u + pb.npde
        else
            idx_uP1 = idx_u + pb.npde
        end

        if pb.npde==1
            interpolant, d_interpolant = interpolation(pb.xmesh[i], u[idx_u], pb.xmesh[i+1], u[idx_uP1], pb.ξ[i], pb.singular, pb.m)
            cr, fr, sr = pb.pdefunction(pb.ξ[i], t, interpolant, d_interpolant)
            if pb.Nr !== nothing && pb.markers[i]
                @views sr = pb.coupling_macro(pb.xmesh[i], t, interpolant, u[idx_u+1 : idx_uP1-1])
            end
        else
            @views interpolation!(interpolant, d_interpolant, pb.xmesh[i], u[idx_u : idx_u + pb.npde - 1], pb.xmesh[i+1], u[idx_uP1 : idx_uP1 + pb.npde - 1], pb.ξ[i], pb.singular, pb.m, pb.npde)
            @views cr, fr, sr = pb.pdefunction(pb.xmesh[i], t, SVector{npde}(interpolant), SVector{npde}(d_interpolant))
            if pb.Nr !== nothing && pb.markers[i]
                @views sr = pb.coupling_macro(pb.xmesh[i], t, SVector{npde}(interpolant), u[idx_u+1 : idx_uP1-1])
            end
        end

        frac1 = (pb.ζ[i]^(pb.m+1) - pb.xmesh[i]^(pb.m+1))/(pb.m+1)
        frac2 = (pb.xmesh[i]^(pb.m+1) - pb.ζ[i-1]^(pb.m+1))/(pb.m+1)
        
        @views assemble_local!(du,u,i,idx_u,pb,cl,fl,sl,cr,fr,sr,frac1,frac2,pr,qr,type_check_c)

        if (pb.Nr !== nothing) && pb.markers[i]
            two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, i)
            cpt_marker += 1
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    frac = (pb.xmesh[end]^(pb.m+1) - pb.ζ[end]^(pb.m+1))/(pb.m+1)

    idx_last =  1 + (pb.Nx-1)*pb.npde + cpt_marker*Nr_tmp
    
    
    @views assemble_right_bd!(du,u,pb.Nx,idx_last,pb,cl,fl,sl,pr,qr,frac,type_check_c)
    
    
    # du[end] = u[end] # du[..] .= 0 --> other option and mass matrix to 1
    

    if (pb.Nr !== nothing) && pb.markers[end]
        idx_u   = idx_last
        idx_uP1 = idx_last + Nr_tmp + pb.npde # doesn't actually exist

        two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pb.Nx)
    end

    return du
end


"""
    mass_matrix(problem)

Assemble the diagonal mass matrix M of the system of differential equations 
when solving a problem with at least one parabolic PDE.
The coefficients from M either take the value 0 or 1 since it is scaled in the 
difference equations in the right-hand side.

The entries of the matrix are set to 0 when the corresponding equation of the system is elliptic
or the boundary condition is pure Dirichlet leading to solve a Differential-Algebraic system of Equations.
In the case where the mass matrix is identity, we solve a system of ODEs.
"""
function mass_matrix(pb::ProblemDefinition{npde}) where {npde}

    inival = pb.inival

    if pb.Nr !== nothing
        n_total = pb.npde*pb.Nx + pb.Nx_marked*pb.Nr
        Nr_tmp = pb.Nr

        # Initialize the mass matrix M
        M = ones(n_total)
        flag_DAE = false

        if pb.npde==1
            pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1], pb.xmesh[end], inival[end - Nr_tmp], pb.tspan[1])

            interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1], pb.xmesh[2], inival[2 + Nr_tmp], pb.ξ[1], pb.singular, pb.m)
            c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)
        else
            @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1 : pb.npde], pb.xmesh[end], inival[end - Nr_tmp - pb.npde+1 : end - Nr_tmp], pb.tspan[1])

            @views interpolation!(pb.interpolant, pb.d_interpolant, pb.xmesh[1], inival[1 : pb.npde], pb.xmesh[2], inival[pb.npde+1 + Nr_tmp : 2*pb.npde + Nr_tmp], pb.ξ[1], pb.singular, pb.m, pb.npde)
            @views c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], SVector{npde}(pb.interpolant), SVector{npde}(pb.d_interpolant))
        end

        # assume here that npde_micro=1 and c_micro≠0
        pl_micro, ql_micro, pr_micro, qr_micro = pb.bdfunction_micro(pb.rmesh[1], inival[pb.npde+1], pb.rmesh[end], inival[pb.npde+Nr_tmp], pb.tspan[1])

        # interpolant_micro, d_interpolant_micro = interpolation(pb.rmesh[1], inival[pb.npde+1], pb.rmesh[2], inival[pb.npde+2], pb.ξ[1], pb.singular, pb.m)
        # c_micro, f_micro, s_micro = pb.pdefunction_micro(pb.ξ[1], pb.tspan[1], interpolant_micro, d_interpolant_micro)

        for j ∈ 1:pb.npde
            if ql[j] == 0 && !pb.singular # For the left boundary, Dirichlet BC(s) leads to a DAE (ignoring singular case)
                M[j] = 0
                flag_DAE = true
            end
        end

        if pb.markers[1]
            if ql_micro == 0 && !pb.singular_micro
                M[pb.npde+1] = 0
                flag_DAE = true
            end

            if qr_micro == 0
                M[pb.npde+Nr_tmp] = 0
                flag_DAE = true
            end
        end

        i = pb.npde + pb.Nr + 1
        for idx_xmesh ∈ 2:pb.Nx-1

            for j ∈ 1:pb.npde
                if c[j] == 0 # elliptic equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    M[i] = 0
                    flag_DAE = true
                end
                i += 1
            end
            
            if pb.markers[idx_xmesh]
                if ql_micro == 0 && !pb.singular_micro
                    M[i] = 0
                    flag_DAE = true
                end
                
                i += Nr_tmp - 1 # assume here that c_micro ≠ 0

                if qr_micro == 0
                    M[i] = 0
                    flag_DAE = true
                end
                i += 1
            end
        end

        for j ∈ 1:pb.npde
            if qr[j] == 0  # For the right boundary, Dirichlet BC(s) leads to a DAE
                M[i] = 0
                flag_DAE = true
            end
            i += 1
        end
        
        if pb.markers[end]
            if ql_micro == 0 && !pb.singular_micro
                M[i] = 0
                flag_DAE = true
                i += Nr_tmp - 1
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

        if pb.npde==1
            pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1], pb.xmesh[end], inival[end], pb.tspan[1])

            interpolant, d_interpolant = interpolation(pb.xmesh[1], inival[1], pb.xmesh[2], inival[2], pb.ξ[1], pb.singular, pb.m)
            c, f, s = pb.pdefunction(pb.ξ[1], pb.tspan[1], interpolant, d_interpolant)
        else
            @views pl, ql, pr, qr = pb.bdfunction(pb.xmesh[1], inival[1 : pb.npde], pb.xmesh[end], inival[end-pb.npde+1 : end], pb.tspan[1])

            @views interpolation!(pb.interpolant, pb.d_interpolant, pb.xmesh[1], inival[1 : pb.npde], pb.xmesh[2], inival[pb.npde+1 : 2*pb.npde], pb.ξ[1], pb.singular, pb.m, pb.npde)
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

        for i in 1:pb.Nx
            for j in 1:pb.npde
                if !pb.markers_macro[i,j]
                    # println("hi")
                    M[j,i] = 0
                    flag_DAE = true
                end
            end
        end

        return Diagonal(vec(M)),flag_DAE
    end
end
