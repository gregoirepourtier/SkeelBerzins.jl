# Assemble the micro-scale equations

"""
"""
function two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, idx_micro)

    N        = pb.Nr
    mesh     = pb.rmesh 
    ξ        = pb.ξ_micro
    ζ        = pb.ζ_micro
    singular = pb.singular_micro
    m        = pb.mr
    npde     = pb.npde_micro

    pdefun   = pb.pdefunction_micro
    bdfun    = pb.bdfunction_micro

    coupling_micro = pb.coupling_micro


    pl, ql, pr, qr = bdfun(mesh[1], u[idx_u + pb.npde], mesh[end], u[idx_uP1 - 1], t)
    if pb.npde == 1
        @views pr = coupling_micro(pb.xmesh[idx_micro], t, u[idx_u], u[idx_u+pb.npde : idx_uP1-1])
    else
        @views pr = coupling_micro(pb.xmesh[idx_micro], t, u[idx_u:idx_u+pb.npde-1], u[idx_u+pb.npde : idx_uP1-1])
    end

    interpolant, d_interpolant = interpolation(mesh[1], u[idx_u + pb.npde], mesh[2], u[idx_u + pb.npde + 1], ξ[1], singular, m)
    cl, fl, sl  = pdefun(ξ[1], t, interpolant, d_interpolant)

    # Left boundary of the domain
    frac = (ζ[1]^(m+1) - mesh[1]^(m+1))/(m+1)


    if singular # ignores the given boundary condition to enforce the symmetry condition
        for i ∈ 1:npde
            if cl[i] ≠ 0
                du[idx_u + pb.npde + i - 1] = ((m+1)*fl[i]/ξ[1] + sl[i]) / cl[i]
            else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[idx_u + pb.npde + i - 1] = ((m+1)*fl[i]/ξ[1] + sl[i])
            end
        end
    else # Regular Case
        for i ∈ 1:npde
            if ql[i] ≠ 0 && cl[i] ≠ 0
                du[idx_u + pb.npde + i - 1] = (pl[i] + ql[i]/mesh[1]^(m) * ((ξ[1]^m)*fl[i] + frac*sl[i])) / (ql[i]/(mesh[1]^(m))*frac*cl[i])
            elseif ql[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[idx_u + pb.npde + i - 1] = (pl[i] + ql[i]/mesh[1]^(m) * ((ξ[1]^m)*fl[i] + frac*sl[i]))
            else # Dirichlet boundary conditions
                du[idx_u + pb.npde + i - 1] = pl[i]
            end
        end
    end

    # Interior meshpoints of the domain
    for i ∈ 2:N-1

        index_local = idx_u + pb.npde + i - 1

        interpolant, d_interpolant = interpolation(mesh[i], u[index_local], mesh[i+1], u[index_local+1], ξ[i], singular, m)
        cr, fr, sr = pdefun(ξ[i], t, interpolant, d_interpolant)

        frac1 = (ζ[i]^(m+1) - mesh[i]^(m+1))/(m+1)
        frac2 = (mesh[i]^(m+1) - ζ[i-1]^(m+1))/(m+1)
        
        if singular
            for j ∈ 1:npde
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[index_local] = (ζ[i]^(m+1)/ξ[i] *fr[j] - ζ[i-1]^(m+1)/ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j]) / (frac1*cr[j] + frac2*cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[index_local] = (ζ[i]^(m+1)/ξ[i] *fr[j] - ζ[i-1]^(m+1)/ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            end
        else # Regular Case
            for j ∈ 1:npde
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[index_local] = (ξ[i]^(m)*fr[j] - ξ[i-1]^(m)*fl[j] + frac1*sr[j] + frac2*sl[j]) / (frac1*cr[j] + frac2*cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[index_local] = (ξ[i]^(m)*fr[j] - ξ[i-1]^(m)*fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            end
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    frac = (mesh[end]^(m+1) - ζ[end]^(m+1))/(m+1)

    index_last = idx_uP1 - 1

    if singular
        for i ∈ 1:npde
            if qr[i] ≠ 0 && cl[i] ≠ 0
                du[index_last] = (pr[i] + qr[i]/mesh[end]^(m) * (ζ[end]^(m+1)/ξ[end] *fl[i] - frac*sl[i])) / (-qr[i]/mesh[end]^(m) * frac*cl[i])
            elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[index_last] = (pr[i] + qr[i]/mesh[end]^(m) * (ζ[end]^(m+1)/ξ[end] *fl[i] - frac*sl[i]))
            else
                du[index_last] = pr[i]
            end
        end
    else # Regular Case
        for i ∈ 1:npde
            if qr[i] ≠ 0 && cl[i] ≠ 0
                du[index_last] = (pr[i] + qr[i]/mesh[end]^(m) * (ξ[end]^m *fl[i] - frac*sl[i])) / (-qr[i]/mesh[end]^(m) * frac*cl[i])
            elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                du[index_last] = (pr[i] + qr[i]/mesh[end]^(m) * (ξ[end]^m *fl[i] - frac*sl[i]))
            else
                du[index_last] = pr[i]
            end
        end
    end
end
