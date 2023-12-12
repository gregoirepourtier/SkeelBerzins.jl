# Assemble the micro-scale equations (working for 1 specie at the moment)

"""
    two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pt_xmesh)

Assemble the system of equations corresponding to the macro-scale.
"""
function two_scale_assembler!(du, u, pb, t, idx_u, idx_uP1, pt_xmesh)

    N = pb.Nr
    mesh = pb.rmesh
    ξ = pb.ξ_micro
    ζ = pb.ζ_micro
    singular = pb.singular_micro
    m = pb.mr
    npde = pb.npde_micro

    pdefun = pb.pdefunction_micro
    bdfun = pb.bdfunction_micro

    coupling_micro = pb.coupling_micro

    pl, ql, pr, qr = bdfun(mesh[1], u[idx_u + pb.npde], mesh[end], u[idx_uP1 - 1], t)
    @views pr = pb.npde == 1 ? coupling_micro(pt_xmesh, t, u[idx_u], u[(idx_u + pb.npde):(idx_uP1 - 1)]) :
                coupling_micro(pt_xmesh, t, u[idx_u:(idx_u + pb.npde - 1)], u[(idx_u + pb.npde):(idx_uP1 - 1)])

    interpolant, d_interpolant = interpolation(mesh[1], u[idx_u + pb.npde], mesh[2], u[idx_u + pb.npde + 1], ξ[1], pb)
    cl, fl, sl = pdefun(ξ[1], t, interpolant, d_interpolant)

    # Left boundary of the domain
    for i ∈ 1:npde
        idx = idx_u + pb.npde + i - 1
        @views assemble_left_bd!(du, u, idx, m, mesh[1], ξ[1], ζ[1], cl[i], fl[i], sl[i], pl[i], ql[i], Val(singular),
                                 Val(true))
    end

    # Interior meshpoints of the domain
    for i ∈ 2:(N - 1)
        index_local = idx_u + pb.npde + i - 1

        interpolant, d_interpolant = interpolation(mesh[i], u[index_local], mesh[i + 1], u[index_local + 1], ξ[i], pb)
        cr, fr, sr = pdefun(ξ[i], t, interpolant, d_interpolant)

        for j ∈ 1:npde # npde = 1
            @views assemble_local!(du, i, index_local, m, mesh, ξ, ζ, cl[j], fl[j], sl[j], cr[j], fr[j], sr[j], pl[j], ql[j],
                                   pr[j], qr[j], Val(singular))
        end

        cl = cr
        fl = fr
        sl = sr
    end

    # Right boundary of the domain
    index_last = idx_uP1 - npde
    for i ∈ 1:npde
        idx = index_last + i - 1
        @views assemble_right_bd!(du, u, idx, m, mesh[end], ξ[end], ζ[end], cl[i], fl[i], sl[i], pr[i], qr[i], Val(singular),
                                  Val(true))
    end

    nothing
end
