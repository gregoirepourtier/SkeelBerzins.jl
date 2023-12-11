# Local assembly for interior points

"""
    assemble_local!(du, u, i, idx_u, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr)

Assemble interior meshpoints of the problem in the "singular" case.
"""
function assemble_local!(du, u, idx_x, idx_u, index_npde, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{true})

    frac1 = (pb.ζ[idx_x]^(pb.m + 1) - pb.xmesh[idx_x]^(pb.m + 1)) / (pb.m + 1)
    frac2 = (pb.xmesh[idx_x]^(pb.m + 1) - pb.ζ[idx_x - 1]^(pb.m + 1)) / (pb.m + 1)

    if !pb.markers_macro[idx_x - 1, index_npde] && pb.markers_macro[idx_x, index_npde]
        @views assemble_left_bd!(du, u, idx_x, idx_u, pb, cr, fr, sr, pl, ql, Val(true), Val(true))
    elseif pb.markers_macro[idx_x, index_npde] && pb.markers_macro[idx_x + 1, index_npde]
        if cl ≠ 0 || cr ≠ 0
            du[idx_u] = (pb.ζ[idx_x]^(pb.m + 1) / pb.ξ[idx_x] * fr -
                         pb.ζ[idx_x - 1]^(pb.m + 1) / pb.ξ[idx_x - 1] * fl +
                         frac1 * sr + frac2 * sl) / (frac1 * cr + frac2 * cl)
        else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
            du[idx_u] = (pb.ζ[idx_x]^(pb.m + 1) / pb.ξ[idx_x] * fr -
                         pb.ζ[idx_x - 1]^(pb.m + 1) / pb.ξ[idx_x - 1] * fl +
                         frac1 * sr + frac2 * sl)
        end
    elseif pb.markers_macro[idx_x, index_npde] && !pb.markers_macro[idx_x + 1, index_npde]
        @views assemble_right_bd!(du, u, idx_x, idx_u, pb, cl, fl, sl, pr, qr, Val(true), Val(true))
    else
        du[idx_u] = u[idx_u]
    end
end

"""
    assemble_local!(du, u, i, idx_u, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr)

Assemble interior meshpoints of the problem in the "regular" case.
"""
function assemble_local!(du, u, idx_x, idx_u, index_npde, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{false})

    frac1 = (pb.ζ[idx_x]^(pb.m + 1) - pb.xmesh[idx_x]^(pb.m + 1)) / (pb.m + 1)
    frac2 = (pb.xmesh[idx_x]^(pb.m + 1) - pb.ζ[idx_x - 1]^(pb.m + 1)) / (pb.m + 1)

    if !pb.markers_macro[idx_x - 1, index_npde] && pb.markers_macro[idx_x, index_npde]
        @views assemble_left_bd!(du, u, idx_x, idx_u, pb, cr, fr, sr, pl, ql, Val(false), Val(true))
    elseif pb.markers_macro[idx_x, index_npde] && pb.markers_macro[idx_x + 1, index_npde]
        if cl !== zero(cl) || cr !== zero(cl)
            du[idx_u] = (pb.ξ[idx_x]^(pb.m) * fr - pb.ξ[idx_x - 1]^(pb.m) * fl + frac1 * sr +
                         frac2 * sl) / (frac1 * cr + frac2 * cl)
        else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
            du[idx_u] = (pb.ξ[idx_x]^(pb.m) * fr - pb.ξ[idx_x - 1]^(pb.m) * fl + frac1 * sr + frac2 * sl)
        end
    elseif pb.markers_macro[idx_x, index_npde] && !pb.markers_macro[idx_x + 1, index_npde]
        @views assemble_right_bd!(du, u, idx_x, idx_u, pb, cl, fl, sl, pr, qr, Val(false), Val(true))
    else
        du[idx_u] = u[idx_u]
    end
end
