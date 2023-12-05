# Local assembly for interior points

"""
"""
function assemble_local!(du, u, i, idx_u, pb, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr)

    frac1 = (pb.ζ[i]^(pb.m + 1) - pb.xmesh[i]^(pb.m + 1)) / (pb.m + 1)
    frac2 = (pb.xmesh[i]^(pb.m + 1) - pb.ζ[i - 1]^(pb.m + 1)) / (pb.m + 1)

    if pb.singular
        for j ∈ 1:(pb.npde)
            if !pb.markers_macro[i - 1, j] && pb.markers_macro[i, j]
                @views assemble_left_bd_singular!(du, i, j, idx_u, pb, cr, fr, sr)
            elseif pb.markers_macro[i, j] && pb.markers_macro[i + 1, j]
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m + 1) / pb.ξ[i] * fr[j] - pb.ζ[i - 1]^(pb.m + 1) / pb.ξ[i - 1] * fl[j] +
                                         frac1 * sr[j] + frac2 * sl[j]) / (frac1 * cr[j] + frac2 * cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m + 1) / pb.ξ[i] * fr[j] - pb.ζ[i - 1]^(pb.m + 1) / pb.ξ[i - 1] * fl[j] +
                                         frac1 * sr[j] + frac2 * sl[j])
                end
            elseif pb.markers_macro[i, j] && !pb.markers_macro[i + 1, j]
                @views assemble_right_bd_singular!(du, i, j, idx_u, i, pb, cl, fl, sl, pr, qr, frac2)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    else # Regular Case
        for j ∈ 1:(pb.npde)
            if !pb.markers_macro[i - 1, j] && pb.markers_macro[i, j]
                @views assemble_left_bd_regular!(du, i, j, idx_u, pb, cr, fr, sr, pl, ql, frac2)
            elseif pb.markers_macro[i, j] && pb.markers_macro[i + 1, j]
                if cl[j] !== zero(cl[j]) || cr[j] !== zero(cl[j])
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m) * fr[j] - pb.ξ[i - 1]^(pb.m) * fl[j] + frac1 * sr[j] + frac2 * sl[j]) /
                                        (frac1 * cr[j] + frac2 * cl[j]) # j + (idx_u-1)*pb.npde
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m) * fr[j] - pb.ξ[i - 1]^(pb.m) * fl[j] + frac1 * sr[j] + frac2 * sl[j])
                end
            elseif pb.markers_macro[i, j] && !pb.markers_macro[i + 1, j]
                @views assemble_right_bd_regular!(du, i, j, idx_u, i, pb, cl, fl, sl, pr, qr, frac2)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    end
end
