# Local assembly for interior points

"""
"""
function assemble_local!(du::T8, u::T9, i, idx_u, pb::ProblemDefinition{m, npde, singular}, cl::T1, fl::T2, sl::T3, cr::T4,
                         fr::T5, sr::T6, frac1, frac2, pl, ql, pr, qr,
                         type_check_c::T7) where {T1, T2, T3, T4, T5, T6, T7, m, npde, singular, T8, T9}
    if pb.singular
        for j ∈ 1:(pb.npde)
            if !pb.markers_macro[i - 1, j] && pb.markers_macro[i, j]
                @views assemble_left_bd_singular!(du, i, j, idx_u, pb, cr, fr, sr, type_check_c)
            elseif pb.markers_macro[i, j] && pb.markers_macro[i + 1, j]
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m + 1) / pb.ξ[i] * fr[j] - pb.ζ[i - 1]^(pb.m + 1) / pb.ξ[i - 1] * fl[j] +
                                         frac1 * sr[j] + frac2 * sl[j]) / (frac1 * cr[j] + frac2 * cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m + 1) / pb.ξ[i] * fr[j] - pb.ζ[i - 1]^(pb.m + 1) / pb.ξ[i - 1] * fl[j] +
                                         frac1 * sr[j] + frac2 * sl[j])
                end
            elseif pb.markers_macro[i, j] && !pb.markers_macro[i + 1, j]
                @views assemble_right_bd_singular!(du, i, j, idx_u, i, pb, cl, fl, sl, pr, qr, frac2, type_check_c)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    else # Regular Case
        for j ∈ 1:(pb.npde)
            if !pb.markers_macro[i - 1, j] && pb.markers_macro[i, j]
                @views assemble_left_bd_regular!(du, i, j, idx_u, pb, cr, fr, sr, pl, ql, frac2, type_check_c)
            elseif pb.markers_macro[i, j] && pb.markers_macro[i + 1, j]
                if cl[j] !== zero(type_check_c) || cr[j] !== zero(type_check_c)
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m) * fr[j] - pb.ξ[i - 1]^(pb.m) * fl[j] + frac1 * sr[j] + frac2 * sl[j]) /
                                        (frac1 * cr[j] + frac2 * cl[j]) # j + (idx_u-1)*pb.npde
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m) * fr[j] - pb.ξ[i - 1]^(pb.m) * fl[j] + frac1 * sr[j] + frac2 * sl[j])
                end
            elseif pb.markers_macro[i, j] && !pb.markers_macro[i + 1, j]
                @views assemble_right_bd_regular!(du, i, j, idx_u, i, pb, cl, fl, sl, pr, qr, frac2, type_check_c)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    end
end