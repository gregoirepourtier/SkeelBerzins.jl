# Local assembly for interior points


"""
"""
function assemble_local!(du, u, i, idx_u, pb, cl, fl, sl, cr, fr, sr, frac1, frac2, pr, qr, type_check_c)
    if pb.singular
        for j ∈ 1:pb.npde
            if pb.markers_macro[i,j] && pb.markers_macro[i+1,j]
                if cl[j] ≠ 0 || cr[j] ≠ 0
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m+1)/pb.ξ[i] *fr[j] - pb.ζ[i-1]^(pb.m+1)/pb.ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j]) / (frac1*cr[j] + frac2*cl[j])
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ζ[i]^(pb.m+1)/pb.ξ[i] *fr[j] - pb.ζ[i-1]^(pb.m+1)/pb.ξ[i-1] *fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            elseif pb.markers_macro[i,j]
                @views assemble_right_bd_singular!(du,i,j,idx_u,pb,cl,fl,sl,pr,qr,frac2,type_check_c)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    else # Regular Case
        for j ∈ 1:pb.npde
            if pb.markers_macro[i,j] && pb.markers_macro[i+1,j]
                if cl[j] !== zero(type_check_c) || cr[j] !== zero(type_check_c)
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m)*fr[j] - pb.ξ[i-1]^(pb.m)*fl[j] + frac1*sr[j] + frac2*sl[j]) / (frac1*cr[j] + frac2*cl[j]) # j + (idx_u-1)*pb.npde
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[j + idx_u - 1] = (pb.ξ[i]^(pb.m)*fr[j] - pb.ξ[i-1]^(pb.m)*fl[j] + frac1*sr[j] + frac2*sl[j])
                end
            elseif pb.markers_macro[i,j]
                @views assemble_right_bd_regular!(du,i,j,idx_u,pb,cl,fl,sl,pr,qr,frac2,type_check_c)
            else
                du[j + idx_u - 1] = u[j + idx_u - 1]
            end
        end
    end
end
