# Assemble boundaries of the problem

"""
"""
function assemble_left_bd!(du, u, idx_x, idx_u, pb, cl, fl, sl, pl, ql, frac, type_check_c)
    if pb.singular # ignores the given boundary condition to enforce the symmetry condition
        for i ∈ 1:pb.npde
            if pb.markers_macro[idx_x,i]
                if cl[i] ≠ 0
                    du[i + idx_u - 1] = ((pb.m+1)*fl[i]/pb.ξ[idx_x] + sl[i]) / cl[i]
                else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[i + idx_u - 1] = ((pb.m+1)*fl[i]/pb.ξ[idx_x] + sl[i])
                end
            else
                du[i + idx_u - 1] = u[i + idx_u - 1]
            end
        end
    else # Regular Case
        for i ∈ 1:pb.npde
            if pb.markers_macro[idx_x,i]
                if ql[i] ≠ 0 && cl[i] !== zero(type_check_c)
                    du[i + idx_u - 1] = (pl[i] + ql[i]/pb.xmesh[idx_x]^(pb.m) * ((pb.ξ[idx_x]^pb.m)*fl[i] + frac*sl[i])) / (ql[i]/(pb.xmesh[idx_x]^(pb.m))*frac*cl[i])
                elseif ql[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[i + idx_u - 1] = (pl[i] + ql[i]/pb.xmesh[idx_x]^(pb.m) * ((pb.ξ[idx_x]^pb.m)*fl[i] + frac*sl[i]))
                else # Dirichlet boundary conditions
                    du[i + idx_u - 1] = pl[i]
                end
            else
                du[i + idx_u - 1] = u[i + idx_u - 1]
            end
        end
    end
end

"""
"""
function assemble_right_bd!(du, u, idx_x, idx_u, pb, cl, fl, sl, pr, qr, frac, type_check_c)
    idx_quad = idx_x==pb.Nx ? idx_x-1 : idx_x

    if pb.singular
        for i ∈ 1:pb.npde
            if pb.markers_macro[idx_x,i]
                if qr[i] ≠ 0 && cl[i] ≠ 0
                    du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ζ[idx_quad]^(pb.m+1)/pb.ξ[idx_quad] *fl[i] - frac*sl[i])) / (-qr[i]/pb.xmesh[idx_x]^(pb.m) * frac*cl[i]) # i + (idx_u-1)*pb.npde
                elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ζ[idx_quad]^(pb.m+1)/pb.ξ[idx_quad] *fl[i] - frac*sl[i]))
                else
                    du[i + idx_u - 1] = pr[i]
                end
            else
                du[i + idx_u - 1] = u[i + idx_u - 1]
            end
        end
    else # Regular Case
        for i ∈ 1:pb.npde
            if pb.markers_macro[idx_x,i]
                if qr[i] ≠ 0 && cl[i] !== zero(type_check_c)
                    du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ξ[idx_quad]^pb.m *fl[i] - frac*sl[i])) / (-qr[i]/pb.xmesh[idx_x]^(pb.m) * frac*cl[i])
                elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
                    du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ξ[idx_quad]^pb.m *fl[i] - frac*sl[i]))
                else
                    du[i + idx_u - 1] = pr[i]
                end
            else
                du[i + idx_u - 1] = u[i + idx_u - 1]
            end
        end
    end
end


"""
"""
function assemble_right_bd_regular!(du, idx_x, i, idx_u, pb, cl, fl, sl, pr, qr, frac, type_check_c)
    # idx_quad = idx_x==pb.Nx ? idx_x-1 : idx_x
    
    if qr[i] ≠ 0 && cl[i] !== zero(type_check_c)
        du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ξ[idx_x]^pb.m *fl[i] - frac*sl[i])) / (-qr[i]/pb.xmesh[idx_x]^(pb.m) * frac*cl[i])
    elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ξ[idx_x]^pb.m *fl[i] - frac*sl[i]))
    else
        du[i + idx_u - 1] = pr[i]
    end
end

"""
"""
function assemble_right_bd_singular!(du, idx_x, i, idx_u, pb, cl, fl, sl, pr, qr, frac, type_check_c)
    # idx_quad = idx_x==pb.Nx ? idx_x-1 : idx_x

    if qr[i] ≠ 0 && cl[i] ≠ 0
        du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ζ[idx_x]^(pb.m+1)/pb.ξ[idx_x] *fl[i] - frac*sl[i])) / (-qr[i]/pb.xmesh[idx_x]^(pb.m) * frac*cl[i])
    elseif qr[i] ≠ 0 && cl[i] == 0 # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[i + idx_u - 1] = (pr[i] + qr[i]/pb.xmesh[idx_x]^(pb.m) * (pb.ζ[idx_x]^(pb.m+1)/pb.ξ[idx_x] *fl[i] - frac*sl[i]))
    else
        du[i + idx_u - 1] = pr[i]
    end
end
