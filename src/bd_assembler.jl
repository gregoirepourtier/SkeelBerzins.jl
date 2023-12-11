# Assemble boundaries of the problem

"""
    assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, ::Val{true}, ::Val{true})

Assemble the left boundary of the problem in the "singular" case.
"""
@inline function assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, ::Val{true}, ::Val{true})
    if cl !== zero(cl)
        du[idx_u] = ((pb.m + 1) * fl / pb.ξ[idx_mesh] + sl) / cl
    else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = ((pb.m + 1) * fl / pb.ξ[idx_mesh] + sl)
    end
end

"""
    assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, ::Val{false}, ::Val{true})

Assemble the left boundary of the problem in the "regular" case.
"""
@inline function assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, ::Val{false}, ::Val{true})

    frac = (pb.ζ[idx_mesh]^(pb.m + 1) - pb.xmesh[idx_mesh]^(pb.m + 1)) / (pb.m + 1)

    if ql ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pl + ql / pb.xmesh[idx_mesh]^(pb.m) * ((pb.ξ[idx_mesh]^pb.m) * fl + frac * sl)) /
                    (ql / (pb.xmesh[idx_mesh]^(pb.m)) * frac * cl)
    elseif ql ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pl + ql / pb.xmesh[idx_mesh]^(pb.m) * ((pb.ξ[idx_mesh]^pb.m) * fl + frac * sl))
    else # Dirichlet boundary conditions
        du[idx_u] = pl
    end
end

"""
    assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, _, ::Val{false})

Assemble the left boundary of the problem as a "dummy" equation.
"""
@inline function assemble_left_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pl, ql, _, ::Val{false})
    du[idx_u] = u[idx_u]
end

"""
    assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, ::Val{true}, ::Val{true})

Assemble the right boundary of the problem in the "singular" case.
"""
@inline function assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, ::Val{true}, ::Val{true})
    idx_quad = idx_mesh == pb.Nx ? idx_mesh - 1 : idx_mesh

    frac = (pb.xmesh[idx_mesh]^(pb.m + 1) - pb.ζ[idx_quad]^(pb.m + 1)) / (pb.m + 1)

    if qr ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pr +
                     qr / pb.xmesh[idx_mesh]^(pb.m) * (pb.ζ[idx_quad]^(pb.m + 1) / pb.ξ[idx_quad] * fl - frac * sl)) /
                    (-qr / pb.xmesh[idx_mesh]^(pb.m) * frac * cl)
    elseif qr ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pr + qr / pb.xmesh[idx_mesh]^(pb.m) * (pb.ζ[idx_quad]^(pb.m + 1) / pb.ξ[idx_quad] * fl - frac * sl))
    else
        du[idx_u] = pr
    end
end

"""
    assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, ::Val{false}, ::Val{true})

Assemble the right boundary of the problem in the "regular" case.
"""
@inline function assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, ::Val{false}, ::Val{true})
    idx_quad = idx_mesh == pb.Nx ? idx_mesh - 1 : idx_mesh

    frac = (pb.xmesh[idx_mesh]^(pb.m + 1) - pb.ζ[idx_quad]^(pb.m + 1)) / (pb.m + 1)

    if qr ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pr + qr / pb.xmesh[idx_mesh]^(pb.m) * (pb.ξ[idx_quad]^pb.m * fl - frac * sl)) /
                    (-qr / pb.xmesh[idx_mesh]^(pb.m) * frac * cl)
    elseif qr ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pr + qr / pb.xmesh[idx_mesh]^(pb.m) * (pb.ξ[idx_quad]^pb.m * fl - frac * sl))
    else
        du[idx_u] = pr # index_npde + (idx_mesh - 1) * pb.npde
    end

end

"""
    assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, _, ::Val{false})

Assemble the right boundary of the problem as a "dummy" equation.
"""
@inline function assemble_right_bd!(du, u, idx_mesh, idx_u, pb, cl, fl, sl, pr, qr, _, ::Val{false})
    du[idx_u] = u[idx_u]
end
