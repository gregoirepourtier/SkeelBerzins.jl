# Assemble boundaries of the problem

"""
    assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, ::Val{true}, ::Val{true})

Assemble the left boundary of the problem in the "singular" case.
"""
@inline function assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, ::Val{true}, ::Val{true})
    if cl !== zero(cl)
        du[idx_u] = ((m + 1) * fl / ξ + sl) / cl
    else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = ((m + 1) * fl / ξ + sl)
    end
end

"""
    assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, ::Val{false}, ::Val{true})

Assemble the left boundary of the problem in the "regular" case.
"""
@inline function assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, ::Val{false}, ::Val{true})

    frac = (ζ^(m + 1) - xmesh^(m + 1)) / (m + 1)

    if ql ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pl + ql / xmesh^m * ((ξ^m) * fl + frac * sl)) / (ql / (xmesh^m) * frac * cl)
    elseif ql ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pl + ql / xmesh^m * ((ξ^m) * fl + frac * sl))
    else # Dirichlet boundary conditions
        du[idx_u] = pl
    end
end

"""
    assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, _, ::Val{false})

Assemble the left boundary of the problem as a "dummy" equation.
"""
@inline function assemble_left_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pl, ql, _, ::Val{false})
    du[idx_u] = u[idx_u]
end

"""
    assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, ::Val{true}, ::Val{true})

Assemble the right boundary of the problem in the "singular" case.
"""
@inline function assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, ::Val{true}, ::Val{true})
    frac = (xmesh^(m + 1) - ζ^(m + 1)) / (m + 1)

    if qr ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pr + qr / xmesh^m * (ζ^(m + 1) / ξ * fl - frac * sl)) / (-qr / xmesh^m * frac * cl)
    elseif qr ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pr + qr / xmesh^m * (ζ^(m + 1) / ξ * fl - frac * sl))
    else
        du[idx_u] = pr
    end
end

"""
    assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, ::Val{false}, ::Val{true})

Assemble the right boundary of the problem in the "regular" case.
"""
@inline function assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, ::Val{false}, ::Val{true})
    frac = (xmesh^(m + 1) - ζ^(m + 1)) / (m + 1)

    if qr ≠ 0 && cl !== zero(cl)
        du[idx_u] = (pr + qr / xmesh^m * (ξ^m * fl - frac * sl)) / (-qr / xmesh^m * frac * cl)
    elseif qr ≠ 0 && cl === zero(cl) # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (pr + qr / xmesh^m * (ξ^m * fl - frac * sl))
    else
        du[idx_u] = pr # index_npde + (idx_mesh - 1) * pb.npde
    end

end

"""
    assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, _, ::Val{false})

Assemble the right boundary of the problem as a "dummy" equation.
"""
@inline function assemble_right_bd!(du, u, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, pr, qr, _, ::Val{false})
    du[idx_u] = u[idx_u]
end
