# Local assembly for interior points

"""
    assemble_local!(du, idx_mesh, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{true})

Assemble interior meshpoints of the problem in the "singular" case.
"""
function assemble_local!(du, idx_mesh, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{true})

    frac1 = (ζ[idx_mesh]^(m + 1) - xmesh[idx_mesh]^(m + 1)) / (m + 1)
    frac2 = (xmesh[idx_mesh]^(m + 1) - ζ[idx_mesh - 1]^(m + 1)) / (m + 1)

    if cl ≠ 0 || cr ≠ 0
        du[idx_u] = (ζ[idx_mesh]^(m + 1) / ξ[idx_mesh] * fr - ζ[idx_mesh - 1]^(m + 1) / ξ[idx_mesh - 1] * fl + frac1 * sr +
                     frac2 * sl) / (frac1 * cr + frac2 * cl)
    else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (ζ[idx_mesh]^(m + 1) / ξ[idx_mesh] * fr - ζ[idx_mesh - 1]^(m + 1) / ξ[idx_mesh - 1] * fl + frac1 * sr +
                     frac2 * sl)
    end
end

"""
    assemble_local!(du, idx_mesh, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{false})

Assemble interior meshpoints of the problem in the "regular" case.
"""
function assemble_local!(du, idx_mesh, idx_u, m, xmesh, ξ, ζ, cl, fl, sl, cr, fr, sr, pl, ql, pr, qr, ::Val{false})

    frac1 = (ζ[idx_mesh]^(m + 1) - xmesh[idx_mesh]^(m + 1)) / (m + 1)
    frac2 = (xmesh[idx_mesh]^(m + 1) - ζ[idx_mesh - 1]^(m + 1)) / (m + 1)

    if cl !== zero(cl) || cr !== zero(cl)
        du[idx_u] = (ξ[idx_mesh]^(m) * fr - ξ[idx_mesh - 1]^(m) * fl + frac1 * sr +
                     frac2 * sl) / (frac1 * cr + frac2 * cl)
    else # stationary equation: set the corresponding coefficient in the mass matrix to 0 to generate a DAE
        du[idx_u] = (ξ[idx_mesh]^(m) * fr - ξ[idx_mesh - 1]^(m) * fl + frac1 * sr + frac2 * sl)
    end
end
