"""
    Luv2Cu(Luv, u)

Companion matrix (over the fraction field `FF`) of the bivariate polynomial
`Luv` regarded as a polynomial in `u`: its eigenvalues are the zeros of `Luv`
with respect to `u`. Used by the matrix-theoretic versions of
[`L1plusL2`](@ref) and [`L1timesL2`](@ref).
"""
function Luv2Cu(Luv, u)
    L = numden(Luv)
    i = _idx(u)
    cs = _ucoeffs(L, i)
    Du = length(cs) - 1
    Du >= 1 || error("Luv2Cu: polynomial has degree 0 in $u")
    lead = FF(cs[end])
    Cu = zero_matrix(FF, Du, Du)
    for k in 1:Du-1
        Cu[k, k+1] = one(FF)
    end
    for Di in 0:Du-1
        Cu[Du, Di+1] = -FF(cs[Di+1]) / lead
    end
    Cu
end
