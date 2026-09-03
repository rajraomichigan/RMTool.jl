"""
    L1timesL2(Luv1, Luv2, u; method = :resultant)

Computes the bivariate polynomial `Luv3` whose zeros (with respect to `u`) are
the products of the zeros of `Luv1` and `Luv2` (with respect to `u`).

As for [`L1plusL2`](@ref), a matrix-theoretic version (`method = :matrix`,
Kronecker product of companion matrices) and the default resultant version
`Res_u(L1(u), u^deg(L2) L2(t/u))` are available.

Reference: N. Raj Rao and Alan Edelman, "The polynomial method for random matrices".
"""
function L1timesL2(Luv1, Luv2, u; method::Symbol = :resultant)
    L1 = numden(Luv1); L2 = numden(Luv2)
    t = RMV.t
    i = _idx(u)
    if method == :matrix
        Cu1 = Luv2Cu(L1, u); Cu2 = Luv2Cu(L2, u)
        C3 = kronecker_product(Cu1, Cu2)
        n = nrows(C3)
        return numden(det(FF(u) * identity_matrix(FF, n) - C3))
    else
        Du2 = degree(L2, i)
        L2t = numden(FF(u)^Du2 * evaluate(L2, _subsvec(u => FF(t) / FF(u))))
        Ltv = _resultant(L1, L2t, i)
        return subs(Ltv, t => u)
    end
end

# vector of FF values for `evaluate` implementing the given substitutions
function _subsvec(prs::Pair...)
    vals = _ffvec()
    for (v, val) in prs
        vals[_idx(v)] = _ff(val)
    end
    vals
end
