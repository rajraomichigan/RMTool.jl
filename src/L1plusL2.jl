"""
    L1plusL2(Luv1, Luv2, u; method = :resultant)

Computes the bivariate polynomial `Luv3` whose zeros (with respect to `u`) are
the sums of the zeros of `Luv1` and `Luv2` (with respect to `u`).

Two implementations are provided (see the reference): a matrix-theoretic one
(`method = :matrix`, Kronecker sums of companion matrices) and one in terms of
the resultant `Res_u(L1(u), L2(t - u))` (`method = :resultant`, the default and
much more efficient one).

Reference: N. Raj Rao and Alan Edelman, "The polynomial method for random matrices".
"""
function L1plusL2(Luv1, Luv2, u; method::Symbol = :resultant)
    L1 = numden(Luv1); L2 = numden(Luv2)
    t = RMV.t
    i = _idx(u)
    if method == :matrix
        Cu1 = Luv2Cu(L1, u); Cu2 = Luv2Cu(L2, u)
        n1 = nrows(Cu1); n2 = nrows(Cu2)
        C3 = kronecker_product(Cu1, identity_matrix(FF, n2)) + kronecker_product(identity_matrix(FF, n1), Cu2)
        n = nrows(C3)
        return numden(det(FF(u) * identity_matrix(FF, n) - C3))
    else
        L2t = subs(L2, u => t - u)
        Ltv = _resultant(L1, L2t, i)
        return subs(Ltv, t => u)
    end
end
