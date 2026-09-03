"""
    AplusBkernel(LmzA1, LmzA2) -> (Lmxz, LmzSum)

Kernel polynomial `L(m, x, z)` for the free additive convolution: the resultant
of the Cauchy-transform polynomials of `A1` and of `A1 + A2` eliminates the
common Cauchy variable, and `f` is replaced by `x - 1/m`. The polynomial of the
sum `LmzSum` is returned as well.
"""
function AplusBkernel(LmzA1, LmzA2)
    r = RMV.r; g = RMV.g; f = RMV.f; x = RMV.x
    LA1 = numden(LmzA1); LA2 = numden(LmzA2)

    LgzA1 = Lmz2Lgz(LA1)
    LrgA1 = Lmz2Lrg(LA1)
    LrgA2 = Lmz2Lrg(LA2)
    LrgSum = L1plusL2(LrgA1, LrgA2, r)
    LmzSum = Lrg2Lmz(LrgSum)
    LgzSum = Lmz2Lgz(LmzSum)

    LgfA1 = subs(LgzA1, z => f)
    Lfz = _resultant(LgfA1, LgzSum, _idx(g))

    Lmxz = subs(Lfz, f => x - _inv(m))
    Lmxz = irreducLuv(Lmxz, m, z)
    return Lmxz, LmzSum
end
