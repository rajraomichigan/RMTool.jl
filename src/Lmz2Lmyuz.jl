"""
    Lmz2Lmyuz(Lmz)

Stieltjes transform polynomial → moment generating function polynomial
`L(myu, z)` where `myu(z) = Σ M_k z^k`, via `z → 1/z`, `m → -myu z`.
"""
function Lmz2Lmyuz(Lmz)
    myu = RMV.myu
    Lmyuz = subs(numden(Lmz), z => _div(1, z), m => -myu * z)
    irreducLuv(Lmyuz, myu, z)
end
