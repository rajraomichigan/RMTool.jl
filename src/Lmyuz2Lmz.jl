"""
    Lmyuz2Lmz(Lmyuz)

Inverse of [`Lmz2Lmyuz`](@ref): `z → 1/z`, `myu → -m z`.
"""
function Lmyuz2Lmz(Lmyuz)
    myu = RMV.myu
    Lmz = subs(numden(Lmyuz), z => _div(1, z), myu => -m * z)
    irreducLuv(Lmz, m, z)
end
