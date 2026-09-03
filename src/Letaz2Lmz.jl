"""
    Letaz2Lmz(Letaz)

Inverse of [`Lmz2Letaz`](@ref): `z → -1/z`, `eta → -z m`.
"""
function Letaz2Lmz(Letaz)
    eta = RMV.eta
    Lmz = subs(numden(Letaz), z => -_div(1, z), eta => -z * m)
    irreducLuv(Lmz, m, z)
end
