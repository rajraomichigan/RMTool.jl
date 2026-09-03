"""
    Lmz2Letaz(Lmz)

Stieltjes transform polynomial → eta-transform polynomial `L(eta, z)`, via
`z → -1/z`, `m → z eta`.
"""
function Lmz2Letaz(Lmz)
    eta = RMV.eta
    Letaz = subs(numden(Lmz), z => -_div(1, z), m => z * eta)
    irreducLuv(Letaz, eta, z)
end
