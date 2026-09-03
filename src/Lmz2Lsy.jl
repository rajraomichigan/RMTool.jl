"""
    Lmz2Lsy(Lmz)

Stieltjes transform polynomial → S-transform polynomial `L(s, y)`, using
`m = -y s` and `z = (y + 1)/(y s)`.
"""
function Lmz2Lsy(Lmz)
    s = RMV.s; y = RMV.y
    Lsy = subs(numden(Lmz), m => -y * s)
    Lsy = subs(Lsy, z => _div(y + 1, y * s))
    irreducLuv(Lsy, s, y)
end
