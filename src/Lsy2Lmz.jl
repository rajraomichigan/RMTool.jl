"""
    Lsy2Lmz(Lsy)

S-transform polynomial → Stieltjes transform polynomial, using
`s = -m/y` followed by `y = -z m - 1`.
"""
function Lsy2Lmz(Lsy)
    s = RMV.s; y = RMV.y
    Lmz = subs(numden(Lsy), s => -_div(m, y))
    Lmz = subs(Lmz, y => -z * m - 1)
    irreducLuv(Lmz, m, z)
end
