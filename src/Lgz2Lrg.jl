"""
    Lgz2Lrg(Lgz)

From the Cauchy transform polynomial `L(g, z)` to the polynomial `L(r, g)`
satisfied by the R-transform, using `z = r + 1/g`.
"""
function Lgz2Lrg(Lgz)
    r = RMV.r; g = RMV.g
    Lrg = subs(numden(Lgz), z => r + _div(1, g))
    irreducLuv(Lrg, r, g)
end
