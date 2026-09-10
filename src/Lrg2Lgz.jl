"""
    Lrg2Lgz(Lrg)

From the R-transform polynomial `L(r, g)` back to the Cauchy transform
polynomial `L(g, z)`, using `r = z - 1/g`.
"""
function Lrg2Lgz(Lrg)
    r = RMV.r; g = RMV.g
    Lgz = subs(numden(Lrg), r => z - _div(1, g))
    irreducLuv(Lgz, g, z)
end
