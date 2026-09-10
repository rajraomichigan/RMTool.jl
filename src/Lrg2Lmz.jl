"""
    Lrg2Lmz(Lrg)

R-transform polynomial → Stieltjes transform polynomial (`Lrg2Lgz` followed by `Lgz2Lmz`).
"""
function Lrg2Lmz(Lrg)
    Lgz = Lrg2Lgz(Lrg)
    Lmz = Lgz2Lmz(Lgz)
    irreducLuv(Lmz, m, z)
end
