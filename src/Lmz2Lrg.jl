"""
    Lmz2Lrg(Lmz)

Stieltjes transform polynomial → R-transform polynomial (`Lmz2Lgz` followed by `Lgz2Lrg`).
"""
Lmz2Lrg(Lmz) = Lgz2Lrg(Lmz2Lgz(Lmz))
