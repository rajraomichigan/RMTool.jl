"""
    addAdtimes(LmzA, d)

Bivariate polynomial for the free additive convolution of `d` copies of the
measure encoded by `LmzA` (R-transform scaled by `d`: `R_B = d R_A`).
"""
function addAdtimes(LmzA, d)
    r = RMV.r
    LrgA = Lmz2Lrg(LmzA)
    LrgB = subs(LrgA, r => FF(r) / _ff(d))
    Lrg2Lmz(LrgB)
end
