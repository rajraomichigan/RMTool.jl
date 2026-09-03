"""
    Lmz2Lgz(Lmz)

Rewrites the bivariate polynomial satisfied by the Stieltjes transform `m(z)`
as the polynomial satisfied by the Cauchy transform `g(z) = -m(z)`.
"""
Lmz2Lgz(Lmz) = subs(numden(Lmz), m => -RMV.g)
