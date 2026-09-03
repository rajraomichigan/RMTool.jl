"""
    Lgz2Lmz(Lgz)

Inverse of [`Lmz2Lgz`](@ref): substitutes `g = -m`.
"""
Lgz2Lmz(Lgz) = subs(numden(Lgz), RMV.g => -m)
