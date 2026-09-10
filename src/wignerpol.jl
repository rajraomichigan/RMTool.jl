"""
    wignerpol()

Bivariate polynomial `m^2 + m z + 1` encoding the Wigner semicircle law
(the limiting spectral measure of [`wigner`](@ref) matrices, support `[-2, 2]`).
"""
wignerpol() = m^2 + m * z + 1
