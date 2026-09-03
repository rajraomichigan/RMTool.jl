"""
    AtimesBkernel(LmzA1, LmzA2) -> Lmxz

Kernel polynomial `L(m, x, z)` for the free multiplicative convolution
`A1 × Q A2 Q'`, obtained by eliminating the Cauchy variable between the
polynomials of `A1` and of the product and substituting `f = (1 + 1/(z m))/x`.
"""
function AtimesBkernel(LmzA1, LmzA2)
    g = RMV.g; f = RMV.f; x = RMV.x
    LmzProd = AtimesB(LmzA1, LmzA2)
    LgzProd = Lmz2Lgz(LmzProd)
    LgzA1 = Lmz2Lgz(numden(LmzA1))

    LgfA2 = irreducLuv(subs(LgzA1, g => g * f, z => _inv(f)), g, f)
    LgzProd2 = irreducLuv(subs(LgzProd, g => FF(g) / FF(z)), g, z)

    Lfz = _resultant(LgfA2, LgzProd2, _idx(g))

    Lmxz = subs(Lfz, f => (1 + _inv(z * m)) / FF(x))
    irreducLuv(Lmxz, m, z)
end
