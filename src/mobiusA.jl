"""
    mobiusA(LmzA, p, q, r, s)

Bivariate polynomial for the Möbius transformation `B = (p*A + q*I)/(r*A + s*I)`
of the matrix `A` encoded by `LmzA`. The parameters may be numbers or symbolic
parameters (elements of `R`).

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = mobiusA(LmzA, 1, 1, 1, 0)
```
"""
function mobiusA(LmzA, p, q, r, s)
    p, q, r, s = _ff(p), _ff(q), _ff(r), _ff(s)
    zf = FF(z); mf = FF(m)
    alpha = (q - s * zf) / (p - r * zf)
    beta  = _inv(p - r * zf)
    temp_pol = subs(numden(LmzA), z => -alpha)
    temp_pol = subs(temp_pol, m => ((mf / beta) - r) / (s - r * alpha))
    irreducLuv(temp_pol, m, z)
end
