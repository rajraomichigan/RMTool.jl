"""
    transposeA(LmzA, c)

Bivariate polynomial for `B = X'*X` where `A = X*X'` and
`c = size(A)/size(B)`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(0 - z))
LmzB = transposeA(LmzA, 2)
```
"""
function transposeA(LmzA, c)
    cf = _ff(c)
    LmzB = subs(numden(LmzA), m => (1 - _inv(cf)) * _inv(-FF(z)) + FF(m) / cf)
    irreducLuv(LmzB, m, z)
end
