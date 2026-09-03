"""
    corrWish(LmzA, LmzB, c)

Bivariate polynomial for a spatio-temporally correlated Wishart matrix
`C = Y*Y'` where `Y = A^{1/2}*G*B^{1/2}`, `G = randn(n, N)` and `c = n/N`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = LmzA
LmzC = corrWish(LmzA, LmzB, 1//2)
```
"""
function corrWish(LmzA, LmzB, c)
    cc = rat(c)
    LmzW = m * (1 - cc - cc * m * z - z) - 1
    LmzWt = transposeA(LmzW, cc)
    LmzT = AtimesB(LmzWt, LmzB)
    LmzTt = transposeA(LmzT, _inv(cc))
    AtimesB(LmzA, LmzTt)
end
