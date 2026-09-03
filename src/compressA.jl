"""
    compressA(LmzA, factor)

Bivariate polynomial for a random compression of a matrix: `Bn` = top
`n × n` block of the `N × N` matrix `Q*A*Q'` where `Q` is Haar unitary and
`factor = n/N` (less than 1).

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(0 - z))
LmzB = compressA(LmzA, 1//2)       # should return the arc-sine law
```
"""
function compressA(LmzA, factor)
    r = RMV.r
    LmzT = scaleA(LmzA, factor)
    LrgT = Lmz2Lrg(LmzT)
    LrgB = subs(LrgT, r => r * _ff(factor))
    Lrg2Lmz(LrgB)
end
