"""
    invA(LmzA)

Bivariate polynomial for `B = inv(A)`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = invA(LmzA)
```
"""
invA(LmzA) = mobiusA(LmzA, 0, 1, 1, 0)
