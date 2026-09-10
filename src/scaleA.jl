"""
    scaleA(LmzA, alpha)

Bivariate polynomial for `B = alpha*A`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = scaleA(LmzA, 2)
```
"""
scaleA(LmzA, alpha) = mobiusA(LmzA, alpha, 0, 0, 1)
