"""
    shiftA(LmzA, alpha)

Bivariate polynomial for `B = A + alpha*I`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = shiftA(LmzA, -1)
```
"""
shiftA(LmzA, alpha) = mobiusA(LmzA, 1, alpha, 0, 1)
