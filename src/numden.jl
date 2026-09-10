"""
    numden(expr)

Numerator of a rational expression, i.e. the expression with its denominator
cleared (MATLAB's `numden`). Accepts elements of `R`, of `FF`, or Julia numbers.

```julia
L1 = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))   # atoms at 1 and 2 of equal weight
```
"""
numden(x::FFElem) = numerator(x)
numden(x::RElem)  = x
numden(x::Number) = rat(x)
