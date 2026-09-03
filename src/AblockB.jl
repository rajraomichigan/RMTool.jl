"""
    AblockB(LmzA, LmzB, c)

Bivariate polynomial transformation for the block diagonal matrix
`C = diag(A, B)` where `c = size(A)/size(C)`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = LmzA
LmzC = AblockB(LmzA, LmzB, 1//2)     # should return the same value as LmzA
```
"""
function AblockB(LmzA, LmzB, c)
    LA = numden(LmzA); LB = numden(LmzB)
    cf = _ff(c)
    if cf == _ff(1)
        LmzC = LA
    elseif iszero(cf)
        LmzC = LB
    elseif LA == LB
        LmzC = LA
    else
        LmzA1 = subs(LA, m => FF(m) / cf)
        LmzB1 = subs(LB, m => FF(m) / (1 - cf))
        LmzC = L1plusL2(LmzA1, LmzB1, m)
    end
    irreducLuv(LmzC, m, z)
end
