"""
    AplusB(LmzA, LmzB)

Bivariate polynomial for the sum of random matrices `C = A + Q*B*Q'` where `Q`
is Haar unitary (free additive convolution).

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(-1 - z))
LmzB = LmzA
LmzC = AplusB(LmzA, LmzB)           # should return the arc-sine law
```
"""
function AplusB(LmzA, LmzB)
    r = RMV.r
    LA = numden(LmzA); LB = numden(LmzB)
    LrgA = Lmz2Lrg(LA)
    LrgB = (LA == LB) ? LrgA : Lmz2Lrg(LB)
    LrgC = L1plusL2(LrgA, LrgB, r)
    Lrg2Lmz(LrgC)
end
