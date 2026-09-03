"""
    AtimesB(LmzA, LmzB)

Bivariate polynomial for the product of random matrices `C = A × Q*B*Q'` where
`Q` is Haar unitary (free multiplicative convolution).
Assumption: `A × B` has real eigenvalues.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(0 - z))
LmzB = LmzA
LmzC = AtimesB(LmzA, LmzB)
LmzD = transposeA(LmzC, 2)          # should return arc-sine law + atom at zero
```
"""
function AtimesB(LmzA, LmzB)
    s = RMV.s
    LsyA = Lmz2Lsy(LmzA)
    LsyB = Lmz2Lsy(LmzB)
    LsyC = L1timesL2(LsyA, LsyB, s)
    Lsy2Lmz(LsyC)
end
