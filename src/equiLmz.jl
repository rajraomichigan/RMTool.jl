"""
    equiLmz(t, M)

Bivariate polynomial for the equilibrium measure corresponding to the
potential `V(x) = t*x^(2M)`.

The constant `a = (M t ∏_{l=1}^{M} (2l-1)/(2l))^(-1/(2M))` is irrational in
general; it is replaced by the closest rational number of its square
(`rationalize`), which is exact when `M = 1`.

```julia
LmzA = equiLmz(1, 1)       # should return a semicircle law
```
"""
function equiLmz(t, M::Integer)
    g = RMV.g
    tq = _rational(rat(t))
    a0 = prod(((2l - 1) // (2l)) for l in 1:M; init = big(1)//1)
    base = a0 * M * tq
    a2 = M == 1 ? 1 // base : rationalize(BigInt, Float64(base)^(-1 / M); tol = 1e-14)
    a2r = rat(a2)

    h1 = z^(2M - 2)
    for j in 1:M-1
        scaling = prod(((2l - 1) // (2l)) for l in 1:j; init = big(1)//1)
        h1 += z^(2M - 2 - 2j) * a2r^j * rat(scaling)
    end

    Lgz = numden((FF(g) / _ff(M * tq) - FF(z^(2M - 1)))^2 - FF((z^2 - a2r) * h1^2))
    Lgz2Lmz(Lgz)
end
