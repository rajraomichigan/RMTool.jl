"""
    AgramWish(LmzA, c, s)

Bivariate polynomial for the Grammian transformation
`B = (X + sqrt(s)*G)*(X + sqrt(s)*G)'` where
`A = X X'`, `G = randn(n, N)` and `c = n/N`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = AgramWish(LmzA, 1//2, 1)
```
"""
function AgramWish(LmzA, c, s)
    cf = _ff(c); sf = _ff(s)
    mf = FF(m); zf = FF(z)
    temp_pol = subs(numden(LmzA), m => mf / (1 + sf * cf * mf),
                                   z => (1 + sf * cf * mf) * (zf * (1 + sf * cf * mf) + sf * (cf - 1)))
    spur = numden(1 + sf * cf * mf)
    ok, quo = divides(temp_pol, spur)
    while ok
        temp_pol = quo
        ok, quo = divides(temp_pol, spur)
    end
    irreducLuv(temp_pol, m, z)
end
