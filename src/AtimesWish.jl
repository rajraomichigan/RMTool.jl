"""
    AtimesWish(LmzA, c)

Bivariate polynomial for a Wishart matrix with covariance `A`:
`B = A × W(c)` where `W(c) = G*G'/N` is Wishart with parameter `c = n/N`,
`G = randn(n, N)`.

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = AtimesWish(LmzA, 1//2)
```
"""
function AtimesWish(LmzA, c)
    cf = _ff(c)
    mf = FF(m); zf = FF(z)
    temp_pol = subs(numden(LmzA), m => mf * (1 - cf - cf * zf * mf),
                                   z => zf / (1 - cf - cf * zf * mf))
    # remove the spurious factor (1 - c - c z m) introduced by the substitution
    spur = numden(1 - cf - cf * zf * mf)
    ok, quo = divides(temp_pol, spur)
    while ok
        temp_pol = quo
        ok, quo = divides(temp_pol, spur)
    end
    irreducLuv(temp_pol, m, z)
end
