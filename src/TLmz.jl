"""
    TLmz(Lmz) -> Matrix

Returns a `(Dm+1) × (Dz+1)` matrix of coefficients of the bivariate polynomial
`Lmz`, whose `(i, j)`-th entry is the coefficient of the `m^(i-1) z^(j-1)` term.
Entries are elements of `R` (they may still contain parameter symbols).

```julia
LmzA = wishartpol(1//2)          # Wishart with parameter 0.5
TLmz(LmzA)
```
"""
function TLmz(Lmz)
    L = numden(Lmz)
    im_ = _idx(m); iz = _idx(z)
    Dm = max(degree(L, im_), 0); Dz = max(degree(L, iz), 0)
    [coeff(L, [im_, iz], [i, j]) for i in 0:Dm, j in 0:Dz]
end
