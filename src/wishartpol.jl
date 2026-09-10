"""
    wishartpol(c)

Bivariate polynomial `m*(1 - c - c*m*z - z) - 1` encoding the Marčenko–Pastur
law with parameter `c = rows/columns` (the limiting spectral measure of
[`wishart`](@ref) matrices). `c` can be a number (`1//2`) or the symbolic
parameter [`c`](@ref RMV).

```julia
b = wishartpol(1//2)
b = wishartpol(c)      # keep c symbolic
```
"""
function wishartpol(cc = c)
    cc = rat(cc)
    m * (1 - cc - cc * m * z - z) - 1
end
