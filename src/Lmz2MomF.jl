"""
    Lmz2MomF(Lmz, number_of_moments = 4; kwargs...)

Fast algorithm for enumerating moments from `Lmz`.

In MATLAB this required the Extended Symbolic Toolbox (Maple's
`gfun[algeqtoseries]`); that package is not available in Julia, so — exactly
as the MuPAD-era MATLAB version did — this function falls back to
[`Lmz2MomS`](@ref) with a warning.

```julia
LmzA = wishartpol(c)
Moments = Lmz2MomF(LmzA, 6)
```
"""
function Lmz2MomF(Lmz, max_moment::Integer = 4; kwargs...)
    @warn "algeqtoseries (Maple gfun) is not available - using Lmz2MomS instead"
    Lmz2MomS(Lmz, max_moment; kwargs...)
end
