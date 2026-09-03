"""
    irreducLuv(Luv, u, v; warn = true) -> Luv
    irreducLuv(Luv, u, v, Val(:minimal)) -> (Luv, minimal)

Makes a bivariate polynomial `Luv` irreducible with respect to `u` and `v`:
the denominator is cleared, the polynomial is factored (Nemo `factor`),
constant and monomial factors are dropped and a single irreducible factor is
returned.

If more than one non-trivial factor remains, the factor of highest total degree
is returned (the MATLAB version returns the last factor found by `factor`), a
warning is printed (unless `warn = false`) and `minimal == false`.

Disclaimer (from the original): this is NOT always the minimal representation.
When the polynomial encodes a compactly supported density, the correct factor
MUST admit a moment expansion about zero; running [`Lmz2MomS`](@ref) on each
factor of `factor(Luv)` identifies the "right factor".

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = irreducLuv(LmzA, m, z)
```
"""
function irreducLuv(Luv, u, v; warn::Bool = true)
    first(irreducLuv(Luv, u, v, Val(:minimal); warn = warn))
end

function irreducLuv(Luv, u, v, ::Val{:minimal}; warn::Bool = true)
    L = numden(Luv)
    iszero(L) && return (L, true)
    fa = factor(L)
    facs = RElem[]
    for (pf, e) in fa
        (is_constant(pf) || length(pf) == 1) && continue   # drop constants and pure monomials
        push!(facs, pf)
    end
    isempty(facs) && return (L, true)
    length(facs) == 1 && return (facs[1], true)
    sort!(facs, by = total_degree)
    warn && @warn "Polynomial can be factored even more; might be useful to identify the correct factor. " *
                  "Type factor(L) to see the factorization. Returning the factor of highest degree." facs
    return (facs[end], false)
end
