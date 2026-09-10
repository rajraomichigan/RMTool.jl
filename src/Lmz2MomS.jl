"""
    Lmz2MomS(Lmz, number_of_moments = 4; numeric = true, warn = true)

Slow (exact, symbolic) algorithm for enumerating the moments of the measure
encoded by `Lmz`. Returns a vector `[M_0 = 1, M_1, ..., M_n]`.

The moment generating polynomial `L(myu, z)` ([`Lmz2Lmyuz`](@ref)) is expanded
about `z = 0` with `myu(0) = 1`, one coefficient at a time. Moments are
returned as exact rationals (`Rational{BigInt}`) when they are numbers and as
elements of `FF` (rational functions of the parameters, e.g. of `c`) otherwise;
pass `numeric = false` to always get `FF` elements.

If more than one seemingly valid moment sequence exists (several factors of
`Lmz` admit an expansion about zero), all are enumerated as the columns of a
matrix with a warning; positivity arguments can be used to isolate the correct
sequence.

```julia
LmzA = wishartpol(c)            # Wishart with symbolic parameter c
Moments = Lmz2MomS(LmzA, 6)     # Narayana polynomials in c
Lmz2MomS(wignerpol(), 8)        # Catalan numbers 1, 0, 1, 0, 2, 0, 5, 0, 14
```

Reference: N. Raj Rao and Alan Edelman, "The polynomial method for random matrices".
"""
function Lmz2MomS(Lmz, max_moment::Integer = 4; numeric::Bool = true, warn::Bool = true)
    myu = RMV.myu
    imyu = _idx(myu); iz = _idx(z)
    Lmyuz = Lmz2Lmyuz(numden(Lmz))

    # candidate polynomials: the whole polynomial, then (if needed) its factors
    candidates = RElem[Lmyuz]
    fa = factor(Lmyuz)
    for (pf, e) in fa
        (is_constant(pf) || pf == Lmyuz) && continue
        push!(candidates, pf)
    end

    sequences = Vector{Vector{FFElem}}()
    for (ci, P) in enumerate(candidates)
        seq = _moment_expansion(P, imyu, iz, max_moment)
        if seq !== nothing
            push!(sequences, seq)
            ci == 1 && break          # the full polynomial gave a unique expansion
        end
    end

    if isempty(sequences)
        warn && @warn "No moments exist. The bivariate polynomial either does not encode a valid probability measure or the measure is non-compact."
        return vcat([1.0], fill(NaN, max_moment))
    end

    if length(sequences) > 1
        warn && @warn "There are $(length(sequences)) expansions about zero with initial coefficient equal to 1. All have been returned. Eliminate some or all of them by positivity arguments."
        M = reduce(hcat, sequences)
        return numeric ? _try_numeric(M) : M
    end
    seq = sequences[1]
    return numeric ? _try_numeric(seq) : seq
end

# Expand myu(z) = 1 + M_1 z + M_2 z^2 + ... from P(myu, z) = 0 by matching the
# coefficient of z^k at every order; returns nothing when P(1, 0) != 0 or when
# the expansion is not unique (∂P/∂myu(1,0) = 0).
function _moment_expansion(P::RElem, imyu::Int, iz::Int, N::Integer)
    myu = RMV.myu
    P10 = subs(P, myu => 1, z => 0)
    iszero(P10) || return nothing
    d = subs(derivative(P, imyu), myu => 1, z => 0)
    iszero(d) && return nothing
    dF = FF(d)
    moms = FFElem[_ff(1)]
    zf = FF(z)
    for k in 1:N
        S = sum(moms[j+1] * zf^j for j in 0:k-1)
        vals = _ffvec(); vals[imyu] = S
        E = evaluate(P, vals)
        nu = numerator(E); de = denominator(E)
        degree(de, iz) == 0 || error("Lmz2MomS: unexpected z-dependence in the denominator")
        ck = coeff(nu, [iz], [k])
        push!(moms, -(FF(ck) / FF(de)) / dF)
    end
    moms
end

function _try_numeric(v::AbstractVector{FFElem})
    all(x -> is_constant(numerator(x)) && is_constant(denominator(x)), v) || return v
    [_rational(x) for x in v]
end
function _try_numeric(M::AbstractMatrix{FFElem})
    all(x -> is_constant(numerator(x)) && is_constant(denominator(x)), M) || return M
    [_rational(x) for x in M]
end
