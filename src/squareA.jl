"""
    squareA(LmzA)

Bivariate polynomial for `B = A^2`.

Implementation: with `w^2 = z`, the Stieltjes transform of `A^2` is
`m_B(z) = (m_A(w) - m_A(-w))/(2w)`, so the polynomials `L_A(2 m w, w)` and
`L_A(-2 m w, -w)` are combined with [`L1plusL2`](@ref) and `w^2` is replaced by `z`.
(The MATLAB version uses `sqrt(z)` symbolically.)

```julia
LmzA = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
LmzB = squareA(LmzA)
MomA = Lmz2MomS(LmzA, 10)
MomB = Lmz2MomS(LmzB, 5)      # compare the moments of A and B = A^2
```
"""
function squareA(LmzA)
    L = numden(LmzA)
    w = RMV.t                     # auxiliary variable playing the role of sqrt(z)
    iw = _idx(w)
    Lmz1 = subs(L, z => w,  m => 2 * m * w)
    Lmz2 = subs(L, z => -w, m => -2 * m * w)
    L3 = L1plusL2(Lmz1, Lmz2, m)       # polynomial in m and w
    # L3 is even (or odd) in w; strip a stray factor of w and replace w^2 -> z
    while degree(L3, iw) > 0 && iszero(coeff(L3, [iw], [0]))
        L3 = divexact(L3, w)
    end
    if all(iseven(e) for e in _exps_of(L3, iw))
        LmzB = sum(coeff(L3, [iw], [e]) * z^(e ÷ 2) for e in _exps_of(L3, iw); init = zero(R))
    else
        LmzB = _resultant(L3, w^2 - z, iw)
    end
    irreducLuv(LmzB, m, z)
end

# set of exponents of the variable of index i occurring in P
_exps_of(P::RElem, i::Int) = Set(k for k in 0:degree(P, i) if !iszero(coeff(P, [i], [k])))
