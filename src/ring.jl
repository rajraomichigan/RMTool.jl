# ---------------------------------------------------------------------------
# The symbolic ring.
#
# MATLAB's RMTool works with symbolic objects created by `syms m z c r g s y ...`.
# Here all polynomials live in one multivariate polynomial ring over QQ whose
# generators are the transform variables used internally by RMTool and a small
# set of free parameter symbols (c, alpha, beta, gamma, delta, a, b, p, q) that
# play the role of MATLAB's symbolic arguments.
# ---------------------------------------------------------------------------

const VARNAMES = ["m", "z", "g", "r", "s", "y", "myu", "eta", "f", "x", "t",
                  "c", "alpha", "beta", "gamma", "delta", "a", "b", "p", "q"]

const _RV = polynomial_ring(QQ, VARNAMES)

"The polynomial ring `QQ[m, z, g, r, s, y, myu, eta, f, x, t, c, alpha, ...]` in which all RMTool polynomials live."
const R = _RV[1]

"The fraction field of [`R`](@ref); use it to write rational expressions such as `m - rat(1//2)//(1 - z)`."
const FF = fraction_field(R)

const RElem = elem_type(R)
const FFElem = elem_type(FF)

"""
    RMV

NamedTuple of all symbolic variables of the ring [`R`](@ref):

* `m`, `z`   – Stieltjes transform variable and its argument (`m(z) = ∫ f(x)/(x-z) dx`)
* `g`        – Cauchy transform (`g = -m`)
* `r`        – R-transform variable
* `s`, `y`   – S-transform variables
* `myu`      – moment generating variable
* `eta`      – eta-transform variable
* `f`, `x`, `t` – auxiliary variables used by the kernel routines and resultants
* `c`, `alpha`, `beta`, `gamma`, `delta`, `a`, `b`, `p`, `q` – free parameters
  (the analogue of MATLAB's `syms c`)

`m`, `z` and `c` are also exported directly, mirroring `syms m z c` in MATLAB.
"""
const RMV = NamedTuple{Tuple(Symbol.(VARNAMES))}(Tuple(_RV[2]))

const m = RMV.m
const z = RMV.z
const c = RMV.c

_idx(v::RElem) = var_index(v)

# ---- conversion of numbers / ring elements into the ring and its fraction field
_qq(x::Integer)       = QQ(x)
_qq(x::Rational)      = QQFieldElem(ZZ(numerator(x)), ZZ(denominator(x)))
_qq(x::AbstractFloat) = _qq(rationalize(BigInt, x))
_qq(x::QQFieldElem)   = x

"""
    rat(x)

Embed a Julia number (`Integer`, `Rational`, `AbstractFloat`) into the ring `R`.
Floats are converted to the closest rational number (`rationalize`), so prefer
exact rationals such as `1//2` where possible.
"""
rat(x::Number) = R(_qq(x))
rat(x::RElem)  = x
rat(x::Union{ZZRingElem,QQFieldElem}) = R(x)
rat(x::FFElem) = (isone(denominator(x)) ? numerator(x) : error("rat: $x is not a polynomial"))

_ff(x::FFElem) = x
_ff(x::RElem)  = FF(x)
_ff(x::Number) = FF(rat(x))
_ff(x::Union{ZZRingElem,QQFieldElem}) = FF(rat(x))

_ffvec() = FFElem[FF(v) for v in _RV[2]]

# convert a constant polynomial / rational to Float64 or Rational{BigInt}
function _rational(q::QQFieldElem)
    BigInt(numerator(q)) // BigInt(denominator(q))
end
function _rational(p::RElem)
    is_constant(p) || error("expected a constant polynomial, got $p")
    _rational(constant_coefficient(p))
end
function _rational(x::FFElem)
    _rational(numerator(x)) / _rational(denominator(x))
end
_f64(x) = Float64(_rational(x))

"""
    subs(L, var => value, ...)

Simultaneous substitution of variables of [`R`](@ref) by arbitrary rational
expressions (elements of `R`, of `FF`, or Julia numbers). The result is the
*numerator* of the substituted expression, i.e. the denominator is cleared
exactly as MATLAB's `numden(subs(...))` does.

```julia
subs(wignerpol(), m => -RMV.g)          # Cauchy transform polynomial
subs(wishartpol(c), c => 1//2)          # fix the parameter
```
"""
function subs(L::RElem, prs::Pair...)
    vals = _ffvec()
    for (v, val) in prs
        vals[_idx(v)] = _ff(val)
    end
    numden(evaluate(L, vals))
end
function subs(L::FFElem, prs::Pair...)
    vals = _ffvec()
    for (v, val) in prs
        vals[_idx(v)] = _ff(val)
    end
    numden(evaluate(numerator(L), vals) / evaluate(denominator(L), vals))
end
subs(L::Number, prs::Pair...) = rat(L)

# univariate coefficient list of P in variable index i (entries are polynomials in the other variables)
function _ucoeffs(P::RElem, i::Int)
    D = degree(P, i)
    D < 0 && return RElem[]
    RElem[coeff(P, [i], [k]) for k in 0:D]
end

# resultant with respect to the variable of index i (Nemo/FLINT, Sylvester-matrix fallback)
function _resultant(A::RElem, B::RElem, i::Int)
    try
        return resultant(A, B, i)
    catch err
        err isa MethodError || rethrow()
    end
    ca = _ucoeffs(A, i); cb = _ucoeffs(B, i)
    dA = length(ca) - 1; dB = length(cb) - 1
    n = dA + dB
    n == 0 && return one(R)
    S = zero_matrix(R, n, n)
    for row in 1:dB, k in 0:dA
        S[row, row + k] = ca[dA - k + 1]
    end
    for row in 1:dA, k in 0:dB
        S[dB + row, row + k] = cb[dB - k + 1]
    end
    det(S)
end

# discriminant with respect to the variable of index i (up to a nonzero constant factor)
function _discriminant(A::RElem, i::Int)
    try
        return discriminant(A, i)
    catch err
        err isa MethodError || rethrow()
    end
    _resultant(A, derivative(A, i), i)
end

# indices of free parameter symbols present in a polynomial
const _PARAM_IDX = Set(_idx(RMV[k]) for k in (:c, :alpha, :beta, :gamma, :delta, :a, :b, :p, :q))
_params_in(L::RElem) = [v for v in vars(L) if _idx(v) in _PARAM_IDX]

# ---- numerics: roots of a real/complex univariate polynomial via the companion matrix
"""
    polyroots(coeffs)

Roots of the polynomial `coeffs[1] + coeffs[2] x + ... + coeffs[end] x^(n-1)`
(ascending order of powers, the opposite of MATLAB's `roots`), computed as the
eigenvalues of the companion matrix.
"""
function polyroots(cs::AbstractVector{<:Number})
    n = findlast(!iszero, cs)
    (n === nothing || n == 1) && return ComplexF64[]
    cs = ComplexF64.(cs[1:n]) ./ cs[n]
    D = n - 1
    C = zeros(ComplexF64, D, D)
    for i in 2:D
        C[i, i-1] = 1
    end
    C[:, D] .= -cs[1:D]
    eigvals(C)
end

_realroots(rs; tol = 1e-8) = sort(real.(filter(r -> abs(imag(r)) < tol * max(1, abs(r)), rs)))

# Float64 coefficient list of a polynomial in one variable of index i (all other variables must be absent)
function _ucoeffs_f64(P::RElem, i::Int)
    cs = _ucoeffs(P, i)
    map(cs) do ck
        is_constant(ck) || error("polynomial still depends on $(vars(ck)); substitute values for these symbols first")
        _f64(ck)
    end
end

# division in the fraction field, accepting ring elements and numbers on either side
_div(a, b) = _ff(a) / _ff(b)
_inv(x) = _div(1, x)
