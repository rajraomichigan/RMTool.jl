# Computing the spectral measure of random matrices

RMTool can compute the spectral measure of random matrix *models*. In other
words, we do not actually generate the matrix. All the computations are done
symbolically. Manipulating polynomials is at the heart of the software. Thus we
have to learn how to represent probability measures as polynomials.

## 2.1 Polynomial representation of probability measures

Polynomials in RMTool are elements of the multivariate polynomial ring
[`R`](@ref) over the rationals provided by Nemo.jl. Typically, a polynomial is
created by using the symbolic variables exported by RMTool and then constructing
it using the algebraic manipulations.

Generically we represent a probability measure by the bivariate polynomial that
encodes its Stieltjes transform. The Stieltjes transform is defined as:

```math
m(z) = \int \frac{1}{x - z} f(x)\,dx \qquad \text{for } \Im z \neq 0.
```

Hence we represent probability measures by the variables `m` and `z`. In Julia
the symbolic variables `m`, `z` (and the parameter `c`) are brought into scope
by `using RMTool` (the analogue of `syms m z c` in MATLAB); all other variables
are fields of [`RMV`](@ref), e.g. `RMV.r`.

Then we construct the polynomial `L(m, z)` as follows:

```julia
julia> L1 = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
```

The above polynomial encodes the probability measure with atoms at 1 and 2 of
equal weight. The command [`numden`](@ref) clears the denominator. (The
function [`atomLmz`](@ref) builds the same polynomial directly:
`atomLmz([1, 2], [1//2, 1//2])`.) Rational numbers must be written as Julia
rationals (`1//2`) or wrapped with [`rat`](@ref); the operator `//` between two
ring elements creates a rational expression whose denominator `numden` clears.

## 2.2 Manipulating the polynomials

Polynomials such as the one created above can be manipulated with the usual
operators: `+`, `-`, `*`. However, since we are interested in manipulating the
probability measures that they encode we do not need to worry about the usual
algebraic operations involving the polynomials.

Instead, we use built-in functions to manipulate these polynomials. In doing so,
we are really manipulating the probability measures encoded by these
polynomials. The functions listed below can be used to this effect. The names
are self-explanatory (see the paper for additional details). The polynomials can
be manipulated using the function (with the arguments indicated). The arguments
can themselves be symbolic (e.g. the parameter `c`).

Typing in `?invA` (for example) will provide you with some concrete examples
to experiment with.

| Function                  | Operation                                                    |
|:--------------------------|:-------------------------------------------------------------|
| `invA(L)`                 | Invert A                                                     |
| `shiftA(L, alpha)`        | A + alpha * I                                                |
| `scaleA(L, alpha)`        | alpha * A                                                    |
| `mobiusA(L, p, q, r, s)`  | (p*A + q*I)/(r*A + s*I)                                      |
| `transposeA(L, c)`        | If A = XX' then B = X'X, c = Size of A/Size of B             |
| `squareA(L)`              | A²                                                           |
| `AtimesWish(L, c)`        | A × Wishart with parameter c = Rows/Columns                  |
| `AgramWish(L, c, s)`      | (X + s G)(X + s G)' with c = Rows/Columns                    |
| `AplusB(La, Lb)`          | A + B (free addition)                                        |
| `AtimesB(La, Lb)`         | A × B (free multiplication)                                  |
| `compressA(La, c)`        | compress A by a factor of c (less than 1)                    |
| `AblockB(La, Lb, c)`      | diag(A, B) with c = size(A)/size(diag(A, B))                 |
| `addAdtimes(L, d)`        | free additive convolution of d copies of A                   |
| `corrWish(La, Lb, c)`     | spatio-temporally correlated Wishart Y Y', Y = A^{1/2} G B^{1/2} |
| `equiLmz(t, M)`           | equilibrium measure for the potential V(x) = t x^{2M}        |

The built-in ensembles are `wignerpol()` (semicircle law) and `wishartpol(c)`
(Marčenko–Pastur law).

## 2.3 Extracting the density from the polynomial

Extracting the density from the polynomial in `m` and `z` involves determining
the "right root". This is, in general, difficult so unless the right root can be
easily spotted, we just return all the real and imaginary roots.

From a bivariate polynomial we can return the roots and some additional
information about the density by typing the command:

```julia
julia> xx = xstart:xstep:xend;
julia> pdfinfo = Lmz2pdf(L, xx);
```

The variable `pdfinfo` returned is a [`PdfInfo`](@ref) structure containing
many fields. Type `?Lmz2pdf` for more details. We plot the (normalized)
imaginary roots by typing the command:

```julia
julia> plot(pdfinfo.range, pdfinfo.density, seriestype = :scatter)
```

Somewhere among all the roots is the "right root". Often it is easy enough to
spot it. The other fields (`poles`, `multipleroots`, `mean`) will help you
determine the region of support as described in the accompanying paper.

## 2.4 Extracting the moments from the polynomial

Even when the density cannot be extracted, the moments can often be. To do so
there are two commands available [`Lmz2MomS`](@ref) and [`Lmz2MomF`](@ref),
where the suffix differentiates between the slow and fast implementations.

### 2.4.1 Slow implementation

As you have seen before you can type in the command:

```julia
julia> Moments = Lmz2MomS(L, number_of_moments)
```

to enumerate the moments. Here `L` is the polynomial in `m` and `z`. The
moments are returned as exact rational numbers, or as rational functions of the
symbolic parameters when the polynomial contains some (e.g. `wishartpol(c)`
gives the Narayana polynomials in `c`).

### 2.4.2 Fast implementation

The fast implementation of the MATLAB toolbox relied on Maple's `gfun` package
(`algeqtoseries`) shipped with the Extended Symbolic Toolbox. No such package is
available in Julia, so

```julia
julia> Moments = Lmz2MomF(L, number_of_moments)
```

issues a warning and calls the slow implementation, exactly as the MuPAD-era
MATLAB version did. In practice the exact expansion used by `Lmz2MomS` (each
moment is obtained by solving one linear equation) is fast enough for the
examples in this manual.

## 2.5 Some other useful commands

One of the great benefits of using symbolic software is being able to reproduce
the result as often as one likes without errors. Of course, one might want to
use the results generated using this package and put it in a technical
publication (please acknowledge use of this software). There are a couple of
built-in functions that are useful in this regard.

Nemo prints polynomials in a legible form; `string(L)` gives the same as text,
and the command [`TLmz`](@ref) rewrites the bivariate polynomial as a matrix of
coefficients. Try the following sequence of commands in Julia:

```julia
julia> b1 = wishartpol(c)
julia> b2 = wignerpol();
julia> b3 = AtimesB(b1, b2)
julia> moments = Lmz2MomS(b3, 6)
julia> TLmz(b3)
```

For LaTeX output, packages such as Latexify.jl can be applied to `string(b3)`.
