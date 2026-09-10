# RMTool.jl — a random matrix and free probability calculator (Julia)

Julia port of the [RMTool](http://www.mit.edu/~raj/rmtool) MATLAB toolbox by
N. Raj Rao (see also the [Python port](https://github.com/rajraomichigan/RMTool.py)).
RMTool computes the *limiting spectral measure* of a large class of random
matrix models **symbolically**, using the polynomial method of
> N. R. Rao and A. Edelman, *The polynomial method for random matrices*,
> Foundations of Computational Mathematics 8 (2008). [arXiv:math/0601389](https://arxiv.org/abs/math/0601389)

A probability measure is represented by the bivariate polynomial `L(m, z) = 0`
satisfied by its Stieltjes transform `m(z) = ∫ f(x)/(x − z) dx`. Operations on
random matrices (free sums, free products, Wishart-type transformations,
Möbius transforms, compressions, …) become operations on these polynomials,
carried out exactly with [Nemo.jl](https://github.com/Nemocas/Nemo.jl).
Densities and moments are then extracted from the resulting polynomial.

## Installation

```julia
julia> ]
pkg> add RMTool
```

(until registration completes, install directly from GitHub: `pkg> add https://github.com/rajraomichigan/RMTool`)

Requires Julia ≥ 1.10; Nemo.jl is installed as a dependency. Add `Plots` for figures:

```julia
pkg> add Plots
```

## Quick start

```julia
using RMTool, Plots

b = wishartpol(1//2)                         # Marcenko-Pastur, c = n/N = 1/2
Lmz2MomS(b, 10)                               # [1, 1, 3//2, 5//2, 37//8, ...]
info = Lmz2pdf(b, -0.05:0.01:5)
plot(info.range, info.density)

b = AplusB(wignerpol(), wishartpol(c))        # free sum, c kept symbolic
b                                              # c*m^3 + c*m^2*z + c*m + m^2 + m*z - m + 1
Lmz2MomS(b)                                    # [1, 1, c + 2, c^2 + 3c + 4, ...]

# theory vs. simulation
using LinearAlgebra
e = Float64[]
for _ in 1:1000
    A = wigner(100); B = wishart(100, 200)
    append!(e, real(eigvals(A + B)))
end
info = Lmz2pdf(AplusB(wignerpol(), wishartpol(1//2)))
yc, h = histw(e, 40)
bar(yc, h); plot!(info.range, info.density, color = :red, lw = 2)
```

Custom atomic measures are built with `atomLmz` (or `numden`):

```julia
L1 = atomLmz([1, 2], [1//2, 1//2])            # half the eigenvalues at 1, half at 2
L1 = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))   # same thing
```

## Function reference

Names follow the MATLAB toolbox (and the Python port) so existing scripts translate line-by-line.

| Function                                             | Operation                                                  |
| ----------------------------------------------------- | ------------------------------------------------------------ |
| `wignerpol()`                                         | semicircle law                                                |
| `wishartpol(c)`                                       | Marcenko-Pastur law, `c` = rows/columns (may be symbolic)      |
| `atomLmz(masses, weights)`                            | atomic measure                                                 |
| `equiLmz(t, M)`                                       | equilibrium measure for the potential `t·x^(2M)`                |
| `invA(L)`                                             | `inv(A)`                                                       |
| `shiftA(L, alpha)`                                    | `A + alpha·I`                                                  |
| `scaleA(L, alpha)`                                    | `alpha·A`                                                      |
| `mobiusA(L, p, q, r, s)`                              | `(p·A + q·I)/(r·A + s·I)`                                       |
| `transposeA(L, c)`                                    | `X'X` given `A = XX'`, `c` = size(A)/size(B)                    |
| `squareA(L)`                                          | `A²`                                                            |
| `AtimesWish(L, c)`                                    | `A × Wishart(c)` (Wishart with covariance A)                    |
| `AgramWish(L, c, s)`                                  | `(A_s + √s·G)(A_s + √s·G)'`                                     |
| `corrWish(La, Lb, c)`                                 | spatio-temporally correlated Wishart                            |
| `AplusB(La, Lb)`                                      | `A + QBQ'` (free additive convolution)                          |
| `AtimesB(La, Lb)`                                     | `A × QBQ'` (free multiplicative convolution)                    |
| `AblockB(La, Lb, c)`                                  | `diag(A, B)`                                                    |
| `compressA(L, c)`                                     | random compression of `A` by a factor `c < 1`                    |
| `addAdtimes(L, d)`                                    | measure with R-transform `d·R_A`                                 |
| `Lmz2MomS(L, n)` / `Lmz2MomF(L, n)`                   | first `n` moments                                                 |
| `Lmz2pdf(L, xx)`                                      | density along the grid `xx` (returns a `PdfInfo`)                  |
| `TLmz(L)`                                             | coefficient matrix of the bivariate polynomial                     |
| `Lmz2Lgz, Lgz2Lrg, Lmz2Lsy, Lmz2Lmyuz, Lmz2Letaz, …`  | transform conversions                                               |
| `irreducLuv, Luv2Cu, L1plusL2, L1timesL2, numden`     | low-level polynomial tools                                          |
| `wigner(n), wishart(n, m), haar(n)`                   | sample random matrices                                               |
| `histw(e, nbins), theoryvssim(...)`                   | normalised histogram / overlay data                                   |

## Notes for MATLAB and Python users

* The symbols `m, z, c` are exported directly (`using RMTool`), replacing
  `syms m z c`; the rest (`g, r, s, y, myu, eta, ...`) live in `RMV`, e.g. `RMV.r`.
* `Lmz2pdf` returns a `PdfInfo` struct with fields `range, rr, density, real,
  poles, multipleroots, mean`. For polynomials of degree > 3 in `m`, `density`
  has one column per root — isolate the right one manually, as in the MATLAB
  version.
* `Lmz2MomF` (which used Maple's `gfun`) simply calls `Lmz2MomS`, with a warning.
* Floating-point constants such as `0.5` are converted to exact rationals
  (`rationalize`) so that factorisation is reliable; prefer Julia's `1//2` for
  other values.
* `AplusBkernel` / `AtimesBkernel` are ported but, as in the original, experimental.
* `pretty`/`latex` → use Nemo's default printing, or `Latexify.jl` for LaTeX.

Full documentation (users guide + function reference) is at
**https://rajraomichigan.github.io/RMTool/**.

## License

GPL-2.0-or-later, as the original toolbox. If you use RMTool in a publication,
please cite the paper above and acknowledge the software.
