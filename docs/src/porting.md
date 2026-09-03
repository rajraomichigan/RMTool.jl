# Porting notes: MATLAB → Julia

The Julia package keeps the structure and the names of the MATLAB toolbox: every
`.m` file became a `.jl` file of the same name in `src/` (the outdated duplicate
`irreducLuv.txt` was not ported). The differences worth knowing are:

* **Symbolic engine.** The MATLAB Symbolic Math Toolbox (MuPAD/Maple) is
  replaced by [Nemo.jl](https://nemocas.github.io/Nemo.jl/stable/). All
  polynomials live in the single ring [`R`](@ref) `= QQ[m, z, g, r, s, y, myu,
  eta, f, x, t, c, alpha, beta, gamma, delta, a, b, p, q]`; `m`, `z`, `c` are
  exported, the others are available as `RMV.g`, `RMV.alpha`, ... The symbols
  `c, alpha, beta, gamma, delta, a, b, p, q` are free parameters and replace
  MATLAB's `syms c`.
* **Exact arithmetic.** Numbers are exact rationals: write `1//2` instead of
  `0.5`. Floating point values are accepted and converted to the nearest
  rational with `rationalize`.
* **Rational expressions.** `numden(expr)` returns the numerator of an element
  of the fraction field [`FF`](@ref); `a // b` between two ring elements creates
  such an element, and `rat(1//2)` embeds a number into the ring.
* **`subs`.** [`subs`](@ref)`(L, var => value, ...)` performs a *simultaneous*
  substitution and clears the denominator (`numden(subs(...))` in MATLAB).
  Where the MATLAB code substitutes sequentially, the Julia code does too.
* **`irreducLuv`.** MATLAB's `factor` returned the factors in an order and the
  toolbox kept the last one. The Julia version drops constant and monomial
  factors and keeps the irreducible factor of highest total degree, with a
  warning when the choice was ambiguous. `irreducLuv(L, u, v, Val(:minimal))`
  returns `(L, minimal)` like the two-output MATLAB call. `AtimesWish` and
  `AgramWish` explicitly divide out the spurious factor introduced by their
  substitution before calling `irreducLuv`.
* **Resultants.** `L1plusL2` and `L1timesL2` use Nemo/FLINT resultants (with a
  Sylvester-matrix fallback); `method = :matrix` selects the Kronecker-product
  companion-matrix implementation described in the paper.
* **`squareA`.** `sqrt(z)` is not a ring element; an auxiliary variable `w`
  with `w^2 = z` is used and eliminated at the end.
* **`Lmz2MomS`.** The moment recursion is implemented as the exact power-series
  expansion of the moment generating function about `z = 0` (one linear solve
  per moment). Several expansions are enumerated only when the polynomial
  factors. `Lmz2MomF` (Maple `gfun`) falls back to `Lmz2MomS`.
* **`Lmz2pdf`.** Returns a [`PdfInfo`](@ref) struct; roots in `m` are computed
  numerically with a companion matrix ([`polyroots`](@ref)). Parameters left in
  the polynomial can be fixed with the keyword `params = Dict(c => 1//3)`
  (default `1//2` with a warning, as in MATLAB).
* **Plotting.** RMTool does not depend on a plotting package. `histw` returns
  `(centers, heights)` and `theoryvssim` returns the data needed for the figure;
  draw them with Plots.jl (`bar`, `plot!`).
* **`TLmz`, `latex`, `pretty`.** `TLmz` returns a `Matrix` of ring elements.
  Nemo's printing replaces `pretty`; LaTeX can be produced with Latexify.jl.
* **`startRMTool`.** The MATLAB start-up script is now a function that prints
  the available variables and returns `wishartpol(c)`.

## File map

| MATLAB file        | Julia file (in `src/`) | Notes |
|:-------------------|:-----------------------|:------|
| `startRMTool.m`    | `startRMTool.jl`       | function instead of script |
| `wignerpol.m`, `wishartpol.m`, `atomLmz.m`, `equiLmz.m` | same names | `equiLmz` rationalizes the constant `a²` |
| `wigner.m`, `wishart.m`, `haar.m` | same names | `LinearAlgebra`/`Random` |
| `histw.m`, `theoryvssim.m` | same names | return data instead of plotting |
| `Lmz2Lgz.m`, `Lgz2Lmz.m`, `Lgz2Lrg.m`, `Lrg2Lgz.m`, `Lmz2Lrg.m`, `Lrg2Lmz.m`, `Lmz2Lsy.m`, `Lsy2Lmz.m`, `Lmz2Lmyuz.m`, `Lmyuz2Lmz.m`, `Lmz2Letaz.m`, `Letaz2Lmz.m` | same names | |
| `irreducLuv.m` (`.txt`) | `irreducLuv.jl` | see above |
| `Luv2Cu.m`, `L1plusL2.m`, `L1timesL2.m`, `TLmz.m` | same names | |
| `mobiusA.m`, `invA.m`, `shiftA.m`, `scaleA.m`, `transposeA.m`, `squareA.m`, `compressA.m`, `addAdtimes.m`, `AblockB.m`, `AplusB.m`, `AtimesB.m`, `AplusBkernel.m`, `AtimesBkernel.m` | same names | |
| `AtimesWish.m`, `AgramWish.m`, `corrWish.m` | same names | |
| `Lmz2pdf.m`, `Lmz2MomS.m`, `Lmz2MomF.m` | same names | |
| —                  | `ring.jl`, `numden.jl` | new: ring definition, `subs`, `numden`, `rat`, `polyroots` |
