# Notes for MATLAB and Python users

RMTool.jl keeps the names of the original MATLAB toolbox — and of the
[Python port](https://github.com/rajraomichigan/RMTool.py) — so scripts
translate line-by-line. The differences worth knowing:

## Symbolic engine

MATLAB's Symbolic Math Toolbox and Python's SymPy are both replaced by
[Nemo.jl](https://github.com/Nemocas/Nemo.jl) (FLINT). All polynomials live in
one ring `R = QQ[m, z, g, r, s, y, myu, eta, f, x, t, c, alpha, beta, gamma,
delta, a, b, p, q]`. `m`, `z`, `c` are exported directly (`using RMTool`,
mirroring `syms m z c` / `from rmtool import *`); the remaining variables are
fields of `RMV`, e.g. `RMV.g`, `RMV.myu`.

## Exact arithmetic

Numbers are exact rationals. Write `1//2` (Julia) where MATLAB used `0.5` and
Python used `sp.Rational(1, 2)`. Floating point values are still accepted and
converted to the nearest rational with `rationalize`, exactly as the Python
port converts floats to `sympy.Rational` for reliable factorisation.

## Rational expressions and substitution

* `numden(expr)` clears the denominator, as in both MATLAB and Python.
* `subs(L, var => value, ...)` performs a simultaneous substitution and
  returns `numden(subs(...))`, matching the MATLAB convention. `a // b`
  between two ring elements creates a fraction; `rat(1//2)` embeds a number
  into the ring `R` (the Julia analogue of `sp.Rational(1,2)`).

## `irreducLuv`

Drops constant and monomial factors and keeps the irreducible factor of
highest total degree, printing a warning (`irreducLuv(L, u, v; warn=false)` to
silence it) when the choice was ambiguous — the same disclaimer as the
original: run `Lmz2MomS` on each factor of `factor(L)` to identify the "right"
one when a measure is compactly supported.

## `Lmz2MomS` / `Lmz2MomF`

`Lmz2MomS` expands the moment generating polynomial exactly, one linear solve
per moment, and returns a `Vector` of exact rationals (or of rational functions
of the free parameters, e.g. of `c`, when the input polynomial is symbolic).
`Lmz2MomF`, which relied on Maple's `gfun` package in MATLAB, simply calls
`Lmz2MomS` with a warning — exactly as the Python port and the MuPAD-era
MATLAB code already did.

## `Lmz2pdf`

Returns a `PdfInfo` struct with fields `range, rr, density, real, poles,
multipleroots, mean` (the Julia analogue of Python's `PdfInfo` dataclass).
For polynomials of degree > 3 in `m`, `density` has one column per root; no
automatic heuristic picks the "right" one — isolate it manually as in the
MATLAB version. Any symbolic parameter left in the polynomial (e.g. `c`) can be
fixed with `params = Dict(c => 1//3)`, defaulting to `1/2` with a warning.

## Plotting

Neither `histw` nor `theoryvssim` draws a figure — `histw` returns
`(centers, heights)` and `theoryvssim` returns the data needed for the
overlay. Draw them with [Plots.jl](https://github.com/JuliaPlots/Plots.jl)
(`bar`, `plot!`), the same separation of concerns as the Python port, which
returns data for Matplotlib rather than calling `plt.plot` itself.

## Output formatting

`pretty` (MATLAB) / `sympy.pretty`, `sympy.latex` (Python) → Nemo's default
`show` methods print polynomials as text; for LaTeX, apply
[Latexify.jl](https://github.com/korsbo/Latexify.jl) to `string(L)`.

## Status of the kernel routines

`AplusBkernel` and `AtimesBkernel` are ported but, as in the original and in
the Python port, experimental.

## File map

| MATLAB file (`.m`)  | Julia (`src/`)        | Python (`rmtool/`) |
|:---------------------|:-----------------------|:---------------------|
| `startRMTool.m`      | `startRMTool.jl`        | — |
| `wignerpol.m`, `wishartpol.m`, `atomLmz.m`, `equiLmz.m` | same names | same names |
| `wigner.m`, `wishart.m`, `haar.m` | same names | same names |
| `histw.m`, `theoryvssim.m` | same names, return data | same names, return data |
| transform conversions (`Lmz2Lgz.m`, …) | same names | same names |
| `irreducLuv.m`, `Luv2Cu.m`, `L1plusL2.m`, `L1timesL2.m`, `TLmz.m` | same names | same names |
| `mobiusA.m`, `invA.m`, `shiftA.m`, `scaleA.m`, `transposeA.m`, `squareA.m`, `compressA.m`, `addAdtimes.m`, `AblockB.m`, `AplusB.m`, `AtimesB.m`, `AplusBkernel.m`, `AtimesBkernel.m` | same names | same names |
| `AtimesWish.m`, `AgramWish.m`, `corrWish.m` | same names | same names |
| `Lmz2pdf.m`, `Lmz2MomS.m`, `Lmz2MomF.m` | same names | same names |
| —                    | `ring.jl`, `numden.jl` (new: ring, `subs`, `numden`, `rat`) | `numden`, module-level ring setup |
