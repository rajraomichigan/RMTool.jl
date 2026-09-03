"""
    RMTool

A Random Matrix and Free Probability Calculator for Julia.

Julia port of N. Raj Rao's RMTool MATLAB toolbox (v1.0, 2006). RMTool computes
the limiting spectral measure of a large class of random matrices using the
*polynomial method* (N. Raj Rao and Alan Edelman, "The polynomial method for
random matrices", arXiv:math.PR/0601389).

Probability measures are encoded by the bivariate polynomial `L(m, z)` satisfied
by their Stieltjes transform `m(z) = ∫ f(x)/(x - z) dx`. All symbolic
computations are performed with Nemo.jl (FLINT) instead of the MATLAB Symbolic
Math Toolbox.

Quick start:

```julia
using RMTool
b = wishartpol(1//2)
Lmz2MomS(b, 10)              # first 10 moments (Narayana polynomials at c = 1/2)
pdfinfo = Lmz2pdf(b, -0.05:0.01:5)
# using Plots; plot(pdfinfo.range, pdfinfo.density)
```
"""
module RMTool

# The symbolic ring below is created at load time and holds FLINT handles, so we
# do not let Julia serialise it into a precompile cache.
__precompile__(false)

using Nemo
using LinearAlgebra
using Random

export R, FF, RMV, m, z, c, rat, polyroots
export numden, subs, irreducLuv, Luv2Cu, L1plusL2, L1timesL2, TLmz
export Lmz2Lgz, Lgz2Lmz, Lgz2Lrg, Lrg2Lgz, Lmz2Lrg, Lrg2Lmz,
       Lmz2Lsy, Lsy2Lmz, Lmz2Lmyuz, Lmyuz2Lmz, Lmz2Letaz, Letaz2Lmz
export invA, shiftA, scaleA, mobiusA, transposeA, squareA, compressA,
       addAdtimes, AblockB, AplusB, AtimesB, AplusBkernel, AtimesBkernel
export wignerpol, wishartpol, atomLmz, equiLmz, corrWish, AtimesWish, AgramWish
export Lmz2pdf, Lmz2MomS, Lmz2MomF, PdfInfo
export wigner, wishart, haar, histw, theoryvssim, startRMTool

include("ring.jl")          # symbolic ring, substitution, helpers
# --- polynomial manipulation kernels
include("numden.jl")
include("irreducLuv.jl")
include("Luv2Cu.jl")
include("L1plusL2.jl")
include("L1timesL2.jl")
include("TLmz.jl")
# --- transform changes (Stieltjes m, Cauchy g, R-transform r, S-transform s, ...)
include("Lmz2Lgz.jl")
include("Lgz2Lmz.jl")
include("Lgz2Lrg.jl")
include("Lrg2Lgz.jl")
include("Lmz2Lrg.jl")
include("Lrg2Lmz.jl")
include("Lmz2Lsy.jl")
include("Lsy2Lmz.jl")
include("Lmz2Lmyuz.jl")
include("Lmyuz2Lmz.jl")
include("Lmz2Letaz.jl")
include("Letaz2Lmz.jl")
# --- deterministic and free-probabilistic operations on measures
include("mobiusA.jl")
include("invA.jl")
include("shiftA.jl")
include("scaleA.jl")
include("transposeA.jl")
include("squareA.jl")
include("compressA.jl")
include("addAdtimes.jl")
include("AblockB.jl")
include("AplusB.jl")
include("AtimesB.jl")
include("AplusBkernel.jl")
include("AtimesBkernel.jl")
# --- canonical ensembles
include("wignerpol.jl")
include("wishartpol.jl")
include("atomLmz.jl")
include("equiLmz.jl")
include("AtimesWish.jl")
include("AgramWish.jl")
include("corrWish.jl")
# --- extracting density and moments
include("Lmz2pdf.jl")
include("Lmz2MomS.jl")
include("Lmz2MomF.jl")
# --- numerical experiments
include("wigner.jl")
include("wishart.jl")
include("haar.jl")
include("histw.jl")
include("theoryvssim.jl")
include("startRMTool.jl")

end # module
