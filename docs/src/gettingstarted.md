# Getting started with RMTool

RMTool is a free, third-party toolbox for computing the limiting spectral
measure of a large class of random matrices. The techniques behind it are based
on the polynomial method (see arXiv: math.PR/0601389). Some applications are
mentioned in the afore-mentioned paper.

In the next two sections, we will provide a quick overview of how to get started.

## 1.1 System requirements and installation instructions

To install and run RMTool you need:

* Julia 1.10 or later
* [Nemo.jl](https://github.com/Nemocas/Nemo.jl) — the computer algebra package
  (wrapping FLINT) which plays the role of MATLAB's Symbolic Math Toolbox.
  It is installed automatically as a dependency.
* (optional) [Plots.jl](https://github.com/JuliaPlots/Plots.jl) to draw the
  densities and histograms.

RMTool can be run on any platform supported by Julia and Nemo (Linux, macOS,
Windows). RMTool is available for free under the GNU General Public License.

To install, extract the package directory (say `RMTool/`) anywhere you like and
type in the Julia package manager (press `]` at the REPL):

```
pkg> dev /path/to/RMTool
```

or, equivalently,

```julia
using Pkg
Pkg.develop(path = "/path/to/RMTool")
```

This resolves and installs Nemo.jl and completes the RMTool installation. The
package is loaded with `using RMTool`.

## 1.2 Other things you need to know

The directory in which you install RMTool contains subdirectories. Three of them are:

* `RMTool/docs`: containing this users manual (`make.jl` builds it with Documenter.jl),
* `RMTool/demos`: containing several demo files (`quickstart.jl`,
  `otherusefulcommands.jl`) implementing the examples in this manual,
* `RMTool/src`: one Julia file per function of the toolbox, with the same names
  as the original MATLAB `.m` files.

Throughout this users manual, we use the typewriter typeface to denote Julia
variables and functions, Julia commands that you should type, and results given
by Julia. Typing `?invA` (for example) at the REPL will provide you with some
concrete examples to experiment with.

## 1.3 Quick start

After installing the package, try the following sequence of commands in Julia.

```julia
julia> using RMTool
julia> startRMTool()
julia> b = wishartpol(1//2);
julia> wishart_moments = Lmz2MomS(b, 10)
julia> pdfinfo = Lmz2pdf(b, -0.05:0.01:5);
julia> using Plots; plot(pdfinfo.range, pdfinfo.density, lw = 2)
```

The above sequence of commands generates the bivariate polynomial
representation of the Marčenko–Pastur density and computes its first 10 moments,
and the density. We can repeat this computation for the Wigner matrix by typing
in the sequence of commands:

```julia
julia> b = wignerpol();
julia> wigner_moments = Lmz2MomS(b, 10)
julia> pdfinfo = Lmz2pdf(b, -4:0.01:4);
julia> plot(pdfinfo.range, pdfinfo.density, lw = 2)
```

Now try the following:

```julia
julia> b = AplusB(wignerpol(), wishartpol(c));   # Let c be symbolic
julia> moments = Lmz2MomS(b)
```

The function [`histw`](@ref) has been provided to compute the normalized
histogram of eigenvalues collected experimentally. Its use is illustrated below
with an example where the limiting spectrum is predicted and compared to an
experimental realization. (Type `?wishart` and `?wigner` for documentation of
these built-in functions.)

```julia
julia> using LinearAlgebra
julia> e = Float64[]; trials = 1000; n = 100;
julia> for idx in 1:trials
           A = wigner(n); B = wishart(n, 2n)
           append!(e, real(eigvals(A + B)))
       end
julia> pdfinfo = Lmz2pdf(AplusB(wignerpol(), wishartpol(1//2)));
julia> yc, h = histw(e, 40);
julia> bar(yc, h); plot!(pdfinfo.range, pdfinfo.density, color = :red, lw = 2)
```

If everything works, then the red line should coincide with the histogram bars.
Congratulations! You are now well on your way to using RMTool!

Note that in Julia the parameter `0.5` is written as the exact rational `1//2`:
RMTool computes with exact rational arithmetic (floating point inputs are
converted to the closest rational number).
