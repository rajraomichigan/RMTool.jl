# Section 2.5 of the users guide: some other useful commands.
using RMTool
using Nemo

b1 = wishartpol(c)
b2 = wignerpol()
b3 = AtimesB(b1, b2)
println(b3)
moments = Lmz2MomS(b3, 6)
println(moments)
println(TLmz(b3))
# LaTeX output: Nemo prints polynomials as text; for LaTeX use Latexify.jl or
# string(b3) and edit by hand.
