# Quick start (Section 1.3 of the RMTool users guide), Julia version.
using RMTool
using LinearAlgebra
# using Plots     # uncomment to draw the figures

startRMTool()

b = wishartpol(1//2)
wishart_moments = Lmz2MomS(b, 10)
pdfinfo = Lmz2pdf(b, -0.05:0.01:5)
# plot(pdfinfo.range, pdfinfo.density, lw = 2)

b = wignerpol()
wigner_moments = Lmz2MomS(b, 10)
pdfinfo = Lmz2pdf(b, -4:0.01:4)
# plot(pdfinfo.range, pdfinfo.density, lw = 2)

# Let c be symbolic
b = AplusB(wignerpol(), wishartpol(c))
moments = Lmz2MomS(b)

# Theory versus simulation
e = Float64[]; trials = 1000; n = 100
for idx in 1:trials
    A = wigner(n); B = wishart(n, 2n)
    append!(e, real(eigvals(A + B)))
end
pdfinfo = Lmz2pdf(AplusB(wignerpol(), wishartpol(1//2)))
d = theoryvssim(e, 40, pdfinfo.density, pdfinfo.range)
# bar(d.centers, d.heights, label = "simulation")
# plot!(d.range, d.density, color = :red, lw = 2, label = "theory")
# If everything works, the red line should coincide with the histogram bars.
