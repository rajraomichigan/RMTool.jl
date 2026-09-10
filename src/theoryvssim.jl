"""
    theoryvssim(e, numberofbins, imagpdf, pdfrange) -> NamedTuple

Compares the empirically collected eigenvalues `e` with the theoretical
density `imagpdf` evaluated on `pdfrange`. Returns
`(centers, heights, range, density)` where `(centers, heights)` is the
normalized histogram from [`histw`](@ref).

The MATLAB version draws the figure directly; here plotting is left to the
user so that RMTool does not depend on a plotting package:

```julia
using Plots
d = theoryvssim(e, 40, pdfinfo.density, pdfinfo.range)
bar(d.centers, d.heights; label = "simulation")
plot!(d.range, d.density; color = :red, lw = 2, label = "theory", xlabel = "x", ylabel = "Probability")
```
"""
function theoryvssim(e::AbstractVector{<:Real}, numberofbins::Integer, imagpdf, pdfrange)
    yc, h = histw(e, numberofbins)
    (centers = yc, heights = h, range = collect(Float64, pdfrange), density = imagpdf)
end
