"""
    histw(y, M) -> (yc, h)

Weighted (normalized) histogram of the data in `y` using `M` boxes. Returns
the bin centers `yc` and the normalized heights `h` (so that the histogram
integrates to one). Plot it with e.g. `bar(yc, h)` from Plots.jl.

```julia
y = randn(100_000)
yc, h = histw(y, 50)
```

Per-Olof Persson, September 2004 (MATLAB original).
"""
function histw(y::AbstractVector{<:Real}, M::Integer)
    ymin, ymax = extrema(y)
    dy = (ymax - ymin) / M
    yy = ymin .+ dy .* (0:M)
    yc = (yy[1:end-1] .+ yy[2:end]) ./ 2
    h = zeros(Float64, M)
    for v in y
        b = clamp(floor(Int, (v - ymin) / dy) + 1, 1, M)
        h[b] += 1
    end
    h ./= (yc[2] - yc[1]) * length(y)
    return collect(yc), h
end
