"""
    PdfInfo

Result of [`Lmz2pdf`](@ref). Fields:

* `range`          – points at which the algebraic curve was evaluated
* `rr`             – all (complex) roots in `m` at every point of `range` (one row per point)
* `density`        – normalized imaginary parts `Im(root)/π` of the roots; one column
                     when `deg_m ≤ 3` (the root with positive imaginary part), all
                     roots otherwise. One of the columns is the probability density.
* `real`           – real parts of the same roots (real part of the algebraic curve)
* `poles`          – locations of poles (real zeros of the leading coefficient in `m`)
* `multipleroots`  – locations of multiple points (real zeros of the discriminant), i.e.
                     the possible boundary points of the support
* `mean`           – first moment of the measure (`NaN` if it could not be computed)

`plot(pdfinfo.range, pdfinfo.density)` (with Plots.jl) plots the normalized
imaginary roots; one of them corresponds to the true density.
"""
struct PdfInfo
    range::Vector{Float64}
    rr::Matrix{ComplexF64}
    density::Matrix{Float64}
    real::Matrix{Float64}
    poles::Vector{Float64}
    multipleroots::Vector{Float64}
    mean::Float64
end

function Base.show(io::IO, p::PdfInfo)
    print(io, "PdfInfo: ", length(p.range), " points on [", isempty(p.range) ? "" : string(first(p.range), ", ", last(p.range)),
          "], ", size(p.density, 2), " root branch(es); poles = ", p.poles,
          ", multipleroots = ", round.(p.multipleroots; digits = 6), ", mean = ", p.mean)
end

"""
    Lmz2pdf(Lmz [, xx]; params = Dict(), tol = 1e-8) -> PdfInfo

Plots the algebraic curve of equation `Lmz(m, z) = 0`: returns the real and
imaginary roots of `Lmz` (in `m`) evaluated over the points in the vector `xx`
(see [`PdfInfo`](@ref)).

If `xx` is omitted the range is chosen automatically from the poles and the
possible boundary points (a warning describes the range used). Any symbolic
parameter left in `Lmz` (e.g. `c`) is set to the value given in `params`
(`params = Dict(c => 1//3)`), or to `1/2` with a warning.

Comment: when the degree of the polynomial in `m` is greater than 3, additional
intervention might be needed to plot the "right root" — all roots are
returned. See the reference for tips on how to isolate it.

```julia
LmzA = wishartpol(1//2)                 # Wishart with parameter 0.5
pdfinfo = Lmz2pdf(LmzA, -3:0.01:3)
# using Plots; plot(pdfinfo.range, pdfinfo.density)
pdfinfo.multipleroots
pdfinfo.mean
```

Reference: N. Raj Rao and Alan Edelman, "The polynomial method for random matrices".
"""
function Lmz2pdf(Lmz, xx = nothing; params = Dict(), tol = 1e-8)
    L = numden(Lmz)
    im_ = _idx(m); iz = _idx(z)

    # --- fix free parameters
    for (v, val) in params
        L = subs(L, v => val)
    end
    extra = _params_in(L)
    if !isempty(extra)
        @warn "Input argument Lmz has symbolic variables other than m and z. Setting extra argument(s) $(extra) to 1/2"
        for v in extra
            L = subs(L, v => 1 // 2)
        end
    end
    other = [v for v in vars(L) if !(_idx(v) in (im_, iz))]
    isempty(other) || error("Lmz2pdf: polynomial depends on non-(m,z) variables $(other)")
    degree(L, im_) >= 1 && degree(L, iz) >= 0 || @warn "Input argument Lmz = $(L) is not a polynomial in m and z"

    Dm = degree(L, im_)
    lDmz = coeff(L, [im_], [Dm])
    poles = degree(lDmz, iz) == 0 ? Float64[] : _realroots(polyroots(_ucoeffs_f64(lDmz, iz)); tol = tol)

    discrimLmz = _discriminant(L, im_)
    if degree(discrimLmz, iz) <= 0
        PossBdryPts = Float64[]
        @warn "Possibly atomic measure encoded by the bivariate polynomial. Examine pdfinfo.poles and pdfinfo.real."
    else
        PossBdryPts = _realroots(polyroots(_ucoeffs_f64(discrimLmz, iz)); tol = tol)
    end

    # --- first moment (best effort)
    Moment1 = try
        mom = Lmz2MomS(L, 1; warn = false)
        mom isa AbstractVector && length(mom) >= 2 ? Float64(mom[2]) : NaN
    catch
        NaN
    end

    # --- range
    if xx === nothing
        pts = vcat(PossBdryPts, poles)
        if isempty(pts)
            xmin, xmax = -3.0, 3.0
        else
            xmin = minimum(pts) - 1; xmax = maximum(pts) + 1
        end
        step_size = min((xmax - xmin) / 100, 0.05)
        xx = collect(xmin:step_size:xmax)
        @warn "Region of support not specified. Range used between $(xmin) and $(xmax) in $(step_size) step increments. " *
              "Inspect pdfinfo.mean and pdfinfo.multipleroots to refine the interval. " *
              "Note: if pdfinfo.mean returned is non-sensical, then the density might be unbounded."
    else
        xx = collect(Float64, xx)
    end

    nx = length(xx)
    rr = zeros(ComplexF64, nx, Dm)
    ncol = Dm > 3 ? Dm : 1
    Dm > 3 && @warn "All roots returned -- isolate correct root manually"
    root = fill(NaN + 0im, nx, ncol)

    for (k, xv) in enumerate(xx)
        Lx = subs(L, z => rationalize(BigInt, xv))
        rts = polyroots(_ucoeffs_f64(Lx, im_))
        rr[k, 1:length(rts)] .= rts
        if Dm > 3
            root[k, 1:Dm] .= real.(rr[k, :]) .+ im .* imag.(rr[k, :]) ./ pi
        else
            pos = filter(rt -> imag(rt) > 0, rts)
            if !isempty(pos)
                rt = pos[argmax(imag.(pos))]
                root[k, 1] = real(rt) + im * imag(rt) / pi
            end
        end
    end

    PdfInfo(xx, rr, imag.(root), real.(root), poles, PossBdryPts, Moment1)
end
