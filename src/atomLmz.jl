"""
    atomLmz(masses, weights)

Bivariate polynomial for an atomic matrix `A` with atoms at the locations in
the vector `masses` with the associated `weights`.

```julia
LmzA = atomLmz([1, 2], [1//2, 1//2])
# an atomic matrix with half of its eigenvalues equal to 1 and the remaining equal to 2
```
"""
function atomLmz(masses::AbstractVector, weights::AbstractVector)
    length(masses) == length(weights) || error("masses and weights must have the same length")
    expr = FF(m)
    for (a, w) in zip(masses, weights)
        expr -= _ff(w) / (_ff(a) - FF(z))
    end
    numden(expr)
end
