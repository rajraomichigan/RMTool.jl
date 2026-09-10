"""
    wishart(n, M = n, isreal = true)

The Wishart matrix: generates an `n × n` Hermitian Wishart matrix `W = G*G'/M`
with `G` an `n × M` Gaussian matrix. If `isreal` is `true` the elements are
real, otherwise complex. The empirical spectral measure converges to the
Marčenko–Pastur law with parameter `c = n/M`, encoded by [`wishartpol`](@ref).

References:
[1] Alan Edelman, Handout 3: Experiments with Classical Matrix Ensembles, Course Notes 18.338, Fall 2004.
[2] J. Wishart, The generalized product moment distribution in samples from a normal multivariate population, Biometrika 20-A (1928).
"""
function wishart(n::Integer, M::Integer = n, isreal::Union{Bool,Integer} = true)
    g = Bool(isreal) ? randn(n, M) : (randn(n, M) + im * randn(n, M)) / sqrt(2)
    (g * g') / M
end
