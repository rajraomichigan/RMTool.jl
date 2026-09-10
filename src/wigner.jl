"""
    wigner(n, isreal = true)

The Wigner matrix: generates an `n × n` symmetric (Hermitian) Wigner matrix
`W = (G + G')/sqrt(2n)`. If `isreal` is `true` the elements are real, otherwise
complex. The empirical spectral measure converges to the semicircle law
encoded by [`wignerpol`](@ref).

References:
[1] Alan Edelman, Handout 3: Experiments with Classical Matrix Ensembles, Course Notes 18.338, Fall 2004.
[2] E. P. Wigner, Characteristic vectors of bordered matrices with infinite dimensions, Ann. Math. 62 (1955).
"""
function wigner(n::Integer, isreal::Union{Bool,Integer} = true)
    g = Bool(isreal) ? randn(n, n) : (randn(n, n) + im * randn(n, n)) / sqrt(2)
    (g + g') / sqrt(2n)
end
