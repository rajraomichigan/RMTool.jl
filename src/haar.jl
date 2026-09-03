"""
    haar(n)

An `n × n` Haar distributed unitary matrix (QR factorisation of a Gaussian
matrix followed by random phases).
"""
function haar(n::Integer)
    Q, _ = qr(randn(n, n))
    Matrix(Q) * Diagonal(exp.(2π * im * rand(n)))
end
