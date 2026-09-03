using Test
using RMTool
using Nemo
using LinearAlgebra

# equality of polynomials up to a nonzero constant factor
equiv(a, b) = divides(a, b)[1] && divides(b, a)[1]

@testset "RMTool" begin
    @testset "polynomial representation" begin
        L1 = numden(m - rat(1//2)//(1 - z) - rat(1//2)//(2 - z))
        @test equiv(L1, atomLmz([1, 2], [1//2, 1//2]))
        @test equiv(wishartpol(1//2), subs(wishartpol(c), c => 1//2))
        @test size(TLmz(wishartpol(1//2))) == (3, 2)
    end

    @testset "moments" begin
        @test Lmz2MomS(wignerpol(), 6) == [1, 0, 1, 0, 2, 0, 5]
        # Narayana polynomials
        mom = Lmz2MomS(wishartpol(c), 4; numeric = false)
        @test mom[3] == FF(c + 1)
        @test mom[4] == FF(c^2 + 3c + 1)
        @test mom[5] == FF(c^3 + 6c^2 + 6c + 1)
        @test Lmz2MomS(wishartpol(1//2), 3) == [1, 1, 3//2, 11//4]
        @test Lmz2MomS(atomLmz([1, 2], [1//2, 1//2]), 3) == [1, 3//2, 5//2, 9//2]
    end

    @testset "transforms round trip" begin
        b = wishartpol(1//2)
        @test equiv(Lrg2Lmz(Lmz2Lrg(b)), b)
        @test equiv(Lsy2Lmz(Lmz2Lsy(b)), b)
        @test equiv(Lmyuz2Lmz(Lmz2Lmyuz(b)), b)
        @test equiv(Letaz2Lmz(Lmz2Letaz(b)), b)
        @test equiv(Lgz2Lmz(Lmz2Lgz(b)), b)
    end

    @testset "free addition" begin
        # two atoms at ±1 -> arc-sine law
        Aa = atomLmz([1, -1], [1//2, 1//2])
        @test equiv(AplusB(Aa, Aa), m^2 * z^2 - 4m^2 - 1)
        S = AplusB(wignerpol(), wishartpol(1//2))
        @test equiv(S, m^3 + m^2 * z + 2m^2 + 2m * z - m + 2)
        @test Lmz2MomS(S, 3) == [1, 1, 5//2, 23//4]
        @test equiv(addAdtimes(wignerpol(), 1), wignerpol())
    end

    @testset "free multiplication and Wishart" begin
        A0 = atomLmz([1, 0], [1//2, 1//2])
        P = AtimesB(A0, A0)
        @test equiv(P, 4m^2 * z^3 - 4m^2 * z^2 + 4m * z^2 - 4m * z - 1)
        @test equiv(transposeA(P, 2), m^2 * z^2 - m^2 * z - 1)
        @test equiv(AtimesWish(atomLmz([1], [1]), c), wishartpol(c))
        @test equiv(AtimesWish(atomLmz([1], [1]), 1//2), wishartpol(1//2))
    end

    @testset "deterministic operations" begin
        L = atomLmz([1, 2], [1//2, 1//2])
        @test equiv(invA(L), atomLmz([1, 1//2], [1//2, 1//2]))
        @test equiv(shiftA(L, -1), atomLmz([0, 1], [1//2, 1//2]))
        @test equiv(scaleA(L, 2), atomLmz([2, 4], [1//2, 1//2]))
        @test equiv(squareA(L), atomLmz([1, 4], [1//2, 1//2]))
        @test equiv(squareA(wignerpol()), wishartpol(1))
        @test Lmz2MomS(squareA(wignerpol()), 4) == [1, 1, 2, 5, 14]
        @test equiv(AblockB(L, L, 1//2), L)
        @test equiv(equiLmz(1, 1), m^2 + 2m * z + 2)
    end

    @testset "density" begin
        info = Lmz2pdf(wishartpol(1//2), -0.05:0.05:5)
        @test length(info.range) == length(-0.05:0.05:5)
        @test size(info.density, 2) == 1
        # support of the Marcenko-Pastur law with c = 1/2: [(1-sqrt(c))^2, (1+sqrt(c))^2]
        @test isapprox(sort(info.multipleroots), [(1 - sqrt(0.5))^2, (1 + sqrt(0.5))^2]; atol = 1e-6)
        @test info.mean ≈ 1
        dx = 0.05
        @test isapprox(sum(filter(!isnan, info.density[:, 1])) * dx, 1; atol = 0.05)
        @test all(x -> isnan(x) || x >= -1e-12, info.density)
    end

    @testset "numerics" begin
        n = 50
        W = wigner(n); @test W ≈ W'
        Wc = wigner(n, false); @test Wc ≈ Wc'
        S = wishart(n, 2n); @test S ≈ S' && all(eigvals(S) .>= -1e-10)
        Q = haar(n); @test Q * Q' ≈ Matrix{ComplexF64}(I, n, n)
        yc, h = histw(randn(10_000), 40)
        @test length(yc) == 40 && isapprox(sum(h) * (yc[2] - yc[1]), 1; atol = 1e-9)
        d = theoryvssim(randn(100), 10, [0.0], [0.0])
        @test length(d.centers) == 10
    end
end
