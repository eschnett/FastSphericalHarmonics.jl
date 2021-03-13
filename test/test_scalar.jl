@testset "Scalar Spherical Harmonic mode indices" begin
    lmax = 10
    N = lmax + 1
    M = 2 * N - 1

    mode_seen = falses(lmax + 1, 2lmax + 1)
    for l in 0:lmax, m in (-l):l
        lm = sph_mode(l, m)
        @test !mode_seen[lm]
        mode_seen[lm] = true
    end
    for l in 0:lmax
        for m in 0:(2 * (lmax - l))
            @test mode_seen[l + 1, m + 1]
        end
        for m in (2 * (lmax - l) + 1):(2 * lmax)
            @test !mode_seen[l + 1, m + 1]
        end
    end
    @test mode_seen == sphones(Bool, N, M)
end

@testset "Scalar Spherical Harmonics for $T" for T in
                                                 [Float64, Complex{Float64}]
    lmax = 100
    atol = 4 / lmax^2

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Test duality
    F = randn(T, N, M)
    C = sph_transform(F)
    @test size(C) == size(F)
    @test eltype(C) == eltype(F)
    F′ = sph_evaluate(C)
    @test size(F′) == size(F)
    @test eltype(F′) == eltype(F)
    C′ = sph_transform(F′)
    @test C′ ≈ C
    F″ = sph_evaluate(C′)
    @test F″ ≈ F′

    # Test l=0, m=0
    F = T[Y_0_0(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C[sph_mode(0, 0)] ≈ 1
    @test sum(abs2.(C)) ≈ 1
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=0
    F = T[Y_1_0(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sph_transform(F)
    @test C[sph_mode(1, 0)] ≈ 1
    @test sum(abs2.(C)) ≈ 1
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=1
    F = T[Y_1p1(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C[sph_mode(1, 1)] ≈ 1
    @test sum(abs2.(C)) ≈ 1
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=-1
    F = T[Y_1m1(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C[sph_mode(1, -1)] ≈ 1
    @test sum(abs2.(C)) ≈ 1
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test orthonormality
    lmax_test = 4
    for l in 0:lmax_test, m in (-l):l
        C = zeros(T, N, M)
        C[sph_mode(l, m)] = 1
        F = sph_evaluate(C)
        for l′ in 0:lmax_test, m′ in (-l′):l′
            C′ = zeros(T, N, M)
            C′[sph_mode(l′, m′)] = 1
            F′ = sph_evaluate(C′)
            int = integrate(F, F′)
            δ = l == l′ && m == m′
            @test isapprox(int, δ; atol=atol)
        end
    end

    # Test Laplacian
    F = [cos(θ) + sin(θ) * cos(ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    ΔC = sph_laplace(C)
    ΔF = sph_evaluate(ΔC)
    @test ΔF ≈ -2F
end
