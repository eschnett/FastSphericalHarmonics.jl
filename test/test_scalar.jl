@testset "Scalar spherical harmonic mode indices" begin
    lmax = 10
    N = lmax + 1
    M = 2 * N - 1

    mode_seen = zeros(Int, N, M)
    for l in 0:lmax, m in (-l):l
        mode_seen[sph_mode(l, m)] += 1
    end
    @test mode_seen == sphones(Bool, N, M)
end

@testset "Scalar spherical harmonics: simple modes ($T)" for T in [Float64,
                                                                   Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Test l=0, m=0
    F = T[Y_0_0(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C ≈ unit(sph_mode(0, 0), size(C))
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=0
    F = T[Y_1_0(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sph_transform(F)
    @test C ≈ unit(sph_mode(1, 0), size(C))
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=1
    F = T[Y_1p1(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C ≈ unit(sph_mode(1, 1), size(C))
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=-1
    F = T[Y_1m1(θ, ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C ≈ unit(sph_mode(1, -1), size(C))
    F′ = sph_evaluate(C)
    @test F′ ≈ F
end

@testset "Scalar spherical harmonics: simple modes ($T)" for T in [Float64,
                                                                   Complex{Float64}]
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4
    for l in 0:lmax_test, m in (-l):l
        F = T[Ylm(l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        C = sph_transform(F)
        @test C ≈ unit(sph_mode(l, m), size(C))
        F′ = sph_evaluate(C)
        @test F′ ≈ F
    end
end

@testset "Scalar spherical harmonics: linearity ($T)" for T in [Float64,
                                                                Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    F = zeros(T, N, M)
    C = sph_transform(F)
    @test iszero(C)

    F1 = randn(T, N, M)
    F2 = randn(T, N, M)
    C1 = sph_transform(F1)
    C2 = sph_transform(F2)
    F = F1 + F2
    C = sph_transform(F)
    @test C ≈ C1 + C2

    F = randn(T, N, M)
    α = randn(T)
    C = sph_transform(F)
    Cα = sph_transform(α * F)
    @test Cα ≈ α * C
end

@testset "Scalar spherical harmonics: duality ($T)" for T in [Float64,
                                                              Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    F = randn(T, N, M)
    C = sph_transform(F)
    @test size(C) == size(F)
    @test eltype(C) == eltype(F)
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    C = randn(T, N, M)
    F = sph_evaluate(C)
    @test size(F) == size(C)
    @test eltype(F) == eltype(C)
    C′ = sph_transform(F)
    @test C′ ≈ C
end

@testset "Scalar spherical harmonics: orthonormality ($T)" for T in [Float64,
                                                                     Complex{Float64}]
    lmax = 100
    atol = 4 / lmax^2

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

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
end

@testset "Scalar spherical harmonics: Laplacian ($T)" for T in [Float64,
                                                                Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # z + 2x
    F = T[cos(θ) + 2 * sin(θ) * cos(ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    ΔC = sph_laplace(C)
    ΔF = sph_evaluate(ΔC)
    @test ΔF ≈ -2F
end
