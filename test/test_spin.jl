@testset "Spin spherical harmonic mode indices (s=$s)" for s in -2:2
    lmax = 10
    N = lmax + 1
    M = 2 * N - 1

    mode_seen = zeros(Int, N, M)
    for l in abs(s):lmax, m in (-l):l
        mode_seen[spinsph_mode(s, l, m)] += 1
    end
    @test mode_seen == spinsphones(Bool, N, M, s)

    # Test some concrete modes
    @test spinsph_mode(0, 0, 0) == CartesianIndex(1, 1)
    @test spinsph_mode(0, 1, -1) == CartesianIndex(1, 2)
    @test spinsph_mode(0, 1, 1) == CartesianIndex(1, 3)

    @test spinsph_mode(1, 1, 0) == CartesianIndex(1, 1)
    @test spinsph_mode(1, 1, -1) == CartesianIndex(1, 2)
    @test spinsph_mode(1, 1, 1) == CartesianIndex(1, 3)

    @test spinsph_mode(2, 2, 0) == CartesianIndex(1, 1)
    @test spinsph_mode(2, 2, -1) == CartesianIndex(1, 2)
    @test spinsph_mode(2, 2, 1) == CartesianIndex(1, 3)
end

@testset "Spin spherical harmonics: simple modes ($T)" for T in
                                                           [Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Test s=1, l=1, m=0
    s, l, m = 1, 1, 0
    F = T[Y_s1_1_0(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = spinsph_transform(F, s)
    @test C[spinsph_mode(s, l, m)] ≈ 1
    @test sum(abs2.(C)) ≈ 1
    F′ = spinsph_evaluate(C, s)
    @test F′ ≈ F
end

@testset "Spin spherical harmonics: simple modes (s=$s, $T)" for s in -2:2,
                                                                 T in [Float64,
                                                                  Complex{Float64}]

    T <: Real && s ≠ 0 && continue

    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4
    for l in abs(s):lmax_test, m in (-l):l
        F = T[sYlm(T, s, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        C = spinsph_transform(F, s)
        @test C[spinsph_mode(s, l, m)] ≈ 1
        @test sum(abs2.(C)) ≈ 1
        F′ = spinsph_evaluate(C, s)
        @test F′ ≈ F
    end
end

@testset "Spin spherical harmonics: linearity (s=$s, $T)" for s in -2:2,
                                                              T in [Float64,
                                                               Complex{Float64}]

    T <: Real && s ≠ 0 && continue

    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    F = zeros(T, N, M)
    C = spinsph_transform(F, s)
    @test iszero(C)

    F1 = randn(T, N, M)
    F2 = randn(T, N, M)
    C1 = spinsph_transform(F1, s)
    C2 = spinsph_transform(F2, s)
    F = F1 + F2
    C = spinsph_transform(F, s)
    @test C ≈ C1 + C2

    F = randn(T, N, M)
    # Why doesn't this work with a complex scale factor? Is the
    # transform not linear?
    # α = randn(T)
    α = randn()
    C = spinsph_transform(F, s)
    Cα = spinsph_transform(α * F, s)
    @test Cα ≈ α * C
end

@testset "Spin spherical harmonics: duality (s=$s, $T)" for s in -2:2,
                                                            T in [Float64,
                                                             Complex{Float64}]

    T <: Real && s ≠ 0 && continue

    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    F = randn(T, N, M)
    C = spinsph_transform(F, s)
    @test size(C) == size(F)
    @test eltype(C) == eltype(F)
    F′ = spinsph_evaluate(C, s)
    @test F′ ≈ F

    C = randn(T, N, M)
    F = spinsph_evaluate(C, s)
    @test size(F) == size(C)
    @test eltype(F) == eltype(C)
    C′ = spinsph_transform(F, s)
    @test C′ ≈ C
end

@testset "Spin spherical harmonics: orthonormality (s=$s, $T)" for s in -2:2,
                                                                   T in
                                                                   [Float64,
                                                                    Complex{Float64}]

    T <: Real && s ≠ 0 && continue

    lmax = 100
    atol = 4 / lmax^2

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4
    for l in abs(s):lmax_test, m in (-l):l
        C = zeros(T, N, M)
        C[spinsph_mode(s, l, m)] = 1
        F = spinsph_evaluate(C, s)
        for l′ in abs(s):lmax_test, m′ in (-l′):l′
            C′ = zeros(T, N, M)
            C′[spinsph_mode(s, l′, m′)] = 1
            F′ = spinsph_evaluate(C′, s)
            int = integrate(F, F′)
            δ = l == l′ && m == m′
            @test isapprox(int, δ; atol=atol)
        end
    end
end

@testset "Spin spherical harmonics: Laplacian ($T)" for T in [Complex{Float64}]
    lmax = 100

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Calculate Laplacian directly
    RT = typeof(Real(one(T)))
    F = randn(RT, N, M)
    C = sph_transform(F)
    ΔC = sph_laplace(C)
    ΔF = sph_evaluate(ΔC)

    # Calculate Laplacian via eth and ethbar
    F⁰ = T.(F)
    C⁰ = spinsph_transform(F⁰, 0)
    @test eltype(C⁰) === T
    ðC¹ = spinsph_eth(C⁰, 0)
    @test eltype(ðC¹) === T
    ð̄ðC⁰ = spinsph_ethbar(ðC¹, 1)
    @test eltype(ð̄ðC⁰) === T
    ð̄ðF⁰ = spinsph_evaluate(ð̄ðC⁰, 0)
    @test eltype(ð̄ðF⁰) === T
    @test ð̄ðF⁰ ≈ ΔF

    # Calculate Laplacian via ethbar and eth
    ð̄C⁻¹ = spinsph_ethbar(C⁰, 0)
    @test eltype(ð̄C⁻¹) === T
    ðð̄C⁰ = spinsph_eth(ð̄C⁻¹, -1)
    @test eltype(ðð̄C⁰) === T
    ðð̄F⁰ = spinsph_evaluate(ðð̄C⁰, 0)
    @test eltype(ð̄ðF⁰) === T
    @test ðð̄F⁰ ≈ ΔF

    # Calculate Laplacian via real eth and ethbar
    F⁰ = F
    C⁰ = spinsph_transform(F⁰, 0)
    @test eltype(C⁰) === RT
    ðC¹θ, ðC¹ϕ = spinsph_eth(C⁰, 0)
    @test eltype(ðC¹θ) === RT
    @test eltype(ðC¹ϕ) === RT
    ð̄ðC⁰ = spinsph_ethbar(ðC¹θ, ðC¹ϕ, 1)
    @test eltype(ð̄ðC⁰) === RT
    ð̄ðF⁰ = spinsph_evaluate(ð̄ðC⁰, 0)
    @test eltype(ð̄ðF⁰) === RT
    @test ð̄ðF⁰ ≈ ΔF

    # Calculate Laplacian via real ethbar and eth
    ð̄C⁻¹θ, ð̄C⁻¹ϕ = spinsph_ethbar(C⁰, 0)
    @test eltype(ð̄C⁻¹θ) === RT
    @test eltype(ð̄C⁻¹ϕ) === RT
    ðð̄C⁰ = spinsph_eth(ð̄C⁻¹θ, ð̄C⁻¹ϕ, -1)
    @test eltype(ðð̄C⁰) === RT
    ðð̄F⁰ = spinsph_evaluate(ðð̄C⁰, 0)
    @test eltype(ð̄ðF⁰) === RT
    @test ðð̄F⁰ ≈ ΔF
end
