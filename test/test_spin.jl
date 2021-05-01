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
    @test C ≈ unit(spinsph_mode(s, l, m), size(C))
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
        @test C ≈ unit(spinsph_mode(s, l, m), size(C))
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

@testset "Spin spherical harmonics: eth and ethbar (s=$s, $T)" for s in -2:2,
                                                                   T in
                                                                   [Complex{Float64}]

    RT = typeof(Real(one(T)))
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4

    for l in max(abs(s), abs(s + 1)):lmax_test, m in (-l):l
        F = T[ðsYlm(s, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        C = spinsph_transform(F, s + 1)
        α = ifelse(m ≥ -s, 1, -1) * sqrt((l - s) * (l + s + 1))
        @test C ≈ unit(α, spinsph_mode(s + 1, l, m), size(C))
        F′ = spinsph_evaluate(C, s + 1)
        @test F′ ≈ F
    end

    for l in max(abs(s), abs(s - 1)):lmax_test, m in (-l):l
        F = T[ð̄sYlm(s, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        C = spinsph_transform(F, s - 1)
        α = -ifelse(m > -s, 1, -1) * sqrt((l + s) * (l - s + 1))
        @test C ≈ unit(α, spinsph_mode(s - 1, l, m), size(C))
        F′ = spinsph_evaluate(C, s - 1)
        @test F′ ≈ F
    end

    for l in abs(s):lmax_test, m in (-l):l
        C = zeros(T, N, M)
        C[spinsph_mode(s, l, m)] = ifelse(m ≥ -s, 1, -1)

        ðC = spinsph_eth(C, s)
        ðF = spinsph_evaluate(ðC, s + 1)
        @test ðF .+ 1 ≈ T[ðsYlm(s, l, m, θ, ϕ) for θ in Θ, ϕ in Φ] .+ 1
    end

    for l in abs(s):lmax_test, m in (-l):l
        C = zeros(T, N, M)
        C[spinsph_mode(s, l, m)] = ifelse(m > -s, 1, -1)

        ð̄C = spinsph_ethbar(C, s)
        ð̄F = spinsph_evaluate(ð̄C, s - 1)
        @test ð̄F .+ 1 ≈ T[ð̄sYlm(s, l, m, θ, ϕ) for θ in Θ, ϕ in Φ] .+ 1
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
end

@testset "Spin spherical harmonics: grad ($T)" for T in [Float64]
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4

    for l in 0:lmax_test, m in (-l):l
        F = Complex{T}[sYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ðF₀ = Complex{T}[ðsYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ð̄F₀ = Complex{T}[ð̄sYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]

        C = spinsph_transform(F, 0)
        @test C ≈ unit(1, spinsph_mode(0, l, m), size(C))

        ðC = spinsph_eth(C, 0)
        ð̄C = spinsph_ethbar(C, 0)

        ðF = spinsph_evaluate(ðC, 1)
        ð̄F = spinsph_evaluate(ð̄C, -1)
        # Different sYlm modes need different signs:
        α = ifelse(m ≥ 0, 1, -1)
        ᾱ = ifelse(m > 0, 1, -1)

        if l == 0
            @test norm(ðF) ≤ sqrt(eps(T))
            @test norm(ð̄F) ≤ sqrt(eps(T))
        else
            @test α * ðF ≈ ðF₀
            @test ᾱ * ð̄F ≈ ð̄F₀
        end
    end

    for l in 0:lmax_test, m in (-l):l
        F = T[Ylm(l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ðF₀ = SVector{2,T}[c2a(ðYlm(l, m, θ, ϕ)) for θ in Θ, ϕ in Φ]
        ð̄F₀ = SVector{2,T}[c2a(ð̄Ylm(l, m, θ, ϕ)) for θ in Θ, ϕ in Φ]
        ðC₀ = spinsph_transform(ðF₀, 1)
        ð̄C₀ = spinsph_transform(ð̄F₀, -1)
        ðC₀::Array{SVector{2,T},2}
        ð̄C₀::Array{SVector{2,T},2}

        C = spinsph_transform(F, 0)
        C::Array{T,2}
        for l′ in 0:lmax, m′ in (-l′):l′
            # Allow ±m modes to mix
            if (l′, abs(m′)) ≠ (l, abs(m))
                @test norm(C[spinsph_mode(0, l′, m′)]) ≤ sqrt(eps(T))
            end
        end

        ðC = spinsph_eth(C, 0)
        ð̄C = spinsph_ethbar(C, 0)
        ðC::Array{SVector{2,T},2}
        ð̄C::Array{SVector{2,T},2}

        if l == 0
            @test norm(ðC) ≤ sqrt(eps(T))
            @test norm(ð̄C) ≤ sqrt(eps(T))
        else
            @test ðC ≈ ðC₀
            @test ð̄C ≈ ð̄C₀
        end

        ðF = spinsph_evaluate(ðC, 1)
        ð̄F = spinsph_evaluate(ð̄C, -1)
        ðF::Array{SVector{2,T},2}
        ð̄F::Array{SVector{2,T},2}
        if l == 0
            @test norm(ðF) ≤ sqrt(eps(T))
            @test norm(ð̄F) ≤ sqrt(eps(T))
        else
            @test ðF ≈ ðF₀
            @test ð̄F ≈ ð̄F₀
        end
    end
end

@testset "Spin spherical harmonics: div ($T)" for T in [Float64]
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    lmax_test = 4

    for l in 0:lmax_test, m in (-l):l
        ðF = Complex{T}[ðsYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ð̄F = Complex{T}[ð̄sYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ΔF₀ = Complex{T}[-l * (l + 1) * sYlm(0, l, m, θ, ϕ) for θ in Θ, ϕ in Φ]

        ðC = spinsph_transform(ðF, 1)
        ð̄C = spinsph_transform(ð̄F, -1)
        ð̄ðC = spinsph_ethbar(ðC, 1)
        ðð̄C = spinsph_eth(ð̄C, -1)

        ð̄ðF = spinsph_evaluate(ð̄ðC, 0)
        ðð̄F = spinsph_evaluate(ðð̄C, 0)
        # Different sYlm modes need different signs:
        α = ifelse(m ≥ 0, 1, -1)
        ᾱ = ifelse(m > 0, 1, -1)

        if l == 0
            @test norm(ð̄ðF) ≤ sqrt(eps(T))
            @test norm(ðð̄F) ≤ sqrt(eps(T))
        else
            @test α * ð̄ðF ≈ ΔF₀
            @test ᾱ * ðð̄F ≈ ΔF₀
        end
    end

    for l in 0:lmax_test, m in (-l):l
        ðF = SVector{2,T}[c2a(ðYlm(l, m, θ, ϕ)) for θ in Θ, ϕ in Φ]
        ð̄F = SVector{2,T}[c2a(ð̄Ylm(l, m, θ, ϕ)) for θ in Θ, ϕ in Φ]
        ΔF₀ = T[-l * (l + 1) * Ylm(l, m, θ, ϕ) for θ in Θ, ϕ in Φ]
        ΔC₀ = spinsph_transform(ΔF₀, 0)
        ΔC₀::Array{T,2}

        ðC = spinsph_transform(ðF, 1)
        ð̄C = spinsph_transform(ð̄F, -1)
        ðC::Array{SVector{2,T},2}
        ð̄C::Array{SVector{2,T},2}
        ð̄ðC = spinsph_ethbar(ðC, 1)
        ðð̄C = spinsph_eth(ð̄C, -1)
        ð̄ðC::Array{T,2}
        ðð̄C::Array{T,2}

        if l == 0
            @test norm(ð̄ðC) ≤ sqrt(eps(T))
            @test norm(ðð̄C) ≤ sqrt(eps(T))
        else
            @test ð̄ðC ≈ ΔC₀
            @test ðð̄C ≈ ΔC₀
        end

        ð̄ðF = spinsph_evaluate(ð̄ðC, 0)
        ðð̄F = spinsph_evaluate(ðð̄C, 0)
        ð̄ðF::Array{T,2}
        ðð̄F::Array{T,2}

        if l == 0
            @test norm(ð̄ðF) ≤ sqrt(eps(T))
            @test norm(ðð̄F) ≤ sqrt(eps(T))
        else
            @test ð̄ðF ≈ ΔF₀
            @test ðð̄F ≈ ΔF₀
        end
    end
end

@testset "Spin spherical harmonics: Laplacian ($T)" for T in [Float64]
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    for iter in 1:10
        F = randn(T, N, M)

        C = spinsph_transform(F, 0)
        ðC = spinsph_eth(C, 0)
        ð̄C = spinsph_ethbar(C, 0)
        ð̄ðC = spinsph_ethbar(ðC, 1)
        ðð̄C = spinsph_eth(ð̄C, -1)

        for l in 0:lmax, m in (-l):l
            if l == 0
                @test norm(ð̄ðC[spinsph_mode(0, l, m)]) ≤ sqrt(eps(T))
                @test norm(ðð̄C[spinsph_mode(0, l, m)]) ≤ sqrt(eps(T))
            else
                @test ð̄ðC[spinsph_mode(0, l, m)] ≈
                      -l * (l + 1) * C[spinsph_mode(0, l, m)]
                @test ðð̄C[spinsph_mode(0, l, m)] ≈
                      -l * (l + 1) * C[spinsph_mode(0, l, m)]
            end
        end

        ð̄ðF = spinsph_evaluate(ð̄ðC, 0)
        ðð̄F = spinsph_evaluate(ðð̄C, 0)

        C₀ = sph_transform(F)
        ΔC₀ = sph_laplace(C₀)
        ΔF₀ = sph_evaluate(ΔC₀)

        @test ð̄ðF ≈ ΔF₀
        @test ðð̄F ≈ ΔF₀
    end
end
