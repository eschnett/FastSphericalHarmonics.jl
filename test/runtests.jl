using FastTransforms
using FastSphericalHarmonics
using Test

# Integrate over the sphere
function integrate(u::Array{T1,2}, v::Array{T2,2}) where {T1,T2}
    @assert size(v) == size(u)
    N, M = size(u)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    Θ, Φ = sph_points(N)
    s = zero(one(eltype(conj(u))) * one(eltype(Θ)) * one(eltype(v)))
    for j in 1:M, i in 1:N
        s += conj(u[i, j]) * sin(Θ[i]) * v[i, j]
    end
    return 4π / (N * M) * s
end

@testset "Points on the sphere" begin
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)
    @test all(0 < θ < π for θ in Θ)
    @test all(0 ≤ ϕ < 2π for ϕ in Φ)
end

################################################################################

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
    lmax = 10

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
    F = ones(T, N, M)
    C = sph_transform(F)
    @test C[sph_mode(0, 0)] ≈ sqrt(4π)
    @test sum(abs2.(C)) ≈ 4π
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=0
    F = T[cos(θ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sph_transform(F)
    @test C[sph_mode(1, 0)] ≈ sqrt(4π / 3)
    @test sum(abs2.(C)) ≈ 4π / 3
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=-1
    F = T[sin(θ) * sin(ϕ) * sqrt(2) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C[sph_mode(1, -1)] ≈ sqrt(8π / 3)
    @test sum(abs2.(C)) ≈ 8π / 3
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test l=1, m=1
    F = T[sin(θ) * cos(ϕ) for θ in Θ, ϕ in Φ]
    C = sph_transform(F)
    @test C[sph_mode(1, 1)] ≈ sqrt(4π / 3)
    @test sum(abs2.(C)) ≈ 4π / 3
    F′ = sph_evaluate(C)
    @test F′ ≈ F

    # Test orthogonality
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
            if l == l′ && m == m′
                @test abs(int) > 1 / 2 - 1 / lmax^2
            else
                @test abs(int) < 2 / lmax^2
            end
        end
    end
end

################################################################################

@testset "Vector Spherical Harmonic mode indices" begin
    lmax = 10
    N = lmax + 1
    M = 2 * N - 1

    for v in 1:2
        mode_seen = falses(lmax + 1, 2 * lmax + 1)
        for l in 1:lmax, m in (-l):l
            lm = sphv_mode(l, m, v)
            @test !mode_seen[lm]
            mode_seen[lm] = true
        end
        @test all(mode_seen .≤ sphvones(Bool, N, M))
    end
end

@testset "Vector Spherical Harmonics for $T" for T in
                                                 [Float64, Complex{Float64}]
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Test duality
    F = randn(T, N, M)
    C = sphv_transform(F)
    @test size(C) == size(F)
    @test eltype(C) == eltype(F)
    F′ = sphv_evaluate(C)
    @test size(F′) == size(F)
    @test eltype(F′) == eltype(F)
    C′ = sphv_transform(F′)
    @test C′ ≈ C
    F″ = sphv_evaluate(C′)
    @test F″ ≈ F′

    # Test gradient field, e_θ, l=1, m=0
    # Test curl field, e_ϕ, l=1, m=0
    F = T[sin(θ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sphv_transform(F)
    @test C[1, 1] ≈ sqrt(8π / 3)
    @test sum(abs2.(C)) ≈ 8π / 3
    F′ = sphv_evaluate(C)
    @test F′ ≈ F

    # Test gradient field, e_θ, l=1, m=-1
    # Test curl field, e_ϕ, l=1, m=-1
    F = T[cos(θ) * sin(ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sphv_transform(F)
    @test C[2, 2] ≈ sqrt(2π / 3)
    @test sum(abs2.(C)) ≈ 2π / 3
    F′ = sphv_evaluate(C)
    @test F′ ≈ F

    # Test gradient field, e_ϕ, l=1, m=-1
    # Test curl field, e_θ, l=1, m=-1
    F = T[sin(ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sphv_transform(F)
    @test C[1, 2] ≈ sqrt(2π)
    @test sum(abs2.(C)) ≈ 2π
    F′ = sphv_evaluate(C)
    @test F′ ≈ F

    # Test gradient field, e_θ, l=1, m=1
    # Test curl field, e_ϕ, l=1, m=1
    F = T[cos(θ) * cos(ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sphv_transform(F)
    @test C[2, 3] ≈ sqrt(2π / 3)
    @test sum(abs2.(C)) ≈ 2π / 3
    F′ = sphv_evaluate(C)
    @test F′ ≈ F

    # Test gradient field, e_ϕ, l=1, m=1
    # Test curl field, e_θ, l=1, m=1
    F = T[cos(ϕ) for θ in Θ, ϕ in Φ]
    @test size(F) == (N, M)
    C = sphv_transform(F)
    @test C[1, 3] ≈ sqrt(2π)
    @test sum(abs2.(C)) ≈ 2π
    F′ = sphv_evaluate(C)
    @test F′ ≈ F

    # Test orthogonality
    lmax_test = 4
    for l in 1:lmax_test, m in (-l):l
        Cθ = zeros(T, N, M)
        Cϕ = zeros(T, N, M)
        Cθ[sphv_mode(l, m, 1)] = 1
        Cϕ[sphv_mode(l, m, 2)] = 1
        Fθ = sphv_evaluate(Cθ)
        Fϕ = sphv_evaluate(Cϕ)
        for l′ in 1:lmax_test, m′ in (-l′):l′
            Cθ′ = zeros(T, N, M)
            Cϕ′ = zeros(T, N, M)
            Cθ′[sphv_mode(l′, m′, 1)] = 1
            Cϕ′[sphv_mode(l′, m′, 2)] = 1
            Fθ′ = sphv_evaluate(Cθ′)
            Fϕ′ = sphv_evaluate(Cϕ′)
            int = integrate(Fθ, Fθ′) + integrate(Fϕ, Fϕ′)
            if l == l′ && m == m′
                @test abs(int) > 0.1
            else
                @test abs(int) < 0.1
            end
        end
    end
end

################################################################################

@testset "Spin Spherical Harmonics s=$s" for s in -2:2
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Test duality
    F = randn(Complex{Float64}, N, M)
    C = spinsph_transform(F, s)
    @test size(C) == size(F)
    @test eltype(C) == eltype(F)
    F′ = spinsph_evaluate(C, s)
    @test size(F′) == size(F)
    @test eltype(F′) == eltype(F)
    C′ = spinsph_transform(F′, s)
    @test C′ ≈ C
    F″ = spinsph_evaluate(C′, s)
    @test F″ ≈ F′
end
