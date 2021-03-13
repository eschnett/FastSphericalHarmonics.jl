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
