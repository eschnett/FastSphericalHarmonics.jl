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

@testset "Laplacian" begin
    # lmax = 100
    lmax = 2

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)

    # Calculate Laplacian directly
    F = randn(N, M)
    C = sph_transform(F)
    ΔC = sph_laplace(C)
    ΔF = sph_evaluate(ΔC)

    # Calculate Laplacian via eth and ethbar
    F⁰ = Complex.(F)
    C⁰ = spinsph_transform(F⁰, 0)
    ðC¹ = spinsph_eth(C⁰, 0)
    ð̄ðC⁰ = spinsph_ethbar(ðC¹, 1)
    ð̄ðF⁰ = spinsph_evaluate(ð̄ðC⁰, 0)
    @test ð̄ðF⁰ ≈ ΔF

    # Calculate Laplacian via ethbar and eth
    ð̄C⁻¹ = spinsph_ethbar(C⁰, 0)
    ðð̄C⁰ = spinsph_eth(ð̄C⁻¹, -1)
    ðð̄F⁰ = spinsph_evaluate(ðð̄C⁰, 0)

    @test ðð̄F⁰ ≈ ΔF
end
