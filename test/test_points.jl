@testset "Points on the sphere" begin
    lmax = 10

    N = lmax + 1
    Θ, Φ = sph_points(N)
    @test length(Θ) == N
    M = length(Φ)
    @test all(0 < θ < π for θ in Θ)
    @test all(0 ≤ ϕ < 2π for ϕ in Φ)
end
