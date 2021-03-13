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
    # lmax = 100
    lmax = 4
    atol = rtol = 4 / lmax^2

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

    # Test gradient field, l=1, m=0
    Fθ = T[grad_Y_1_0_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_1_0_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test size(Fθ) == (N, M)
    @test size(Fϕ) == (N, M)
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(1, 0, 1)] ≈ sqrt(2)
    @test sum(abs2.(Cθ)) ≈ 2
    @test sum(abs2.(Cϕ)) ≈ 0
    Cgrad = Cθ                  # l=0

    # Test curl field, l=1, m=0
    Fθ = T[curl_Y_1_0_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[curl_Y_1_0_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test sum(abs2.(Cθ)) ≈ 0
    @test Cϕ[sphv_mode(1, 0, 2)] ≈ -sqrt(2)
    @test sum(abs2.(Cϕ)) ≈ 2
    Ccurl = Cϕ                  # l=0

    # Test gradient field, l=1, m=1
    Fθ = T[grad_Y_1p1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_1p1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(1, 1, 1)] ≈ sqrt(1 / 2)
    @test Cϕ[sphv_mode(1, 1, 2)] ≈ sqrt(3 / 2)
    @test sum(abs2.(Cθ)) ≈ 1 / 2
    @test sum(abs2.(Cϕ)) ≈ 3 / 2

    # Test curl field, l=1, m=1
    Fθ = T[curl_Y_1p1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[curl_Y_1p1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    # TODO: Using converse `v` value
    @test Cθ[sphv_mode(1, 1, 2)] ≈ sqrt(3 / 2)
    @test Cϕ[sphv_mode(1, 1, 1)] ≈ -sqrt(1 / 2)
    @test sum(abs2.(Cθ)) ≈ 3 / 2
    @test sum(abs2.(Cϕ)) ≈ 1 / 2

    # Test gradient field, l=1, m=-1
    Fθ = T[grad_Y_1m1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_1m1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(1, -1, 1)] ≈ sqrt(1 / 2)
    @test Cϕ[sphv_mode(1, -1, 2)] ≈ sqrt(3 / 2)
    @test sum(abs2.(Cθ)) ≈ 1 / 2
    @test sum(abs2.(Cϕ)) ≈ 3 / 2

    # Test curl field, l=1, m=-1
    Fθ = T[curl_Y_1m1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[curl_Y_1m1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 1 * 2; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    # TODO: Using converse `v` value
    @test Cθ[sphv_mode(1, -1, 2)] ≈ sqrt(3 / 2)
    @test Cϕ[sphv_mode(1, -1, 1)] ≈ -sqrt(1 / 2)
    @test sum(abs2.(Cθ)) ≈ 3 / 2
    @test sum(abs2.(Cϕ)) ≈ 1 / 2

    # Test gradient field, l=2, m=0
    Fθ = T[grad_Y_2_0_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_2_0_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 2 * 3; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(2, 0, 1)] ≈ sqrt(6)
    @test sum(abs2.(Cθ)) ≈ 6
    @test sum(abs2.(Cϕ)) ≈ 0

    # Test gradient field, l=2, m=1
    Fθ = T[grad_Y_2p1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_2p1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 2 * 3; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(2, 1, 1)] ≈ sqrt(8 / 3)
    @test Cϕ[sphv_mode(2, 1, 2)] ≈ sqrt(10 / 4)
    @show chop.(Cθ)
    @show chop.(Cϕ)
    @test sum(abs2.(Cθ)) ≈ 8 / 3
    @test sum(abs2.(Cϕ)) ≈ 10 / 4
    @assert false

    # Test gradient field, l=3, m=0
    Fθ = T[grad_Y_3_0_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_3_0_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 3 * 4; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(3, 0, 1)] ≈ sqrt(12)
    @test sum(abs2.(Cθ)) ≈ 12
    @test sum(abs2.(Cϕ)) ≈ 0

    # Test gradient field, l=3, m=1
    Fθ = T[grad_Y_3p1_θ(θ, ϕ) for θ in Θ, ϕ in Φ]
    Fϕ = T[grad_Y_3p1_ϕ(θ, ϕ) for θ in Θ, ϕ in Φ]
    @test isapprox(integrate(Fθ, Fϕ, Fθ, Fϕ), 3 * 4; atol=atol, rtol=rtol)
    Cθ = sphv_transform(Fθ)
    Cϕ = sphv_transform(Fϕ)
    @test Cθ[sphv_mode(3, 1, 1)] ≈ sqrt(27 / 4)
    @test sum(abs2.(Cθ)) ≈ 27 / 4
    @test sum(abs2.(Cϕ)) ≈ 0

    # Test orthogonality
    lmax_test = 4
    for gc in 1:2, l in 1:lmax_test, m in (-l):l
        Cθ = zeros(T, N, M)
        Cϕ = zeros(T, N, M)
        if gc == 1
            Cθ[sphv_mode(l, m, 1)] = sqrt(1 / 2)
            Cϕ[sphv_mode(l, m, 2)] = sqrt(3 / 2)
        else
            Cθ[sphv_mode(l, m, 2)] = sqrt(3 / 2)
            Cϕ[sphv_mode(l, m, 1)] = -sqrt(1 / 2)
        end
        Fθ = sphv_evaluate(Cθ)
        Fϕ = sphv_evaluate(Cϕ)
        for gc′ in 1:2, l′ in 1:lmax_test, m′ in (-l′):l′
            gc′ == gc || continue
            Cθ′ = zeros(T, N, M)
            Cϕ′ = zeros(T, N, M)
            if gc′ == 1
                Cθ′[sphv_mode(l′, m′, 1)] = sqrt(1 / 2)
                Cϕ′[sphv_mode(l′, m′, 2)] = sqrt(3 / 2)
            else
                Cθ′[sphv_mode(l′, m′, 2)] = sqrt(3 / 2)
                Cϕ′[sphv_mode(l′, m′, 1)] = -sqrt(1 / 2)
            end
            Fθ′ = sphv_evaluate(Cθ′)
            Fϕ′ = sphv_evaluate(Cϕ′)
            int = integrate(Fθ, Fϕ, Fθ′, Fϕ′)
            δ = l * (l + 1) * (gc == gc′ && l == l′ && m == m′)
            if !isapprox(int, δ; atol=atol, rtol=rtol)
                @show gc gc′ l l′ m m′ int δ
                @assert false
            end
            @test isapprox(int, δ; atol=atol, rtol=rtol)
        end
    end
end
