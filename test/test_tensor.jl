@testset "Tensor representations (scalars)" begin
    T = Complex{Float64}

    lmax = 10
    N = lmax + 1
    Θ, Φ = sph_points(N)
    M = length(Φ)

    for fun in
        [(θ, ϕ) -> 1.0, (θ, ϕ) -> sin(θ) * cos(ϕ), (θ, ϕ) -> sin(θ) * sin(ϕ),
         (θ, ϕ) -> cos(θ)]
        F = Tensor{0}([fun(θ, ϕ) for θ in Θ, ϕ in Φ])
        C = SpinTensor{0}(F)
        F′ = Tensor{0}(C)
        C′ = SpinTensor{0}(F′)
        @test F′ ≈ F
        @test C′ ≈ C
    end
end

@testset "Tensor representations (vectors)" begin
    T = Complex{Float64}

    lmax = 10
    N = lmax + 1
    Θ, Φ = sph_points(N)
    M = length(Φ)

    for fun in
        [(θ, ϕ) -> 1.0, (θ, ϕ) -> sin(θ) * cos(ϕ), (θ, ϕ) -> sin(θ) * sin(ϕ),
         (θ, ϕ) -> cos(θ)]
        F⁰ = T[fun(θ, ϕ) for θ in Θ, ϕ in Φ]::AbstractArray{<:Complex}
        C⁰ = spinsph_transform(F⁰, 0)
        C¹ = spinsph_eth(C⁰, 0)
        F¹ = spinsph_evaluate(C¹, 1)::AbstractArray{<:Complex}

        F = Tensor{1}([SVector(real(f¹), imag(f¹)) for f¹ in F¹])
        C = SpinTensor{1}(F)
        @test C.coeffs[1] ≈ C¹
        @test C.coeffs[2] ≈ C¹ # because F⁰ is real
        F′ = Tensor{1}(C)
        C′ = SpinTensor{1}(F′)
        @test F′ ≈ F
        @test C′ ≈ C
    end
end

@testset "Tensor representations (tensors)" begin
    T = Complex{Float64}

    lmax = 10
    N = lmax + 1
    Θ, Φ = sph_points(N)
    M = length(Φ)

    m = SVector(1, im)
    m̄ = conj(m)

    for fun in
        [(θ, ϕ) -> 1.0, (θ, ϕ) -> sin(θ) * cos(ϕ), (θ, ϕ) -> sin(θ) * sin(ϕ),
         (θ, ϕ) -> cos(θ)]
        F⁰ = T[fun(θ, ϕ) for θ in Θ, ϕ in Φ]::AbstractArray{<:Complex}
        C⁰ = spinsph_transform(F⁰, 0)
        C¹ = spinsph_eth(C⁰, 0)
        C² = spinsph_eth(C¹, 1)
        F² = spinsph_evaluate(C², 2)::AbstractArray{<:Complex}

        F = Tensor{2}([SMatrix{2,2}(f² * m̄[a] * m̄[b] / 4
                                    for a in 1:2, b in 1:2) for f² in F²])
        C = SpinTensor{2}(F)
        @test C.coeffs[1, 1] ≈ C²
        @test isapprox(C.coeffs[1, 2], zero(C²); atol=1000eps()) # because F⁰ is real
        @test isapprox(C.coeffs[2, 1], zero(C²); atol=1000eps()) # because F⁰ is real
        @test isapprox(C.coeffs[2, 2], C²; atol=1000eps()) # because F⁰ is real
        F′ = Tensor{2}(C)
        C′ = SpinTensor{2}(F′)
        @test F′ ≈ F
        @test C′ ≈ C
    end
end

@testset "Tensor gradients" begin
    T = Complex{Float64}

    lmax = 2 #TODO 10
    N = lmax + 1
    Θ, Φ = sph_points(N)
    M = length(Φ)

    m = SVector(1, im)
    m̄ = conj(m)

    for fun in [(name="1", f=(θ, ϕ) -> 1, fθ=(θ, ϕ) -> 0, fϕ=(θ, ϕ) -> 0),
                (name="x", f=(θ, ϕ) -> sin(θ) * cos(ϕ), fθ=(θ, ϕ) -> cos(θ) * cos(ϕ),
                 fϕ=(θ, ϕ) -> -sin(ϕ)),
                (name="y", f=(θ, ϕ) -> sin(θ) * sin(ϕ), fθ=(θ, ϕ) -> cos(θ) * sin(ϕ),
                 fϕ=(θ, ϕ) -> cos(ϕ)),
                (name="z", f=(θ, ϕ) -> cos(θ), fθ=(θ, ϕ) -> -sin(θ), fϕ=(θ, ϕ) -> 0)]
        @show fun.name
        F⁰ = T[fun.f(θ, ϕ) for θ in Θ, ϕ in Φ]::AbstractArray{<:Complex}
        C⁰ = spinsph_transform(F⁰, 0)
        C¹ = spinsph_eth(C⁰, 0)
        C⁻¹ = spinsph_ethbar(C⁰, 0)
        F¹ = spinsph_evaluate(C¹, 1)
        F⁻¹ = spinsph_evaluate(C⁻¹, -1)

        F¹′ = T[-sqrt(3 / 16π) * (1 - cos(θ)) * cis(ϕ) for θ in Θ, ϕ in Φ]::AbstractArray{<:Complex}

        println("F⁰:")
        display(chop.(F⁰))
        println()
        println("C⁰:")
        display(chop.(C⁰))
        println()
        println("C¹:")
        display(chop.(C¹))
        println()
        println("F¹:")
        display(chop.(F¹))
        println()
        println("F¹′:")
        display(chop.(F¹′))
        println()
        println("F⁻¹:")
        display(chop.(F⁻¹))
        println()

        F = Tensor{0}(F⁰)
        C = SpinTensor{0}(F)
        dC = tensor_gradient(C)
        dF = Tensor{1}(dC)

        println("F:")
        display(chop.(map(x -> x[], F.values)))
        println()
        println("dF (θ,ϕ):")
        display(chop.(map(x -> x[1], dF.values)))
        println()
        display(chop.(map(x -> x[2], dF.values)))
        println()

        # expected = [SVector{2}(m̄[a] * f¹ + m[a] * f⁻¹ for a in 1:2)
        #             for (f¹, f⁻¹) in zip(F¹, F⁻¹)] / 2
        expected = [SVector{2,T}(fun.fθ(θ, ϕ), fun.fϕ(θ, ϕ)) for θ in Θ, ϕ in Φ]

        println("expected (θ,ϕ):")
        display(chop.(map(x -> x[1], expected)))
        println()
        display(chop.(map(x -> x[2], expected)))
        println()

        for a in 1:2
            # ES 2022-01-25: I think the complex spherical harmonics
            # in this package are not the "standard" complex spherical
            # harmonics.
            @test map(x -> x[a], dF.values) ≈ map(x -> x[a], expected)
        end

        # F = Tensor{1}([SVector(real(f¹), imag(f¹)) for f¹ in F¹])
        # C = SpinTensor{1}(F)
        # dC = tensor_gradient(C)
        # dF = Tensor{2}(dC)
    end
end
