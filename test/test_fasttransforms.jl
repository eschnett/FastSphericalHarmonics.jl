using FastTransforms
using Test

@testset "FastTransforms spin spherical harmonics" begin
    lmax = 4
    s = 1

    C1 = zeros(Complex{Float64}, lmax + 1, 2lmax + 1)
    C1[1, 2] = 1                 # l=1, m=-1

    C2 = zeros(Complex{Float64}, lmax + 1, 2lmax + 1)
    C2[2, 2] = 1                 # l=2, m=-1

    # Evalute
    P = plan_spinsph2fourier(C1, s)
    PS = plan_spinsph_synthesis(C1, s)
    F1 = PS * (P * C1)
    F2 = PS * (P * C2)

    # Calculate scalar product kernel
    F = conj(F1) .* F2

    # Integrate (by calculating l=0 mode)
    P = plan_sph2fourier(C1)
    PA = plan_sph_analysis(F)
    C = P \ (PA * F)

    # Should be zero
    @test abs(C[1, 1]) < 10eps()
end
