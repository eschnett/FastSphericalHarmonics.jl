# spinsph_mode

export spinsph_transform
"""
    spinsph_transform(F::Array{Complex{Float64},2}, s::Int)

Calculate the spin spherical harmonic transformation with spin weight
`s`.
"""
function spinsph_transform(F::Array{Complex{Float64},2}, s::Int)
    N, M = size(F)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_spinsph2fourier(F, s)
    PA = plan_spinsph_analysis(F, s)
    C = P \ (PA * F)
    return C
end

export spinsph_evaluate
"""
    spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)

Evaluate the spin spherical harmonic transformation with spin weight
`s` on points on the sphere.
"""
function spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_spinsph2fourier(C, s)
    PS = plan_spinsph_synthesis(C, s)
    F = PS * (P * C)
    return F
end

export spinsph_eth
"""
    spinsph_eth(C::Array{Complex{Float64},2}, s::Int)

spinsph_eth raises the spin weight by one.

Note: ð of spin-weight 0 is gradient.
"""
function spinsph_eth(C::Array{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    # ðC = zero(C)
    # for l in 1:lmax
    #     ðC[l, 1] = sqrt(l * (l + 1)) * C[l + 1, 1]
    # end
    # for m in 1:mmax, l in 0:lmax
    #     ðC[l + 1, 2m + 0] = -sqrt((l + m) * (l + m + 1)) * C[l + 1, 2m + 0]
    #     ðC[l + 1, 2m + 1] = +sqrt((l + m) * (l + m + 1)) * C[l + 1, 2m + 1]
    # end

    ðC = zero(C)
    for l in 1:(lmax + mmax), m in 0:l
        if l - lmax ≤ m ≤ mmax
            # sgn = m ≥ 0 ? 1 : -1
            if m == 0
                ðC[sph_mode(l, m)] = sqrt((l - s) * (l + s + 1)) * #= - 1 =#
                                     C[sph_mode(l, m)]
            else
                ðC[sph_mode(l, -m)] = -sqrt((l - s) * (l + s + 1)) *
                                      C[sph_mode(l, -m)]
                ðC[sph_mode(l, +m)] = +sqrt((l - s) * (l + s + 1)) *
                                      C[sph_mode(l, +m)]
            end
        end
    end

    return ðC
end

export spinsph_ethbar
"""
    spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)

spinsph_ethbar lowers the spin weight by one.

ð̄ of spin-weight 0 is curl?
ð̄ of spin-weight 1 is divergence?
"""
function spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ð̄C = zero(C)
    for l in 1:(lmax + mmax), m in 0:l
        if l - lmax ≤ m ≤ mmax
            # sgn = m ≥ 0 ? 1 : -1
            if m == 0
                ð̄C[sph_mode(l, m)] = -sqrt((l + s) * (l - s + 1)) *
                                      C[sph_mode(l, m)] #= - 1 =#
            else
                ð̄C[sph_mode(l, -m)] = +sqrt((l + s) * (l - s + 1)) *
                                       C[sph_mode(l, -m)]
                ð̄C[sph_mode(l, +m)] = -sqrt((l + s) * (l - s + 1)) *
                                       C[sph_mode(l, +m)]
            end
        end
    end

    return ð̄C
end
