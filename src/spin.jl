# spinsph_mode

export spinsph_transform!
"""
    spinsph_transform!(F::Array{Complex{Float64},2}, s::Int)

Calculate the spin spherical harmonic transformation with spin weight
`s`. This is an in-place transform, i.e. the array `F` will be
overwritten by the coefficients. Use [`spinsph_transform`](@ref) for a
non-mutating function.

Use [`sph_points`](@ref) to caluclate the location of the points on
the sphere for the input array `F`.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_transform`](@ref), [`spinsph_evaluate!`](@ref),
[`sph_points`](@ref), [`sph_mode`](@ref)
"""
function spinsph_transform!(F::Array{Complex{Float64},2}, s::Int)
    N, M = size(F)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_spinsph2fourier(F, s)
    PA = plan_spinsph_analysis(F, s)
    C = F
    lmul!(PA, C)
    ldiv!(P, C)
    return C
end
export spinsph_transform
"""
    C = spinsph_transform!(F::Array{Complex{Float64},2}, s::Int)
    C::Array{Complex{Float64},2}

Calculate the spin spherical harmonic transformation with spin weight
`s`. You can use [`spinsph_transform!`](@ref) for more efficient a
mutating function that overwrites its argument `F`.

Use [`sph_points`](@ref) to caluclate the location of the points on
the sphere for the input array `F`.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_transform!`](@ref), [`spinsph_evaluate`](@ref),
[`sph_points`](@ref), [`sph_mode`](@ref)
"""
function spinsph_transform(F::Array{Complex{Float64},2}, s::Int)
    return spinsph_transform!(copy(F), s)
end

export spinsph_evaluate!
"""
    spinsph_evaluate!(C::Array{Complex{Float64},2}, s::Int)

Evaluate the spin spherical harmonic transformation with spin weight
`s` on points on the sphere. This is an in-place transform, i.e. the
array `C` will be overwritten by the point values. Use
[`spinsph_evaluate`](@ref) for a non-mutating function.

Use [`sph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`spinsph_evaluate`](@ref), [`spinsph_transform!`](@ref),
[`sph_mode`](@ref), [`sph_points`](@ref)
"""
function spinsph_evaluate!(C::Array{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_spinsph2fourier(C, s)
    PS = plan_spinsph_synthesis(C, s)
    F = C
    lmul!(P, F)
    lmul!(PS, F)
    return F
end
export spinsph_evaluate
"""
    F = spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)
    F::Array{Complex{Float64},2}

Evaluate the spin spherical harmonic transformation with spin weight
`s` on points on the sphere. You can use [`spinsph_evaluate!`](@ref)
for more efficient a mutating function that overwrites its argument
`C`.

Use [`sph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`spinsph_evaluate!`](@ref), [`spinsph_transform`](@ref),
[`sph_mode`](@ref), [`sph_points`](@ref)
"""
function spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)
    return spinsph_evaluate!(copy(C), s)
end

export spinsph_eth!
"""
    spinsph_eth!(C::Array{Complex{Float64},2}, s::Int)

Apply the differential operator ð ("eth") to the coefficients `C`.
This raises the spin weight `s` by 1. For a real function of spin
weight 0, this is equivalent to calculating the gradient.

This is an in-place transform, i.e. the array `C` will be overwritten
by the coefficients. Use [`spinsph_eth`](@ref) for a non-mutating
function.

This function is the converse of [`spinsph_ethbar!`](@ref), which is a
derivative operator lowering the spin weight.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_eth`](@ref), [`spinsph_ethbar!`](@ref),
[`sph_mode`](@ref)
"""
function spinsph_eth!(C::Array{Complex{Float64},2}, s::Int)
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

    ðC = C
    for l in 0:(lmax + mmax), m in 0:l
        if l - lmax ≤ m ≤ mmax
            if m == 0
                ðC[sph_mode(l, m)] *= sqrt((l - s) * (l + s + 1))
            else
                ðC[sph_mode(l, -m)] *= -sqrt((l - s) * (l + s + 1))
                ðC[sph_mode(l, +m)] *= +sqrt((l - s) * (l + s + 1))
            end
        end
    end

    return ðC
end
export spinsph_eth
"""
    ðC = spinsph_eth(C::Array{Complex{Float64},2}, s::Int)
    ðC::Array{Complex{Float64},2}

Apply the differential operator ð ("eth") to the coefficients `C`.
This raises the spin weight `s` by 1. For a real function of spin
weight 0, this is equivalent to calculating the gradient.

You can use [`spinsph_eth!`](@ref) for more efficient a mutating
function that overwrites its argument `C`.

This function is the converse of [`spinsph_ethbar`](@ref), which is a
derivative operator lowering the spin weight.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_eth!`](@ref), [`spinsph_ethbar`](@ref),
[`sph_mode`](@ref)
"""
spinsph_eth(C::Array{Complex{Float64},2}, s::Int) = spinsph_eth!(copy(C), s)

export spinsph_ethbar!
"""
    spinsph_ethbar!(C::Array{Complex{Float64},2}, s::Int)

Apply the differential operator ð̄ ("eth-bar") to the coefficients `C`.
This lowers the spin weight `s` by 1. For a function of spin weight 1,
this is equivalent to calculating the divergence.

This is an in-place transform, i.e. the array `C` will be overwritten
by the coefficients. Use [`spinsph_ethbar`](@ref) for a non-mutating
function.

This function is the converse of [`spinsph_eth!`](@ref), which is a
derivative operator raising the spin weight.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_ethvar`](@ref), [`spinsph_eth!`](@ref),
[`sph_mode`](@ref)
"""
function spinsph_ethbar!(C::Array{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ð̄C = C
    for l in 0:(lmax + mmax), m in 0:l
        if l - lmax ≤ m ≤ mmax
            if m == 0
                ð̄C[sph_mode(l, m)] *= -sqrt((l + s) * (l - s + 1))
            else
                ð̄C[sph_mode(l, -m)] *= +sqrt((l + s) * (l - s + 1))
                ð̄C[sph_mode(l, +m)] *= -sqrt((l + s) * (l - s + 1))
            end
        end
    end

    return ð̄C
end
export spinsph_ethbar
"""
    ðC = spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)
    ðC::Array{Complex{Float64},2}

Apply the differential operator ð̄ ("eth-bar") to the coefficients `C`.
This lowers the spin weight `s` by 1. For a function of spin weight 1,
this is equivalent to calculating the divergence.

You can use [`spinsph_ethbar!`](@ref) for more efficient a mutating
function that overwrites its argument `C`.

This function is the converse of [`spinsph_eth`](@ref), which is a
derivative operator raising the spin weight.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`spinsph_ethbar!`](@ref), [`spinsph_eth`](@ref),
[`sph_mode`](@ref)
"""
function spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)
    return spinsph_ethbar!(copy(C), s)
end
