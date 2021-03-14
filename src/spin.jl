export spinsph_mode
"""
    idx = spinsph_mode(l::Integer, m::Integer, s::Integer)
    idx::CartesianIndex{2}

Calculate the Cartesian index `idx` for the `l`,`m` mode for spin
weight `s`. This index can be used to access the coefficients, i.e.
the result of [`spinsph_transform`](@ref) or the input to
[`spinsph_evaluate`](@ref).

Coefficients are stored in a two-dimensional array. Not all array
elements are used. See [this
page](https://mikaelslevinsky.github.io/FastTransforms/transforms.html),
section "spinsph2fourier", for details.

See also: [`spinsph_transform!`](@ref), [`spinsph_transform`](@ref),
[`spinsph_evaluate!`](@ref), [`spinsph_evaluate`](@ref)
"""
function spinsph_mode(s::Int, l::Int, m::Int)
    @assert l ≥ abs(s)
    @assert abs(m) ≤ l
    return CartesianIndex(l - max(abs(s), abs(m)) + 1, 2 * abs(m) + (m ≥ 0))
end
function spinsph_mode(s::Integer, l::Integer, m::Integer)
    return spinsph_mode(Int(s), Int(l), Int(m))
end

################################################################################

export spinsph_transform!
"""
    spinsph_transform!(F::Array{Complex{Float64},2}, s::Int)

Calculate the spin spherical harmonic transformation with spin weight
`s`. This is an in-place transform, i.e. the array `F` will be
overwritten by the coefficients. Use [`spinsph_transform`](@ref) for a
non-mutating function.

Use [`sph_points`](@ref) to caluclate the location of the points on
the sphere for the input array `F`.

Use [`spinsph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode with spin weight `s`.

See also: [`spinsph_transform`](@ref), [`spinsph_evaluate!`](@ref),
[`sph_points`](@ref), [`spinsph_mode`](@ref)
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

Use [`spinsph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode with spin weight `s`.

See also: [`spinsph_transform!`](@ref), [`spinsph_evaluate`](@ref),
[`sph_points`](@ref), [`spinsph_mode`](@ref)
"""
function spinsph_transform(F::Array{Complex{Float64},2}, s::Int)
    return spinsph_transform!(copy(F), s)
end

################################################################################

export spinsph_evaluate!
"""
    spinsph_evaluate!(C::Array{Complex{Float64},2}, s::Int)

Evaluate the spin spherical harmonic transformation with spin weight
`s` on points on the sphere. This is an in-place transform, i.e. the
array `C` will be overwritten by the point values. Use
[`spinsph_evaluate`](@ref) for a non-mutating function.

Use [`spinsph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode with spin weight `s`.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`spinsph_evaluate`](@ref), [`spinsph_transform!`](@ref),
[`spinsph_mode`](@ref), [`sph_points`](@ref)
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

Use [`spinsph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode with spin weight `s`.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`spinsph_evaluate!`](@ref), [`spinsph_transform`](@ref),
[`spinsph_mode`](@ref), [`sph_points`](@ref)
"""
function spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)
    return spinsph_evaluate!(copy(C), s)
end

################################################################################

export spinsph_eth
"""
    ðC = spinsph_eth(C::Array{Complex{Float64},2}, s::Int)
    ðC::Array{Complex{Float64},2}

Apply the differential operator ð ("eth") to the coefficients `C`.
This raises the spin weight `s` by 1. For a real function of spin
weight 0, this is equivalent to calculating the gradient.

This function is the converse of [`spinsph_ethbar`](@ref), which is a
derivative operator lowering the spin weight.

Use [`spinsph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode with spin weight `s`.

See also: [`spinsph_ethbar`](@ref), [`spinsph_mode`](@ref)
"""
function spinsph_eth(C::Array{Complex{Float64},2}, s::Int)
    # spinsph_eth!(copy(C), s)

    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ðC = zero(C)
    for l in max(abs(s + 1), abs(s)):(lmax + mmax), m in (-l):l
        if l - lmax ≤ abs(m) ≤ mmax
            sgn = m ≥ 0 ? 1 : -1
            ðC[spinsph_mode(s + 1, l, m)] = sgn *
                                            sqrt((l - s) * (l + s + 1)) *
                                            C[spinsph_mode(s, l, m)]
        end
    end

    return ðC
end

################################################################################

export spinsph_ethbar
"""
    ðC = spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)
    ðC::Array{Complex{Float64},2}

Apply the differential operator ð̄ ("eth-bar") to the coefficients `C`.
This lowers the spin weight `s` by 1. For a function of spin weight 1,
this is equivalent to calculating the divergence.

This function is the converse of [`spinsph_eth`](@ref), which is a
derivative operator raising the spin weight.

Use [`spinsph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode with spin weight `s`.

See also: [`spinsph_eth`](@ref), [`spinsph_mode`](@ref)
"""
function spinsph_ethbar(C::Array{Complex{Float64},2}, s::Int)
    # return spinsph_ethbar!(copy(C), s)

    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ð̄C = zero(C)
    for l in max(abs(s - 1), abs(s)):(lmax + mmax), m in (-l):l
        if l - lmax ≤ abs(m) ≤ mmax
            sgn = m ≥ 0 ? 1 : -1
            ð̄C[spinsph_mode(s - 1, l, m)] = -sgn *
                                             sqrt((l + s) * (l - s + 1)) *
                                             C[spinsph_mode(s, l, m)]
        end
    end

    return ð̄C
end
