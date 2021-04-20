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

function coeff_complex2real(C::AbstractArray{Complex{Float64},2}, s::Int)
    @assert s == 0
    C′ = Array{Float64}(undef, size(C))
    C′[:, 1] = real.(C[:, 1])
    for col in 2:2:size(C′, 2)
        for row in 1:size(C′, 1)
            avg = (C[row, col] + conj(C[row, col + 1])) / sqrt(2)
            C′[row, col] = imag(avg)
            C′[row, col + 1] = real(avg)
        end
    end
    return C′
end

function coeff_real2complex(C::AbstractArray{Float64,2}, s::Int)
    @assert s == 0
    C′ = Array{Complex{Float64}}(undef, size(C))
    C′[:, 1] = C[:, 1]
    for col in 2:2:size(C, 2)
        for row in 1:size(C, 1)
            val = Complex{Float64}(C[row, col + 1], C[row, col]) / sqrt(2)
            C′[row, col] = val
            C′[row, col + 1] = conj(val)
        end
    end
    return C′
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
    C = spinsph_transform(F::AbstractArray{Complex{Float64},2}, s::Int)
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
function spinsph_transform(F::AbstractArray{Complex{Float64},2}, s::Int)
    return spinsph_transform!(Array(F), s)
end
function spinsph_transform(F::AbstractArray{Float64,2}, s::Int)
    F′ = Array{Complex{Float64}}(F)
    C′ = spinsph_transform(F′, s)
    C = coeff_complex2real(C′, s)
    return C
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
    F = spinsph_evaluate(C::AbstractArray{Complex{Float64},2}, s::Int)
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
function spinsph_evaluate(C::AbstractArray{Complex{Float64},2}, s::Int)
    return spinsph_evaluate!(Array(C), s)
end
function spinsph_evaluate(C::AbstractArray{Float64,2}, s::Int)
    C′ = coeff_real2complex(C, s)
    F′ = spinsph_evaluate(C′, s)
    F = real.(F′)
    return F
end

################################################################################

export spinsph_eth
"""
    ðC = spinsph_eth(C::AbstractArray{Complex{Float64},2}, s::Int)
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
function spinsph_eth(C::AbstractArray{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ðC = zeros(Complex{Float64}, size(C))
    for l in max(abs(s + 1), abs(s)):(lmax + mmax), m in (-l):l
        if l - lmax ≤ abs(m) ≤ mmax
            ðC[spinsph_mode(s + 1, l, m)] = sqrt((l - s) * (l + s + 1)) *
                                            C[spinsph_mode(s, l, m)]
        end
    end

    return ðC
end

export spinsph_grad
function spinsph_grad(C::AbstractArray{Float64,2}, s::Int)
    C′ = coeff_real2complex(C, s)
    ðC = spinsph_eth(C′, s)
    return ðC
end

export spinsph_divbar
function spinsph_divbar(C::AbstractArray{Complex{Float64},2}, s::Int)
    ð̄̄C′ = spinsph_eth(C, s)
    ð̄̄C = coeff_complex2real(ð̄̄C′, s + 1)
    return ð̄̄C
end

################################################################################

export spinsph_ethbar
"""
    ðC = spinsph_ethbar(C::AbstractArray{Complex{Float64},2}, s::Int)
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
function spinsph_ethbar(C::AbstractArray{Complex{Float64},2}, s::Int)
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2

    ð̄C = zeros(Complex{Float64}, size(C))
    for l in max(abs(s - 1), abs(s)):(lmax + mmax), m in (-l):l
        if l - lmax ≤ abs(m) ≤ mmax
            ð̄C[spinsph_mode(s - 1, l, m)] = -sqrt((l + s) * (l - s + 1)) *
                                             C[spinsph_mode(s, l, m)]
        end
    end

    return ð̄C
end

export spinsph_gradbar
function spinsph_gradbar(C::AbstractArray{Float64,2}, s::Int)
    C′ = coeff_real2complex(C, s)
    ð̄C = spinsph_ethbar(C′, s)
    return ð̄C
end

export spinsph_div
function spinsph_div(C::AbstractArray{Complex{Float64},2}, s::Int)
    ð̄C′ = spinsph_ethbar(C, s)
    ð̄C = coeff_complex2real(ð̄C′, s - 1)
    return ð̄C
end
