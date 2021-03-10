module FastSphericalHarmonics

using FastTransforms
using LinearAlgebra

################################################################################

export SpHTypes
"""
    const SpHTypes = Union{Float64,Complex{Float64}}

The types supported by FastSphericalHarmonics (the same types as for
FastTransforms).
"""
const SpHTypes = Union{Float64,Complex{Float64}}

export sph_points
"""
    Θ, Φ = sph_points(N::Integer)
    Θ::Vector{Float64}
    Φ::Vector{Float64}

Calculate the locations of points on the sphere when using `N` points
in the θ (latitudinal) direction.

It is `length(Θ) = N` and `length(Φ) = M` where `M = 2N-1`.
"""
function sph_points(N::Int)
    @assert N > 0
    M = 2 * N - 1
    return π / N * (0.5:(N - 0.5)), 2π / M * (0:(M - 1))
end
sph_points(N::Integer) = sph_points(Int(N))

export sph_lmax
"""
    lmax = sph_lmax(N::Integer)
    lmax::Int

Calculate the maximum `l` mode that can be represented with `N`
points. It is `lmax = N - 1`.
"""
function sph_lmax(N::Int)
    @assert N > 0
    return N - 1
end
sph_lmax(N::Integer) = sph_lmax(Int(N))

################################################################################

export sph_mode
"""
    idx = sph_mode(l::Integer, m::Integer)
    idx::CartesianIndex{2}

Calculate the Cartesian index `idx` for the `l`,`m` mode. This index
can be used to access the coefficients, i.e. the result of
[`sph_transform`](@ref) or the input to [`sph_evaluate`](@ref).

Coefficients are stored in a two-dimensional array. Not all array
elements are used.

See also: [`sph_transform!`](@ref), [`sph_transform`](@ref),
[`sph_evaluate!`](@ref), [`sph_evaluate`](@ref)
"""
function sph_mode(l::Int, m::Int)
    @assert l ≥ 0
    @assert -l ≤ m ≤ l
    # (0,0) (1,-1) (1,1) (2,-2) (2,2)
    # (1,0) (2,-1) (2,1)
    # (2,0)
    return CartesianIndex(l - abs(m) + 1, 2 * abs(m) + (m ≥ 0))
end
sph_mode(l::Integer, m::Integer) = sph_mode(Int(l), Int(m))

export sph_transform!
"""
    sph_transform!(F::Array{T,2}) where {T<:SpHTypes}

Transform an array of points `F` into spherical harmonics. This is an
in-place transform, i.e. the array `F` will be overwritten by the
coefficients. Use [`sph_transform`](@ref) for a non-mutating function.

Use [`sph_points`](@ref) to caluclate the location of the points on
the sphere for the input array `F`.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`sph_transform`](@ref), [`sph_evaluate!`](@ref),
[`sph_points`](@ref), [`sph_mode`](@ref)
"""
function sph_transform!(F::Array{T,2}) where {T<:SpHTypes}
    N, M = size(F)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_sph2fourier(F)
    PA = plan_sph_analysis(F)
    C = F
    lmul!(PA, C)
    ldiv!(P, C)
    return C
end
export sph_transform
"""
    C = sph_transform(F::Array{T,2}) where {T<:SpHTypes}
    C::Array{T,2}

Transform an array of points `F` into spherical harmonics. You can use
[`sph_transform!`](@ref) for more efficient a mutating function that
overwrites its argument `F`.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere for the input array `F`.

Use [`sph_mode`](@ref) to calculate the location in the output
coefficient array for a particular `l`,`m` mode.

See also: [`sph_transform!`](@ref), [`sph_evaluate`](@ref),
[`sph_points`](@ref), [`sph_mode`](@ref)
"""
sph_transform(F::Array{T,2}) where {T<:SpHTypes} = sph_transform!(copy(F))

export sph_evaluate!
"""
    sph_evaluate!(C::Array{T,2}) where {T<:SpHTypes}

Evaluate an array of coefficients `C` on the grid points on a sphere.
This is an in-place transform, i.e. the array `C` will be overwritten
by the point values. Use [`sph_evaluate`](@ref) for a non-mutating
function.

Use [`sph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`sph_evaluate`](@ref), [`sph_transform!`](@ref),
[`sph_mode`](@ref), [`sph_points`](@ref)
"""
function sph_evaluate!(C::Array{T,2}) where {T<:SpHTypes}
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_sph2fourier(C)
    PS = plan_sph_synthesis(C)
    F = C
    lmul!(P, F)
    lmul!(PS, F)
    return F
end
export sph_evaluate
"""
    F = sph_evaluate(C::Array{T,2}) where {T<:SpHTypes}
    F::Array{T,2}

Evaluate an array of coefficients `C` on the grid points on a sphere.
You can use [`sph_evaluate!`](@ref) for more efficient a mutating
function that overwrites its argument `C`.

Use [`sph_mode`](@ref) to calculate the location in the input
coefficient array for a particular `l`,`m` mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`sph_evaluate!`](@ref), [`sph_transform`](@ref),
[`sph_mode`](@ref), [`sph_points`](@ref)
"""
sph_evaluate(C::Array{T,2}) where {T<:SpHTypes} = sph_evaluate!(copy(C))

################################################################################

export sphv_mode
"""
    idx = sphv_mode(l::Integer, m::Integer, v::Integer)
    idx::CartesianIndex{2}

Calculate the Cartesian index `idx` for the `v` components of `l`,`m`
mode. `v=1` is for the `θ` component, `v=2` for the `ϕ` component.
This index can be used to access the coefficients, i.e. the result of
[`sphv_transform`](@ref) or the input to [`sphv_evaluate`](@ref).

Coefficients are stored in a two-dimensional array. Not all array
elements are used.

See also: [`sphv_transform!`](@ref), [`sphv_transform`](@ref),
[`sphv_evaluate!`](@ref), [`sphv_evaluate`](@ref)
"""
function sphv_mode(l::Int, m::Int, v::Int)
    @assert l ≥ 1
    @assert -l ≤ m ≤ l
    @assert 1 ≤ v ≤ 2

    # e,l,m  -> lm

    # θ,1,0  -> 1,1

    # θ,1,-1 -> 2,2
    # ϕ,1,-1 -> 1,2
    # θ,1,1  -> 2,3
    # ϕ,1,1  -> 1,3

    # θ,2,0  -> 2,1

    # θ,2,-2 -> 2,4
    # ϕ,2,-2 -> 1,4
    # θ,2,-1 -> 3,2   (and 1,2)
    # ϕ,2,-1 -> 2,2
    # θ,2,1  -> 3,3   (and 1,3)
    # ϕ,2,1  -> 2,3
    # θ,2,2  -> 2,5
    # ϕ,2,2  -> 1,5

    # (0,0) (1,-1) (1,1) (2,-2) (2,2)
    # (1,0) (2,-1) (2,1)
    # (2,0)

    # e_θ:
    # (1,0) (*,**) (*,*) (*,**) (*,*)
    # (2,0) (1,-1) (1,1) (2,-2) (2,2)
    # (*,*) (2,-1) (2,1)

    # e_ϕ:
    # (*,*) (1,-1) (1,1) (2,-2), (2,2)
    # (*,*) (2,-1) (2,1)
    # (*,*)

    return m == 0 ? CartesianIndex(l, 1) :
           CartesianIndex(l - abs(m) + (3 - v), 2 * abs(m) + (m ≥ 0))
end
function sphv_mode(l::Integer, m::Integer, v::Integer)
    return sphv_mode(Int(l), Int(m), Int(v))
end

export sphv_transform!
"""
    sphv_transform!(F::Array{T,2}) where {T<:SpHTypes}

Transform an array of points `F` into vector spherical harmonics. The
`θ` and `ϕ` components of a vector field are transformed
independently, both by calling this function.

This is an in-place transform, i.e. the array `F` will be overwritten
by the coefficients. Use [`sphv_transform`](@ref) for a non-mutating
function.

Use [`sph_points`](@ref) to caluclate the location of the points on
the sphere for the input array `F`.

Use [`sphv_mode`](@ref) to calculate the location in the output
coefficient array for a particular component of a particular `l`,`m`
mode.

See also: [`sphv_transform`](@ref), [`sphv_evaluate!`](@ref),
[`sph_points`](@ref), [`sphv_mode`](@ref)
"""
function sphv_transform!(F::Array{T,2}) where {T<:SpHTypes}
    N, M = size(F)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_sphv2fourier(F)
    PA = plan_sphv_analysis(F)
    C = F
    lmul!(PA, C)
    ldiv!(P, C)
    return C
end
export sphv_transform
"""
    C = sphv_transform(F::Array{T,2}) where {T<:SpHTypes}
    C::Array{T,2}

Transform an array of points `F` into vector spherical harmonics. The
`θ` and `ϕ` components of a vector field are transformed
independently, both by calling this function.

You can use [`sphv_transform!`](@ref) for more efficient a mutating
function that overwrites its argument `F`.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere for the input array `F`.

Use [`sphv_mode`](@ref) to calculate the location in the output
coefficient array for a particular component of a particular `l`,`m`
mode.

See also: [`sphv_transform!`](@ref), [`sphv_evaluate`](@ref),
[`sph_points`](@ref), [`sphv_mode`](@ref)
"""
sphv_transform(F::Array{T,2}) where {T<:SpHTypes} = sphv_transform!(copy(F))

export sphv_evaluate!
"""
    sphv_evaluate!(C::Array{T,2}) where {T<:SpHTypes}

Evaluate an array of vector spherical harmonic coefficients `C` on the
grid points on a sphere. The `θ` and `ϕ` components of a vector field
are evaluated independently, both by calling this function for the
respective sets of coefficients.

This is an in-place transform, i.e. the array `C` will be overwritten
by the point values. Use [`sphv_evaluate`](@ref) for a non-mutating
function.

Use [`sphv_mode`](@ref) to calculate the location in the input
coefficient array for a particular component of a particular `l`,`m`
mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`sphv_evaluate`](@ref), [`sphv_transform!`](@ref),
[`sphv_mode`](@ref), [`sph_points`](@ref)
"""
function sphv_evaluate!(C::Array{T,2}) where {T<:SpHTypes}
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    P = plan_sphv2fourier(C)
    PS = plan_sphv_synthesis(C)
    F = C
    lmul!(P, F)
    lmul!(PS, F)
    return F
end
export sphv_evaluate
"""
    F = sphv_evaluate(C::Array{T,2}) where {T<:SpHTypes}
    F::Array{T,2}

Evaluate an array of vector spherical harmonic coefficients `C` on the
grid points on a sphere. The `θ` and `ϕ` components of a vector field
are evaluated independently, both by calling this function for the
respective sets of coefficients.

You can use [`sphv_evaluate!`](@ref) for more efficient a mutating
function that overwrites its argument `C`.

Use [`sphv_mode`](@ref) to calculate the location in the input
coefficient array for a particular component of a particular `l`,`m`
mode.

Use [`sph_points`](@ref) to caluclate the location of the points on the
sphere in the output array `F`.

See also: [`sphv_evaluate!`](@ref), [`sphv_transform`](@ref),
[`sphv_mode`](@ref), [`sph_points`](@ref)
"""
sphv_evaluate(C::Array{T,2}) where {T<:SpHTypes} = sphv_evaluate!(copy(C))

################################################################################

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

end
