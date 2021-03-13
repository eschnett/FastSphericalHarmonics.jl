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
