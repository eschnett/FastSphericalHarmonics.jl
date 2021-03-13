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
