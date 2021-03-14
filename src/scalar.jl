export sph_mode
"""
    idx = sph_mode(l::Integer, m::Integer)
    idx::CartesianIndex{2}

Calculate the Cartesian index `idx` for the `l`,`m` mode. This index
can be used to access the coefficients, i.e. the result of
[`sph_transform`](@ref) or the input to [`sph_evaluate`](@ref).

Coefficients are stored in a two-dimensional array. Not all array
elements are used. See [this
page](https://mikaelslevinsky.github.io/FastTransforms/transforms.html),
section "sph2fourier", for details.

See also: [`sph_transform!`](@ref), [`sph_transform`](@ref),
[`sph_evaluate!`](@ref), [`sph_evaluate`](@ref)
"""
function sph_mode(l::Int, m::Int)
    @assert l ≥ 0
    @assert abs(m) ≤ l
    # (0,0) (1,-1) (1,1) (2,-2) (2,2)
    # (1,0) (2,-1) (2,1)
    # (2,0)
    return CartesianIndex(l - abs(m) + 1, 2 * abs(m) + (m ≥ 0))
end
sph_mode(l::Integer, m::Integer) = sph_mode(Int(l), Int(m))

################################################################################

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

################################################################################

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

export sph_laplace!
"""
    sph_laplace!(C::Array{T,2}) where {T<:SpHTypes}

Calculate the Laplacian of a set of coefficients `C`. This is an
in-place transform, i.e. the array `C` will be overwritten by the
result. Use [`sph_laplace`](@ref) for a non-mutating function.

See also: [`sph_transform!`](@ref), [`sph_evaluate!`](@ref),
[`sph_laplace`](@ref)
"""
function sph_laplace!(C::Array{T,2}) where {T<:SpHTypes}
    N, M = size(C)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1
    mmax = M ÷ 2
    for l in 0:(lmax + mmax), m in (-l):l
        if l - lmax ≤ abs(m) ≤ mmax
            C[sph_mode(l, m)] *= -l * (l + 1)
        end
    end
    return C
end

export sph_laplace
"""
    ΔC = sph_laplace(C::Array{T,2}) where {T<:SpHTypes}
    ΔC::Array{T,2}

Calculate the Laplacian of a set of coefficients `C`. You can use
[`sph_laplace!`](@ref) for more efficient a mutating function that
overwrites its argument `C`.

See also: [`sph_transform`](@ref), [`sph_evaluate`](@ref),
[`sph_laplace!`](@ref)
"""
sph_laplace(C::Array{T,2}) where {T<:SpHTypes} = sph_laplace!(copy(C))
