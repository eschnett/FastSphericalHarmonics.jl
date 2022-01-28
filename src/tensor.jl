stensor(D::Integer) = SArray{Tuple{ntuple(d -> 2, D)...}}
stensor(D::Integer, T::Type) = stensor(D){T}

export Tensor
@computed struct Tensor{D,T} # <: AbstractArray{D,T}
    values::AbstractArray{stensor(D, T),2}
end

@doc """
    struct Tensor{D,T}

Tensor of rank `D` with elements of type `T`, stored component-wise.

The components are stored "normalized", i.e. projected onto the dyad
vectors `eθ^a` and `eϕ^a`. This removes singularities at the poles by
multiplying the `ϕ` components with respective powers of `sin θ`. This
also means that covariant (index down) and contravariant (index up)
indices are represented in the same way.

eθ^a = [1, 0]           eθ_a = [1, 0]
eϕ^a = [0, 1/(sin θ)]   eϕ_a = [1, sin θ]

eθ and eϕ are orthogonal, and both have length 1.

g_ab = diag[1, (sin θ)^2]
g_ab = eθ_a eϕ_b + eθ_b eϕ_a = m_a m̄_b + m̄_a m_b

Example:

Scalar `s`: stored as is
Vector `v^a`: store `[eθ_a v^a, eϕ_a v_a]`
Tensor `t_ab`: store `t[1,1] = eθ^a eθ^b t_ab`
                     `t[1,2] = eθ^a eϕ^b t_ab`
                     `t[2,1] = eϕ^a eθ^b t_ab`
                     `t[2,2] = eϕ^a eϕ^b t_ab`
""" Tensor

export SpinTensor
@computed struct SpinTensor{D,T} # <: AbstractArray{D,T}
    # We really want `<:AbstractArray{T,2}` instead
    coeffs::stensor(D, AbstractArray{T,2})
end

@doc """
    struct SpinTensor{D,T}

Tensor of rank `D` with elements of type `T`, stored separated by spin
weights.

Grouping tensor components by their spin weight allows representing
tensors via spin-weighted spherical harmonics, and allows calculating
covariant derivatives (covariant with respect to the unit sphere).
""" SpinTensor

# Convenience constructors
function Tensor{D}(values::AbstractArray{<:SArray{X,T} where {X},2}) where {D,T}
    return Tensor{D,T}(values)
end
Tensor{0}(values::AbstractArray{<:Number,2}) = Tensor{0}(stensor(0).(values))

function SpinTensor{D}(coeffs::SArray{X,<:AbstractArray{T,2}} where {X}) where {D,
                                                                                T}
    return SpinTensor{D,T}(coeffs)
end
function SpinTensor{0}(coeffs::AbstractArray{<:Number,2})
    return SpinTensor{0}(stensor(0)(coeffs))
end

# Basic operations
Base.eltype(::Tensor{D,T}) where {D,T} = T
Base.:(==)(t1::Tensor{D}, t2::Tensor{D}) where {D} = t1.values == t2.values
function Base.isapprox(t1::Tensor{D}, t2::Tensor{D}; kws...) where {D}
    return isapprox(t1.values, t2.values; kws...)
end
Base.zero(t::Tensor{D}) where {D} = Tensor{D}(zero(t.values))
Base.:-(t::Tensor{D}) where {D} = Tensor{D}(-t.values)
Base.conj(t::Tensor{D}) where {D} = Tensor{D}(conj(t.values))
function Base.:+(t1::Tensor{D}, t2::Tensor{D}) where {D}
    return Tensor{D}(t1.values + t2.values)
end
function Base.:-(t1::Tensor{D}, t2::Tensor{D}) where {D}
    return Tensor{D}(t1.values - t2.values)
end
Base.:*(a::Number, t::Tensor{D}) where {D} = Tensor{D}(a * t.values)
Base.:*(t::Tensor{D}, a::Number) where {D} = Tensor{D}(t.values * a)
Base.:/(t::Tensor{D}, a::Number) where {D} = Tensor{D}(t.values / a)

Base.eltype(::SpinTensor{D,T}) where {D,T} = T
function Base.:(==)(t1::SpinTensor{D}, t2::SpinTensor{D}) where {D}
    return t1.coeffs == t2.coeffs
end
function Base.isapprox(t1::SpinTensor{D}, t2::SpinTensor{D}; kws...) where {D}
    return isapprox(t1.coeffs, t2.coeffs; kws...)
end
Base.zero(t::SpinTensor{D}) where {D} = SpinTensor{D}(zero(t.coeffs))
Base.:-(t::SpinTensor{D}) where {D} = SpinTensor{D}(-t.coefffs)
Base.conj(t::SpinTensor{D}) where {D} = SpinTensor{D}(conj(t.coeffs))
function Base.:+(t1::SpinTensor{D}, t2::SpinTensor{D}) where {D}
    return SpinTensor{D}(t1.coeffs + t2.coeffs)
end
function Base.:-(t1::SpinTensor{D}, t2::SpinTensor{D}) where {D}
    return SpinTensor{D}(t1.coeffs - t2.coeffs)
end
Base.:*(a::Number, t::SpinTensor{D}) where {D} = SpinTensor{D}(a * t.coeffs)
Base.:*(t::SpinTensor{D}, a::Number) where {D} = SpinTensor{D}(t.coeffs * a)
Base.:/(t::SpinTensor{D}, a::Number) where {D} = SpinTensor{D}(t.coeffs / a)

"Convert `Tensor` to `SpinTensor`"
function SpinTensor{D}(tensor::Tensor{D}) where {D}
    T = eltype(tensor)
    @assert T <: Number
    CT = typeof(Complex(zero(T)))
    tensor::Tensor{D,T}
    # This represents `m` in terms of `eθ` and `eϕ`
    #     m^a = eθ^a + im eϕ^a
    #     m^a m_a = 0
    #     m^a m̄_a = 2
    m = SVector{2}(1, im)
    m̄ = conj(m)
    if D == 0
        # Avoid real-valued spin spherical harmonics
        values = [Complex(v[]) for v in tensor.values]
        coeffs = spinsph_transform(values, 0)
        return SpinTensor{D}(Scalar(coeffs))::SpinTensor{D,CT}
    end
    if D == 1
        values_m = [sum(m[a] * v[a] for a in 1:2) for v in tensor.values]::Array{CT,
                                                                                 2}
        values_m̄ = [sum(m̄[a] * v[a] for a in 1:2) for v in tensor.values]::Array{CT,
                                                                                   2}
        coeffs_m = spinsph_transform(values_m, +1)
        coeffs_m̄ = spinsph_transform(values_m̄, -1)
        return SpinTensor{D}(SVector{2}(coeffs_m, coeffs_m̄))::SpinTensor{D,CT}
    end
    if D == 2
        values_mm = [sum(m[a] * m[b] * v[a, b] for a in 1:2, b in 1:2)
                     for v in tensor.values]
        values_mm̄ = [sum(m[a] * m̄[b] * v[a, b] for a in 1:2, b in 1:2)
                      for v in tensor.values]
        values_m̄m = [sum(m̄[a] * m[b] * v[a, b] for a in 1:2, b in 1:2)
                      for v in tensor.values]
        values_m̄m̄ = [sum(m̄[a] * m̄[b] * v[a, b] for a in 1:2, b in 1:2)
                       for v in tensor.values]
        coeffs_mm = spinsph_transform(values_mm, +2)
        coeffs_mm̄ = spinsph_transform(values_mm̄, 0)
        coeffs_m̄m = spinsph_transform(values_m̄m, 0)
        coeffs_m̄m̄ = spinsph_transform(values_m̄m̄, -2)
        return SpinTensor{D}(SMatrix{2,2}(coeffs_mm, coeffs_m̄m, coeffs_mm̄,
                                          coeffs_m̄m̄))::SpinTensor{D,CT}
    end
    @assert false
end
SpinTensor(tensor::Tensor{D}) where {D} = SpinTensor{D}(tensor)

"Convert `SpinTensor` to `Tensor`"
function Tensor{D}(spintensor::SpinTensor{D}) where {D}
    T = eltype(spintensor)
    @assert T <: Complex
    spintensor::SpinTensor{D,T}
    # See above
    m = SVector{2}(1, im)
    m̄ = conj(m)
    if D == 0
        coeffs = spintensor.coeffs[]
        values = spinsph_evaluate(coeffs, 0)
        return Tensor{D}(Scalar.(values))::Tensor{D,T}
    end
    if D == 1
        values_m = spinsph_evaluate(spintensor.coeffs[1], +1)
        values_m̄ = spinsph_evaluate(spintensor.coeffs[2], -1)
        values = [SVector{2}((vm * m̄[a] + vm̄ * m[a]) / 2 for a in 1:2)
                  for (vm, vm̄) in zip(values_m, values_m̄)]
        return Tensor{D}(values)::Tensor{D,T}
    end
    if D == 2
        values_mm = spinsph_evaluate(spintensor.coeffs[1, 1], +2)
        values_mm̄ = spinsph_evaluate(spintensor.coeffs[1, 2], 0)
        values_m̄m = spinsph_evaluate(spintensor.coeffs[2, 1], 0)
        values_m̄m̄ = spinsph_evaluate(spintensor.coeffs[2, 2], -2)
        values = [SMatrix{2,2}((vmm * m̄[a] * m̄[b] +
                                vmm̄ * m̄[a] * m[b] +
                                vm̄m * m[a] * m̄[b] +
                                vm̄m̄ * m[a] * m[b]) / 4 for a in 1:2, b in 1:2)
                  for (vmm, vmm̄, vm̄m, vm̄m̄) in
                      zip(values_mm, values_mm̄, values_m̄m, values_m̄m̄)]
        return Tensor{D}(values)::Tensor{D,T}
    end
    @assert false
end
Tensor(spintensor::SpinTensor{D}) where {D} = Tensor{D}(spintensor)

export tensor_gradient
"Calculate gradient"
function tensor_gradient(spintensor::SpinTensor{D}) where {D}
    T = eltype(spintensor)
    @assert T <: Complex
    spintensor::SpinTensor{D,T}
    if D == 0
        coeffs = spintensor.coeffs[]
        dcoeffs_m = spinsph_eth(coeffs, 0)
        dcoeffs_m̄ = spinsph_ethbar(coeffs, 0)
        return SpinTensor{D + 1}(stensor(D + 1)(dcoeffs_m, dcoeffs_m̄))::SpinTensor{D +
                                                                                    1,
                                                                                    T}
    end
    if D == 1
        coeffs_m = spintensor.coeffs[1]
        coeffs_m̄ = spintensor.coeffs[2]
        dcoeffs_mm = spinsph_eth(coeffs_m, 1)
        dcoeffs_mm̄ = spinsph_ethbar(coeffs_m, 1)
        dcoeffs_m̄m = spinsph_eth(coeffs_m̄, -1)
        dcoeffs_m̄m̄ = spinsph_ethbar(coeffs_m̄, -1)
        return SpinTensor{D + 1}(stensor(D + 1)(dcoeffs_mm, dcoeffs_m̄m,
                                                dcoeffs_mm̄, dcoeffs_m̄m̄))::SpinTensor{D +
                                                                                        1,
                                                                                        T}
    end
    @assert false
end
