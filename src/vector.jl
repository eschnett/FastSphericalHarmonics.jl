function mul2(A::Array{T,4} where {T}, x::Array{T,2} where {T})
    m = size(A)[1:2]
    n = size(A)[3:4]
    @assert size(x) == n
    y = reshape(A, prod(m), prod(n)) * reshape(x, prod(n))
    return reshape(y, m)
end

function ldiv2(A::Array{T,4} where {T}, x::Array{T,2} where {T})
    m = size(A)[1:2]
    n = size(A)[3:4]
    @assert size(x) == m
    y = reshape(A, prod(m), prod(n)) \ reshape(x, prod(m))
    r = reshape(y, n)
    if !(mul2(A, r) ≈ x)
        @show A r mul2(A, r) x
    end
    @assert mul2(A, r) ≈ x
    return reshape(y, n)
end

################################################################################

# Julien Molina, Richard Mikael Slevinsky, "A rapid and
# well-conditioned algorithm for the Helmholtz--Hodge decomposition of
# vector fields on the sphere", https://arxiv.org/abs/1809.04555

function sphv_mode(lmax::Int, l::Int, m::Int)
    # eq. (11)
    @assert l ≥ abs(abs(m) - 1)
    @assert l ≤ lmax
    if m == 0
        @assert 1 ≤ l ≤ lmax
        return CartesianIndex(l, 1)
    end
    @assert abs(m) ≤ lmax
    if !(1 ≤ l - abs(m) + 1 ≤ lmax + 1)
        @show lmax l m
    end
    @assert 1 ≤ l - abs(m) + 1 ≤ lmax + 1
    @assert 1 ≤ 2 * abs(m) + (m ≥ 0) ≤ 2lmax + 1
    return CartesianIndex(l - abs(m) + 1, 2 * abs(m) + (m ≥ 0))
end

export sphv_transform
function sphv_transform(Fθ::Array{T,2}, Fϕ::Array{T,2}) where {T<:SpHTypes}
    @assert size(Fθ) == size(Fϕ)
    N, M = size(Fθ)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    lmax = N - 1

    # Convert from position basis to Z basis
    P = plan_sphv2fourier(Fθ)
    PA = plan_sphv_analysis(Fθ)
    Zθ = P \ (PA * Fθ)
    Zϕ = P \ (PA * Fϕ)
    @show Zθ Zϕ

    # Convert from Z basis to cscθY basis
    # Z_lm = α_lm cscθY_(l−1)m + β_lm cscθY_(l+1)m
    # m ∈ -∞:∞, l ∈ ||m|-1|:∞
    ZfromV = zeros(N, M, N, M)
    for m in (-lmax):lmax, l in abs(abs(m) - 1):lmax
        @show :ZfromV l m
        if m == 0
            # eq. (9)
            # Pl1(cos(θ)) = ∂θ Pl0(cos(θ))
            ZfromV[sphv_mode(lmax, l, m), sphv_mode(lmax, l, m)] = 1
        else
            # eq. (15)
            α(l, m) = -sqrt((l - m) * (l - m + 1) / ((2l - 1) * (2l + 1)))
            β(l, m) = -sqrt((l + m) * (l + m + 1) / ((2l + 1) * (2l + 3)))
            # eq. (14)
            if abs(m) ≤ l - 1 && l - 1 ≥ 0
                @show α(l, m)
                ZfromV[sphv_mode(lmax, l, m), sphv_mode(lmax, l - 1, m)] = α(l,
                                                                             m)
            end
            if abs(m) ≤ l && l + 1 ≤ lmax
                @show β(l, m)
                ZfromV[sphv_mode(lmax, l, m), sphv_mode(lmax, l + 1, m)] = β(l,
                                                                             m)
            end
        end
    end
    @show ZfromV
    Vθ = ldiv2(ZfromV, Zθ)
    Vϕ = ldiv2(ZfromV, Zϕ)
    @show Vθ Vϕ

    # Convert from cscθY basis to ∇Y basis
    # ∇θY_lm = ∂θ Y_lm = γ_lm csc(θ) Y_(l-1)m + δ_lm Y_(l+1)m
    # l ∈ 0:lmax, m ∈ -l:l
    CθfromVθ = zeros(N, M, N, M)
    CϕfromVϕ = zeros(N, M, N, M)
    for l in 1:lmax, m in (-l):l
        if m == 0
            CθfromVθ[sphv_mode(lmax, l, m), sphv_mode(lmax, l, m)] = 1
            CϕfromVϕ[sphv_mode(lmax, l, m), sphv_mode(lmax, l, m)] = 1
        else
            # eq. (18)
            γ(l, m) = -(l + 1) * sqrt((l - m) * (l + m) / ((2l - 1) * (2l + 1)))
            function δ(l, m)
                return l *
                       sqrt((l - m + 1) * (l + m + 1) / ((2l + 1) * (2l + 3)))
            end
            # eq. (17)
            if l - 1 ≥ 0
                CθfromVθ[sphv_mode(lmax, l, m), sphv_mode(lmax, l - 1, m)] = γ(l,
                                                                               m)
            end
            if l + 1 ≤ lmax
                CθfromVθ[sphv_mode(lmax, l, m), sphv_mode(lmax, l + 1, m)] = δ(l,
                                                                               m)
            end
            # eq. (22)
            CϕfromVϕ[sphv_mode(lmax, l, m), sphv_mode(lmax, l, -m)] = -m
        end
    end
    Cθ = mul2(CθfromVθ, Vθ)
    Cϕ = mul2(CϕfromVϕ, Vϕ)
    @show Cθ Cϕ

    return Cθ, Cϕ
end

export sphv_evaluate
function sphv_evaluate(Cθ::Array{T,2}, Cϕ::Array{T,2}) where {T<:SpHTypes}
    @assert size(Cθ) == size(Cϕ)
    N, M = size(Cθ)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1

    @assert false
end
