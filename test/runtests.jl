using FastTransforms
using FastSphericalHarmonics
using SpecialFunctions
using Test

bitsign(b::Bool) = b ? -1 : 1
bitsign(n::Integer) = bitsign(isodd(n))

chop(x) = abs2(x) < 100eps(x) ? zero(x) : x
chop(x::Complex) = Complex(chop(real(x)), chop(imag(x)))

# [Generalized binomial coefficient](https://en.wikipedia.org/wiki/Binomial_coefficient#Generalization_and_connection_to_the_binomial_series)
function Base.binomial(α::Number, k::Integer)
    k == 0 && return one(α)
    return prod((α - i) / (k - i) for i in 0:(k - 1))
end

# [Jacobi
# Polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials)
function JacobiP(α, β, n, z)
    return gamma(α + n + 1) / (factorial(n) * gamma(α + β + n + 1)) *
           sum(binomial(n, m) * gamma(α + β + n + m + 1) / gamma(α + m + 1) *
               ((z - 1) / 2)^m for m in 0:n)
end

# [Associated Legendre
# Polynomials](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials)
function LegendreP(l, m, x)
    return bitsign(m) *
           2^l *
           sqrt(1 - x^2)^m *
           sum(factorial(k) / factorial(k - m) *
               x^(k - m) *
               binomial(l, k) *
               binomial((l + k - 1) / 2, l) for k in m:l)
end

# [Spherical
# harmonics](https://mikaelslevinsky.github.io/FastTransforms/transforms.html)
# (section "sph2fourier")
function Ylm(l, m, θ, ϕ)
    return bitsign(abs(m)) *
           sqrt((l + 1 / 2) * factorial(l - abs(m)) / factorial(l + abs(m))) *
           LegendreP(l, abs(m), cos(θ)) *
           sqrt((2 - (m == 0)) / 2π) *
           (m ≥ 0 ? cos(abs(m) * ϕ) : sin(abs(m) * ϕ))
end

# [Spin-weighted spherical
# harmonics](https://mikaelslevinsky.github.io/FastTransforms/transforms.html)
# (section "spinsph2fourier")
function sYlm(s, l, m, θ, ϕ)
    l0 = max(abs(m), abs(s))
    l1 = min(abs(m), abs(s))
    return cis(m * ϕ) / sqrt(2π) *
           sqrt((l + 1 / 2) * factorial(l + l0) * factorial(l - l0) /
                (factorial(l + l1) * factorial(l - l1))) *
           sin(θ / 2)^abs(m + s) *
           cos(θ / 2)^abs(m - s) *
           JacobiP(abs(m + s), abs(m - s), l - l0, cos(θ))
end
sYlm(::Type{<:Complex}, s, l, m, θ, ϕ) = sYlm(s, l, m, θ, ϕ)
function sYlm(::Type{<:Real}, s, l, m, θ, ϕ)
    @assert s == 0
    if m == 0
        return real(sYlm(s, l, abs(m), θ, ϕ))
    elseif m > 0
        return sqrt(2) * real(sYlm(s, l, abs(m), θ, ϕ))
    else
        return sqrt(2) * imag(sYlm(s, l, abs(m), θ, ϕ))
    end
end

# Real spherical harmonics
Y_0_0(θ, ϕ) = sqrt(1 / 4π)

Y_1_0(θ, ϕ) = sqrt(3 / 4π) * cos(θ)
Y_1p1(θ, ϕ) = sqrt(3 / 4π) * sin(θ) * cos(ϕ)
Y_1m1(θ, ϕ) = sqrt(3 / 4π) * sin(θ) * sin(ϕ)

# Complex spin-weighted spherical harmonics

# sYlm[s_, l_, m_, q_, f_] = 
#  Module[{l0 = Max[Abs[m], Abs[s]], l1 = Min[Abs[m], Abs[s]]}, 
#      (Exp[I*m*f]/Sqrt[2*Pi])*
#    Sqrt[(l + 1/2)*(l + l0)!*((l - l0)!/((l + l1)!*(l - l1)!))]*
#        Sin[q/2]^Abs[m + s]*Cos[q/2]^Abs[m - s]*
#    JacobiP[l - l0, Abs[m + s], Abs[m - s], 
#          Cos[q]]]

Y_s0_0_0(θ, ϕ) = 1 / (2 * sqrt(π))

Y_s0_1_0(θ, ϕ) = 1 / 2 * sqrt(3 / π) * cos(θ)
Y_s0_1p1(θ, ϕ) = cis(ϕ) * sqrt(3 / (2 * π)) * cos(θ / 2) * sin(θ / 2)
Y_s0_1p1(θ, ϕ) = cis(-f) * sqrt(3 / (2 * π)) * cos(θ / 2) * sin(θ / 2)

Y_s1_1_0(θ, ϕ) = sqrt(3 / (2 * π)) * cos(θ / 2) * sin(θ / 2)
Y_s1_1p1(θ, ϕ) = 1 / 2 * cis(f) * sqrt(3 / π) * sin(θ / 2)^2
Y_s1_1m1(θ, ϕ) = 1 / 2 * cis(-f) * sqrt(3 / π) * cos(θ / 2)^2

Y_s1_2_0(θ, ϕ) = sqrt(15 / (2 * π)) * cos(θ / 2) * cos(θ) * sin(θ / 2)
Y_s1_2p1(θ, ϕ) = 1 / 4 * cis(f) * sqrt(5 / π) * (2 + 4 * cos(θ)) * sin(θ / 2)^2
Y_s1_2p2(θ, ϕ) = cis(2 * f) * sqrt(5 / π) * cos(θ / 2) * sin(θ / 2)^3

# # Real gradient spherical harmonics
# # [∂θ, 1/sin(θ) ∂ϕ]
# grad_Y_1_0_θ(θ, ϕ) = sqrt(3 / 4π) * sin(θ)
# grad_Y_1_0_ϕ(θ, ϕ) = 0.0
# grad_Y_1p1_θ(θ, ϕ) = sqrt(3 / 4π) * cos(θ) * cos(ϕ)
# grad_Y_1p1_ϕ(θ, ϕ) = sqrt(3 / 4π) * cos(ϕ)
# grad_Y_1m1_θ(θ, ϕ) = sqrt(3 / 4π) * cos(θ) * sin(ϕ)
# grad_Y_1m1_ϕ(θ, ϕ) = sqrt(3 / 4π) * sin(ϕ)
# 
# grad_Y_2_0_θ(θ, ϕ) = sqrt(45 / 4π) * cos(θ) * sin(θ)
# grad_Y_2_0_ϕ(θ, ϕ) = 0.0
# grad_Y_2p1_θ(θ, ϕ) = sqrt(15 / 4π) * cos(2θ) * cos(ϕ)
# grad_Y_2p1_ϕ(θ, ϕ) = sqrt(15 / 4π) * cos(θ) * cos(ϕ)
# 
# grad_Y_3_0_θ(θ, ϕ) = sqrt(63 / 256π) * (5 * sin(3θ) + sin(θ))
# grad_Y_3_0_ϕ(θ, ϕ) = 0.0
# grad_Y_3p1_θ(θ, ϕ) = sqrt(42 / 256π) * (15 * cos(2θ) - 7) * cos(θ) * cos(ϕ)
# grad_Y_3p1_ϕ(θ, ϕ) = sqrt(42 / 256π) * (5 * cos(2θ) + 3) * cos(ϕ)
# 
# # Real curl spherical harmonics
# # [1/sin(θ) ∂ϕ, -∂θ]
# curl_Y_1_0_θ(θ, ϕ) = grad_Y_1_0_ϕ(θ, ϕ)
# curl_Y_1_0_ϕ(θ, ϕ) = -grad_Y_1_0_θ(θ, ϕ)
# curl_Y_1p1_θ(θ, ϕ) = grad_Y_1p1_ϕ(θ, ϕ)
# curl_Y_1p1_ϕ(θ, ϕ) = -grad_Y_1p1_θ(θ, ϕ)
# curl_Y_1m1_θ(θ, ϕ) = grad_Y_1m1_ϕ(θ, ϕ)
# curl_Y_1m1_ϕ(θ, ϕ) = -grad_Y_1m1_θ(θ, ϕ)

# Integrate over the sphere
function integrate(u::Array{T1,2}, v::Array{T2,2}) where {T1,T2}
    @assert size(u) == size(v)
    N, M = size(u)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    Θ, Φ = sph_points(N)
    s = zero(conj(one(T1)) * one(T2) * one(eltype(Θ)))
    for j in 1:M, i in 1:N
        s += conj(u[i, j]) * v[i, j] * sin(Θ[i])
    end
    return 2π^2 / (N * M) * s
end

function integrate(uθ::Array{T1,2}, uϕ::Array{T1,2}, vθ::Array{T2,2},
                   vϕ::Array{T2,2}) where {T1,T2}
    @assert size(uθ) == size(uϕ) == size(vθ) == size(vϕ)
    N, M = size(uθ)
    @assert M > 0 && N > 0
    @assert M == 2 * N - 1
    Θ, Φ = sph_points(N)
    s = zero(conj(one(T1)) * one(T2) * one(eltype(Θ)))
    for j in 1:M, i in 1:N
        s += (conj(uθ[i, j]) * vθ[i, j] + conj(uϕ[i, j]) * vϕ[i, j]) * sin(Θ[i])
    end
    return 2π^2 / (N * M) * s
end

################################################################################

include("test_fasttransforms.jl")
include("test_points.jl")
include("test_scalar.jl")
# include("test_vector.jl")
include("test_spin.jl")
