using FastTransforms
using FastSphericalHarmonics
using Test

chop(x) = abs2(x) < 10eps(x) ? zero(x) : x

# Real Spherical Harmonics
Y_0_0(θ, ϕ) = sqrt(1 / 4π)

Y_1_0(θ, ϕ) = sqrt(3 / 4π) * cos(θ)
Y_1p1(θ, ϕ) = sqrt(3 / 4π) * sin(θ) * cos(ϕ)
Y_1m1(θ, ϕ) = sqrt(3 / 4π) * sin(θ) * sin(ϕ)

Y_2_0(θ, ϕ) = sqrt(5 / 64π) * (3 * cos(2θ) - 1)
Y_2p1(θ, ϕ) = sqrt(15 / 16π) * sin(2θ) * cos(ϕ)
Y_2p2(θ, ϕ) = sqrt(15 / 64π) * (-cos(2θ) + 1) * cos(2ϕ)

Y_3_0(θ, ϕ) = sqrt(7 / 16π) * (5 * cos(θ)^3 - 3 * cos(θ))
Y_3p1(θ, ϕ) = sqrt(21 / 64π) * (5 * cos(θ)^2 - 1) * sin(θ)

# Real Gradient Spherical Harmonics
# [∂θ, 1/sin(θ) ∂ϕ]
grad_Y_1_0_θ(θ, ϕ) = sqrt(3 / 4π) * sin(θ)
grad_Y_1_0_ϕ(θ, ϕ) = 0.0
grad_Y_1p1_θ(θ, ϕ) = sqrt(3 / 4π) * cos(θ) * cos(ϕ)
grad_Y_1p1_ϕ(θ, ϕ) = sqrt(3 / 4π) * cos(ϕ)
grad_Y_1m1_θ(θ, ϕ) = sqrt(3 / 4π) * cos(θ) * sin(ϕ)
grad_Y_1m1_ϕ(θ, ϕ) = sqrt(3 / 4π) * sin(ϕ)

grad_Y_2_0_θ(θ, ϕ) = sqrt(45 / 4π) * cos(θ) * sin(θ)
grad_Y_2_0_ϕ(θ, ϕ) = 0.0
grad_Y_2p1_θ(θ, ϕ) = sqrt(15 / 4π) * cos(2θ) * cos(ϕ)
grad_Y_2p1_ϕ(θ, ϕ) = sqrt(15 / 4π) * cos(θ) * cos(ϕ)

grad_Y_3_0_θ(θ, ϕ) = sqrt(63 / 256π) * (5 * sin(3θ) + sin(θ))
grad_Y_3_0_ϕ(θ, ϕ) = 0.0
grad_Y_3p1_θ(θ, ϕ) = sqrt(42 / 256π) * (15 * cos(2θ) - 7) * cos(θ) * cos(ϕ)
grad_Y_3p1_ϕ(θ, ϕ) = sqrt(42 / 256π) * (5 * cos(2θ) + 3) * cos(ϕ)

# Real Curl Spherical Harmonics
# [1/sin(θ) ∂ϕ, -∂θ]
curl_Y_1_0_θ(θ, ϕ) = grad_Y_1_0_ϕ(θ, ϕ)
curl_Y_1_0_ϕ(θ, ϕ) = -grad_Y_1_0_θ(θ, ϕ)
curl_Y_1p1_θ(θ, ϕ) = grad_Y_1p1_ϕ(θ, ϕ)
curl_Y_1p1_ϕ(θ, ϕ) = -grad_Y_1p1_θ(θ, ϕ)
curl_Y_1m1_θ(θ, ϕ) = grad_Y_1m1_ϕ(θ, ϕ)
curl_Y_1m1_ϕ(θ, ϕ) = -grad_Y_1m1_θ(θ, ϕ)

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

include("test_points.jl")
include("test_scalar.jl")
# include("test_vector.jl")
include("test_spin.jl")
