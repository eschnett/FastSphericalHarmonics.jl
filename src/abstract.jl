export ash_grid_size, ash_nmodes
function ash_grid_size(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = Int(lmax) + 1
    M = 2 * N - 1
    return (N, M)::NTuple{2,Int}
end
function ash_nmodes(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = Int(lmax) + 1
    M = 2 * N - 1
    return (N, M)::NTuple{2,Int}
end

export ash_ntheta, ash_nphi, ash_thetas, ash_phis, ash_point_coord,
       ash_point_delta, ash_grid_as_phi_theta
ash_ntheta(lmax) = ash_grid_size(lmax)[1]::Int
ash_nphi(lmax) = ash_grid_size(lmax)[2]::Int
function ash_thetas(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = Int(lmax) + 1
    return sph_points(N)[1]
end
function ash_phis(lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax))
    N = Int(lmax) + 1
    return sph_points(N)[2]
end
function ash_point_coord(ij::Union{CartesianIndex{2},NTuple{2,Int}},
                         lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = Int(lmax) + 1
    M = 2 * N - 1
    t, p = Tuple(ij)
    theta = (π / N * (0.5:(N - 0.5)))[t]
    phi = (2π / M * (0:(M - 1)))[p]
    return theta, phi
end
function ash_point_delta(ij::Union{CartesianIndex{2},NTuple{2,Int}},
                         lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = Int(lmax) + 1
    M = 2 * N - 1
    dtheta = π / N
    dphi = 2π / M
    return dtheta, dphi
end
ash_grid_as_phi_theta(grid::AbstractMatrix) = transpose(grid)

export ash_mode_index
function ash_mode_index(s::Integer, l::Integer, m::Integer, lmax::Integer)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    abs(s) ≤ l ≤ lmax || throw(DomainError(l, "Need abs(s) ≤ l ≤ lmax"))
    -l ≤ m ≤ l || throw(DomainError(m, "Need -l ≤ m ≤ l"))
    return spinsph_mode(s, l, m)::CartesianIndex{2}
end
export ash_mode_numbers
function ash_mode_numbers(s::Int, ind::CartesianIndex{2}, lmax::Int)
    0 ≤ lmax || throw(DomainError(lmax, "Need 0 ≤ lmax"))
    N = lmax + 1
    M = 2 * N - 1
    # ind[2] = 2 * abs(m) + (m ≥ 0)
    msign = isodd(ind[2]) ? +1 : -1
    mabs = ind[2] ÷ 2
    m = msign * mabs
    # ind[1] = l - max(abs(s), abs(m)) + 1
    l = ind[1] + max(abs(s), abs(m)) - 1
    @assert ash_mode_index(s, l, m, lmax) == ind
    return (l, m)::NTuple{2,Int}
end
function ash_mode_numbers(s::Integer, ind::CartesianIndex{2}, lmax::Integer)
    return ash_mode_numbers(Int(s), ind, Int(lmax))
end
function ash_mode_numbers(s::Integer, ind::NTuple{2,<:Integer}, lmax::Integer)
    return ash_mode_numbers(Int(s), CartesianIndex{2}(ind...), Int(lmax))
end

function change_signs!(flm::AbstractArray{<:Complex}, s::Integer)
    for ind in CartesianIndices(flm)
        row, col = Tuple(ind)
        m = col == 1 ? 0 : (col % 2 == 0 ? -1 : 1) * (col ÷ 2)
        l = row + max(abs(s), abs(m)) - 1
        sign = (isodd(s) && m < -s) || (m >= -s && isodd(m)) ? -1 : +1
        flm[ind] = flipsign(flm[ind], sign)
    end
    return flm
end

export ash_transform!, ash_transform, ash_evaluate!, ash_evaluate
function ash_transform!(flm::AbstractArray{<:Complex},
                        f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    flm .= f
    spinsph_transform!(flm, s)
    change_signs!(flm, s)
    return flm
end
function ash_transform(f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    return ash_transform!(similar(f), f, s, lmax)
end

function ash_evaluate!(f::AbstractMatrix{<:Complex},
                       flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    f .= flm
    change_signs!(f, s)
    spinsph_evaluate!(f, s)
    return f
end
function ash_evaluate(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_evaluate!(similar(flm), flm, s, lmax)
end

export ash_eth!, ash_eth, ash_ethbar!, ash_ethbar
function ash_eth!(ðflm::AbstractArray{<:Complex}, flm::AbstractArray{<:Complex},
                  s::Integer, lmax::Integer)
    flm′ = ðflm
    flm′ .= flm
    ðflm .= spinsph_eth(flm′, s)
    return ðflm
end
function ash_eth(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_eth!(similar(flm), flm, s, lmax)
end

function ash_ethbar!(ð̄flm::AbstractArray{<:Complex},
                     flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    flm′ = ð̄flm
    flm′ .= flm
    ð̄flm .= spinsph_ethbar(flm′, s)
    return ð̄flm
end
function ash_ethbar(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_ethbar!(similar(flm), flm, s, lmax)
end
