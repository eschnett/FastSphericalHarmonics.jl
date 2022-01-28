export ash_grid_size, ash_nmodes
function ash_grid_size(lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    M = 2 * N - 1
    return N, M
end
function ash_nmodes(lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    M = 2 * N - 1
    return N, M
end

export ash_ntheta, ash_nphi, ash_thetas, ash_phis, ash_point_coord,
       ash_point_delta, ash_grid_as_phi_theta
ash_ntheta(lmax) = ash_grid_size(lmax)[1]
ash_nphi(lmax) = ash_grid_size(lmax)[2]
function ash_thetas(lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    return sph_points(N)[1]
end
function ash_phis(lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    return sph_points(N)[2]
end
function ash_point_coord(ij::CartesianIndex{2}, lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    M = 2 * N - 1
    t, p = Tuple(ij)
    theta = (π / N * (0.5:(N - 0.5)))[t]
    phi = (2π / M * (0:(M - 1)))[p]
    return theta, phi
end
function ash_point_delta(ij::CartesianIndex{2}, lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    N = Int(lmax) + 1
    M = 2 * N - 1
    dtheta = π / N
    dphi = 2π / M
    return dtheta, dphi
end
ash_grid_as_phi_theta(grid::AbstractMatrix) = transpose(grid)

function ash_mode_index(s::Integer, l::Integer, m::Integer, lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    abs(s) ≤ l ≤ lmax || throw(DomainError())
    -l ≤ m ≤ l || throw(DomainError())
    return spinsph_mode(s, l, m)
end

export ash_transform!, ash_transform, ash_evaluate!, ash_evaluate
function ash_transform!(flm::AbstractArray{<:Complex},
                        f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    flm .= f
    return spinsph_transform!(flm, s)
end
function ash_transform(f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    return spinsph_transform(f, s)
end

function ash_evaluate!(f::AbstractMatrix{<:Complex},
                       flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    f .= flm
    return spinsph_evaluate!(f, s)
end
function ash_evaluate(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return spinsph_evaluate(flm, s)
end

export ash_eth!, ash_eth, ash_ethbar!, ash_ethbar
function ash_eth!(ðflm::AbstractArray{<:Complex}, flm::AbstractArray{<:Complex},
                  s::Integer, lmax::Integer)
    ðflm .= spinsph_eth(flm, s)
    return ðflm
end
function ash_eth(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return spinsph_eth(flm, s)
end

function ash_ethbar!(ðflm::AbstractArray{<:Complex},
                     flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    ð̄flm .= spinsph_ethbar(flm, s)
    return ð̄flm
end
function ash_ethbar(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return spinsph_ethbar(flm, s)
end
