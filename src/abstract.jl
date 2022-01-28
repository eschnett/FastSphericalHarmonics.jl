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

export ash_mode_index
function ash_mode_index(s::Integer, l::Integer, m::Integer, lmax::Integer)
    0 ≤ lmax || throw(DomainError())
    abs(s) ≤ l ≤ lmax || throw(DomainError())
    -l ≤ m ≤ l || throw(DomainError())
    return spinsph_mode(s, l, m)
end

function change_signs!(flm::AbstractArray{<:Complex}, s::Integer)
    for (row, col) in CartesianIndices(flm)
        m = col == 1 ? 0 : (col % 2 == 0 ? -1 : 1) * (col ÷ 2)
        l = row + max(abs(s), abs(m)) - 1
        sign = (isodd(s) && m < -s) || (m >= -s && isodd(m)) ? -1 : +1
        flm[row, col] = flipsign(flm[row, col], sign)
    end
    return flm
end

export ash_transform!, ash_transform, ash_evaluate!, ash_evaluate
function ash_transform!(flm::AbstractArray{<:Complex},
                        f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    flm .= f
    spinsph_transform!(flm, s)
    return change_signs!(flm, s)
end
function ash_transform(f::AbstractMatrix{<:Complex}, s::Integer, lmax::Integer)
    return ash_transform!(similar(f), f, s)
end

function ash_evaluate!(f::AbstractMatrix{<:Complex},
                       flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    f .= flm
    change_signs!(f, s)
    return spinsph_evaluate!(f, s)
end
function ash_evaluate(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_evaluate!(similar(flm), flm, s)
end

export ash_eth!, ash_eth, ash_ethbar!, ash_ethbar
function ash_eth!(ðflm::AbstractArray{<:Complex}, flm::AbstractArray{<:Complex},
                  s::Integer, lmax::Integer)
    flm′ = ðflm
    flm′ .= flm
    change_signs!(flm′, s)
    ðflm .= spinsph_eth(flm′, s)
    return change_signs!(ðflm, s + 1)
end
function ash_eth(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_eth(similar(flm), flm, s)
end

function ash_ethbar!(ð̄flm::AbstractArray{<:Complex},
                     flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    flm′ = ð̄flm
    flm′ .= flm
    change_signs!(flm′, s)
    ð̄flm .= spinsph_ethbar(flm′, s)
    return change_signs!(ð̄flm, s - 1)
end
function ash_ethbar(flm::AbstractArray{<:Complex}, s::Integer, lmax::Integer)
    return ash_ethbar(similar(flm), flm, s)
end
