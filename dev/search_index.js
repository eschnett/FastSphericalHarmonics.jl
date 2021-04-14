var documenterSearchIndex = {"docs":
[{"location":"#FastSphericalHarmonics.jl","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"The FastSphericalHarmonics.jl package wraps the FastTransforms.jl Julia package to calculate Spherical Harmonics.","category":"page"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"FastSphericalHarmonics.jl is a powerful, efficient, and well thought out package. Unfortunately, its user interface is difficult to use for beginners, and its documentation is very technical. FastSphericalHarmonics.jl provides functions and documentation that are easier to use. It would be worthwhile to fold FastSphericalHarmonics.jl into FastTransforms.jl at some point.","category":"page"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"Most functions come in two versions, one that mutates its arguments and one that allocates its output.","category":"page"},{"location":"#Helper-functions-and-types","page":"FastSphericalHarmonics.jl","title":"Helper functions and types","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"SpHTypes\nsph_points\nsph_lmax","category":"page"},{"location":"#FastSphericalHarmonics.SpHTypes","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.SpHTypes","text":"const SpHTypes = Union{Float64,Complex{Float64}}\n\nThe types supported by FastSphericalHarmonics (the same types as for FastTransforms).\n\n\n\n\n\n","category":"type"},{"location":"#FastSphericalHarmonics.sph_points","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_points","text":"Θ, Φ = sph_points(N::Integer)\nΘ::Vector{Float64}\nΦ::Vector{Float64}\n\nCalculate the locations of points on the sphere when using N points in the θ (latitudinal) direction.\n\nIt is length(Θ) = N and length(Φ) = M where M = 2N-1.\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_lmax","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_lmax","text":"lmax = sph_lmax(N::Integer)\nlmax::Int\n\nCalculate the maximum l mode that can be represented with N points. It is lmax = N - 1.\n\n\n\n\n\n","category":"function"},{"location":"#Scalar-Spherical-Harmonics","page":"FastSphericalHarmonics.jl","title":"Scalar Spherical Harmonics","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"sph_mode\nsph_transform!\nsph_transform\nsph_evaluate!\nsph_evaluate\nsph_laplace!\nsph_laplace","category":"page"},{"location":"#FastSphericalHarmonics.sph_mode","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_mode","text":"idx = sph_mode(l::Integer, m::Integer)\nidx::CartesianIndex{2}\n\nCalculate the Cartesian index idx for the l,m mode. This index can be used to access the coefficients, i.e. the result of sph_transform or the input to sph_evaluate.\n\nCoefficients are stored in a two-dimensional array. Not all array elements are used. See this page, section \"sph2fourier\", for details.\n\nSee also: sph_transform!, sph_transform, sph_evaluate!, sph_evaluate\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_transform!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_transform!","text":"sph_transform!(F::Array{T,2}) where {T<:SpHTypes}\n\nTransform an array of points F into spherical harmonics. This is an in-place transform, i.e. the array F will be overwritten by the coefficients. Use sph_transform for a non-mutating function.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sph_mode to calculate the location in the output coefficient array for a particular l,m mode.\n\nSee also: sph_transform, sph_evaluate!, sph_points, sph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_transform","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_transform","text":"C = sph_transform(F::AbstractArray{T,2}) where {T<:SpHTypes}\nC::Array{T,2}\n\nTransform an array of points F into spherical harmonics. You can use sph_transform! for more efficient a mutating function that overwrites its argument F.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sph_mode to calculate the location in the output coefficient array for a particular l,m mode.\n\nSee also: sph_transform!, sph_evaluate, sph_points, sph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_evaluate!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_evaluate!","text":"sph_evaluate!(C::Array{T,2}) where {T<:SpHTypes}\n\nEvaluate an array of coefficients C on the grid points on a sphere. This is an in-place transform, i.e. the array C will be overwritten by the point values. Use sph_evaluate for a non-mutating function.\n\nUse sph_mode to calculate the location in the input coefficient array for a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sph_evaluate, sph_transform!, sph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_evaluate","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_evaluate","text":"F = sph_evaluate(C::AbstractArray{T,2}) where {T<:SpHTypes}\nF::Array{T,2}\n\nEvaluate an array of coefficients C on the grid points on a sphere. You can use sph_evaluate! for more efficient a mutating function that overwrites its argument C.\n\nUse sph_mode to calculate the location in the input coefficient array for a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sph_evaluate!, sph_transform, sph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_laplace!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_laplace!","text":"sph_laplace!(C::Array{T,2}) where {T<:SpHTypes}\n\nCalculate the Laplacian of a set of coefficients C. This is an in-place transform, i.e. the array C will be overwritten by the result. Use sph_laplace for a non-mutating function.\n\nSee also: sph_transform!, sph_evaluate!, sph_laplace\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_laplace","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_laplace","text":"ΔC = sph_laplace(C::AbstractArray{T,2}) where {T<:SpHTypes}\nΔC::Array{T,2}\n\nCalculate the Laplacian of a set of coefficients C. You can use sph_laplace! for more efficient a mutating function that overwrites its argument C.\n\nSee also: sph_transform, sph_evaluate, sph_laplace!\n\n\n\n\n\n","category":"function"},{"location":"#Spin-Spherical-Harmonics","page":"FastSphericalHarmonics.jl","title":"Spin Spherical Harmonics","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"spinsph_mode\nspinsph_transform!\nspinsph_transform\nspinsph_evaluate!\nspinsph_evaluate\nspinsph_eth\nspinsph_ethbar","category":"page"},{"location":"#FastSphericalHarmonics.spinsph_mode","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_mode","text":"idx = spinsph_mode(l::Integer, m::Integer, s::Integer)\nidx::CartesianIndex{2}\n\nCalculate the Cartesian index idx for the l,m mode for spin weight s. This index can be used to access the coefficients, i.e. the result of spinsph_transform or the input to spinsph_evaluate.\n\nCoefficients are stored in a two-dimensional array. Not all array elements are used. See this page, section \"spinsph2fourier\", for details.\n\nSee also: spinsph_transform!, spinsph_transform, spinsph_evaluate!, spinsph_evaluate\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_transform!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_transform!","text":"spinsph_transform!(F::Array{Complex{Float64},2}, s::Int)\n\nCalculate the spin spherical harmonic transformation with spin weight s. This is an in-place transform, i.e. the array F will be overwritten by the coefficients. Use spinsph_transform for a non-mutating function.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse spinsph_mode to calculate the location in the output coefficient array for a particular l,m mode with spin weight s\u0001.\n\nSee also: spinsph_transform, spinsph_evaluate!, sph_points, spinsph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_transform","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_transform","text":"C = spinsph_transform(F::AbstractArray{Complex{Float64},2}, s::Int)\nC::Array{Complex{Float64},2}\n\nCalculate the spin spherical harmonic transformation with spin weight s. You can use spinsph_transform! for more efficient a mutating function that overwrites its argument F.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse spinsph_mode to calculate the location in the output coefficient array for a particular l,m mode with spin weight s.\n\nSee also: spinsph_transform!, spinsph_evaluate, sph_points, spinsph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_evaluate!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_evaluate!","text":"spinsph_evaluate!(C::Array{Complex{Float64},2}, s::Int)\n\nEvaluate the spin spherical harmonic transformation with spin weight s on points on the sphere. This is an in-place transform, i.e. the array C will be overwritten by the point values. Use spinsph_evaluate for a non-mutating function.\n\nUse spinsph_mode to calculate the location in the input coefficient array for a particular l,m mode with spin weight s.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: spinsph_evaluate, spinsph_transform!, spinsph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_evaluate","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_evaluate","text":"F = spinsph_evaluate(C::AbstractArray{Complex{Float64},2}, s::Int)\nF::Array{Complex{Float64},2}\n\nEvaluate the spin spherical harmonic transformation with spin weight s on points on the sphere. You can use spinsph_evaluate! for more efficient a mutating function that overwrites its argument C.\n\nUse spinsph_mode to calculate the location in the input coefficient array for a particular l,m mode with spin weight s.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: spinsph_evaluate!, spinsph_transform, spinsph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_eth","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_eth","text":"ðC = spinsph_eth(C::AbstractArray{Complex{Float64},2}, s::Int)\nðC::Array{Complex{Float64},2}\n\nApply the differential operator ð (\"eth\") to the coefficients C. This raises the spin weight s by 1. For a real function of spin weight 0, this is equivalent to calculating the gradient.\n\nThis function is the converse of spinsph_ethbar, which is a derivative operator lowering the spin weight.\n\nUse spinsph_mode to calculate the location in the output coefficient array for a particular l,m mode with spin weight s.\n\nSee also: spinsph_ethbar, spinsph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_ethbar","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_ethbar","text":"ðC = spinsph_ethbar(C::AbstractArray{Complex{Float64},2}, s::Int)\nðC::Array{Complex{Float64},2}\n\nApply the differential operator ð̄ (\"eth-bar\") to the coefficients C. This lowers the spin weight s by 1. For a function of spin weight 1, this is equivalent to calculating the divergence.\n\nThis function is the converse of spinsph_eth, which is a derivative operator raising the spin weight.\n\nUse spinsph_mode to calculate the location in the output coefficient array for a particular l,m mode with spin weight s.\n\nSee also: spinsph_eth, spinsph_mode\n\n\n\n\n\n","category":"function"}]
}
