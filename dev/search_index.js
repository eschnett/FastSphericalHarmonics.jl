var documenterSearchIndex = {"docs":
[{"location":"#FastSphericalHarmonics.jl","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"The FastSphericalHarmonics.jl package wraps the FastTransforms.jl Julia package to calculate Spherical Harmonics.","category":"page"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"FastSphericalHarmonics.jl is a powerful, efficient, and well thought out package. Unfortunately, its user interface is difficult to use for beginners, and its documentation is very technical. FastSphericalHarmonics.jl provides functions and documentation that are easier to use. It would be worthwhile to fold FastSphericalHarmonics.jl into FastTransforms.jl at some point.","category":"page"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"Most functions come in two versions, one that mutates its arguments and one that allocates its output.","category":"page"},{"location":"#Helper-functions-and-types","page":"FastSphericalHarmonics.jl","title":"Helper functions and types","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"SpHTypes\nsph_points\nsph_lmax","category":"page"},{"location":"#FastSphericalHarmonics.SpHTypes","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.SpHTypes","text":"const SpHTypes = Union{Float64,Complex{Float64}}\n\nThe types supported by FastSphericalHarmonics (the same types as for FastTransforms).\n\n\n\n\n\n","category":"type"},{"location":"#FastSphericalHarmonics.sph_points","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_points","text":"Θ, Φ = sph_points(N::Integer)\nΘ::Vector{Float64}\nΦ::Vector{Float64}\n\nCalculate the locations of points on the sphere when using N points in the θ (latitudinal) direction.\n\nIt is length(Θ) = N and length(Φ) = M where M = 2N-1.\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_lmax","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_lmax","text":"lmax = sph_lmax(N::Integer)\nlmax::Int\n\nCalculate the maximum l mode that can be represented with N points. It is lmax = N - 1.\n\n\n\n\n\n","category":"function"},{"location":"#Scalar-Spherical-Harmonics","page":"FastSphericalHarmonics.jl","title":"Scalar Spherical Harmonics","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"sph_mode\nsph_transform!\nsph_transform\nsph_evaluate!\nsph_evaluate","category":"page"},{"location":"#FastSphericalHarmonics.sph_mode","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_mode","text":"idx = sph_mode(l::Integer, m::Integer)\nidx::CartesianIndex{2}\n\nCalculate the Cartesian index idx for the l,m mode. This index can be used to access the coefficients, i.e. the result of sph_transform or the input to sph_evaluate.\n\nCoefficients are stored in a two-dimensional array. Not all array elements are used.\n\nSee also: sph_transform!, sph_transform, sph_evaluate!, sph_evaluate\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_transform!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_transform!","text":"sph_transform!(F::Array{T,2}) where {T<:SpHTypes}\n\nTransform an array of points F into spherical harmonics. This is an in-place transform, i.e. the array F will be overwritten by the coefficients. Use sph_transform for a non-mutating function.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sph_mode to calculate the location in the output coefficient array for a particular l,m mode.\n\nSee also: sph_transform, sph_evaluate!, sph_points, sph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_transform","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_transform","text":"C = sph_transform(F::Array{T,2}) where {T<:SpHTypes}\nC::Array{T,2}\n\nTransform an array of points F into spherical harmonics. You can use sph_transform! for more efficient a mutating function that overwrites its argument F.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sph_mode to calculate the location in the output coefficient array for a particular l,m mode.\n\nSee also: sph_transform!, sph_evaluate, sph_points, sph_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_evaluate!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_evaluate!","text":"sph_evaluate!(C::Array{T,2}) where {T<:SpHTypes}\n\nEvaluate an array of coefficients C on the grid points on a sphere. This is an in-place transform, i.e. the array C will be overwritten by the point values. Use sph_evaluate for a non-mutating function.\n\nUse sph_mode to calculate the location in the input coefficient array for a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sph_evaluate, sph_transform!, sph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sph_evaluate","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sph_evaluate","text":"F = sph_evaluate(C::Array{T,2}) where {T<:SpHTypes}\nF::Array{T,2}\n\nEvaluate an array of coefficients C on the grid points on a sphere. You can use sph_evaluate! for more efficient a mutating function that overwrites its argument C.\n\nUse sph_mode to calculate the location in the input coefficient array for a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sph_evaluate!, sph_transform, sph_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#Vector-Spherical-Harmonics","page":"FastSphericalHarmonics.jl","title":"Vector Spherical Harmonics","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"sphv_mode\nsphv_transform!\nsphv_transform\nsphv_evaluate!\nsphv_evaluate","category":"page"},{"location":"#FastSphericalHarmonics.sphv_mode","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sphv_mode","text":"idx = sphv_mode(l::Integer, m::Integer, v::Integer)\nidx::CartesianIndex{2}\n\nCalculate the Cartesian index idx for the v components of l,m mode. v=1 is for the θ component, v=2 for the ϕ component. This index can be used to access the coefficients, i.e. the result of sphv_transform or the input to sphv_evaluate.\n\nCoefficients are stored in a two-dimensional array. Not all array elements are used.\n\nSee also: sphv_transform!, sphv_transform, sphv_evaluate!, sphv_evaluate\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sphv_transform!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sphv_transform!","text":"sphv_transform!(F::Array{T,2}) where {T<:SpHTypes}\n\nTransform an array of points F into vector spherical harmonics. The θ and ϕ components of a vector field are transformed independently, both by calling this function.\n\nThis is an in-place transform, i.e. the array F will be overwritten by the coefficients. Use sphv_transform for a non-mutating function.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sphv_mode to calculate the location in the output coefficient array for a particular component of a particular l,m mode.\n\nSee also: sphv_transform, sphv_evaluate!, sph_points, sphv_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sphv_transform","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sphv_transform","text":"C = sphv_transform(F::Array{T,2}) where {T<:SpHTypes}\nC::Array{T,2}\n\nTransform an array of points F into vector spherical harmonics. The θ and ϕ components of a vector field are transformed independently, both by calling this function.\n\nYou can use sphv_transform! for more efficient a mutating function that overwrites its argument F.\n\nUse sph_points to caluclate the location of the points on the sphere for the input array F.\n\nUse sphv_mode to calculate the location in the output coefficient array for a particular component of a particular l,m mode.\n\nSee also: sphv_transform!, sphv_evaluate, sph_points, sphv_mode\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sphv_evaluate!","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sphv_evaluate!","text":"sphv_evaluate!(C::Array{T,2}) where {T<:SpHTypes}\n\nEvaluate an array of vector spherical harmonic coefficients C on the grid points on a sphere. The θ and ϕ components of a vector field are evaluated independently, both by calling this function for the respective sets of coefficients.\n\nThis is an in-place transform, i.e. the array C will be overwritten by the point values. Use sphv_evaluate for a non-mutating function.\n\nUse sphv_mode to calculate the location in the input coefficient array for a particular component of a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sphv_evaluate, sphv_transform!, sphv_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.sphv_evaluate","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.sphv_evaluate","text":"F = sphv_evaluate(C::Array{T,2}) where {T<:SpHTypes}\nF::Array{T,2}\n\nEvaluate an array of vector spherical harmonic coefficients C on the grid points on a sphere. The θ and ϕ components of a vector field are evaluated independently, both by calling this function for the respective sets of coefficients.\n\nYou can use sphv_evaluate! for more efficient a mutating function that overwrites its argument C.\n\nUse sphv_mode to calculate the location in the input coefficient array for a particular component of a particular l,m mode.\n\nUse sph_points to caluclate the location of the points on the sphere in the output array F.\n\nSee also: sphv_evaluate!, sphv_transform, sphv_mode, sph_points\n\n\n\n\n\n","category":"function"},{"location":"#Spin-Spherical-Harmonics","page":"FastSphericalHarmonics.jl","title":"Spin Spherical Harmonics","text":"","category":"section"},{"location":"","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.jl","text":"spinsph_transform\nspinsph_evaluate","category":"page"},{"location":"#FastSphericalHarmonics.spinsph_transform","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_transform","text":"spinsph_transform(F::Array{Complex{Float64},2}, s::Int)\n\nCalculate the spin spherical harmonic transformation with spin weight s.\n\n\n\n\n\n","category":"function"},{"location":"#FastSphericalHarmonics.spinsph_evaluate","page":"FastSphericalHarmonics.jl","title":"FastSphericalHarmonics.spinsph_evaluate","text":"spinsph_evaluate(C::Array{Complex{Float64},2}, s::Int)\n\nEvaluate the spin spherical harmonic transformation with spin weight s on points on the sphere.\n\n\n\n\n\n","category":"function"}]
}
