# FastSphericalHarmonics.jl

The
[`FastSphericalHarmonics.jl`](https://github.com/eschnett/FastSphericalHarmonics.jl)
package wraps the
[`FastTransforms.jl`](https://github.com/JuliaApproximation/FastTransforms.jl)
Julia package to calculate Spherical Harmonics.

`FastSphericalHarmonics.jl` is a powerful, efficient, and well thought
out package. Unfortunately, its user interface is difficult to use for
beginners, and its documentation is very technical.
`FastSphericalHarmonics.jl` provides functions and documentation that
are easier to use. It would be worthwhile to fold
`FastSphericalHarmonics.jl` into `FastTransforms.jl` at some point.

Most functions come in two versions, one that mutates its arguments
and one that allocates its output.

## Helper functions and types

```@docs
SpHTypes
sph_points
sph_lmax
```

## Scalar Spherical Harmonics

```@docs
sph_mode
sph_transform!
sph_transform
sph_evaluate!
sph_evaluate
```
