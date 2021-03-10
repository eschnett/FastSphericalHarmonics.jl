# FastSphericalHarmonics.jl

East-to-use Spherical Harmonics, based on
[`FastTransforms.jl`](https://github.com/JuliaApproximation/FastTransforms.jl)

[![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/FastSphericalHarmonics.jl/dev)
* [![GitHub
  CI](https://github.com/eschnett/FastSphericalHarmonics.jl/workflows/CI/badge.svg)](https://github.com/eschnett/FastSphericalHarmonics.jl/actions)
* [![Codecov](https://codecov.io/gh/eschnett/FastSphericalHarmonics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eschnett/FastSphericalHarmonics.jl)

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

The documentation lists the implemented functions as well as their
Julia signatures. Most functions come in two versions, one that
mutates its arguments and one that allocates its output.
