# FastSphericalHarmonics.jl

East-to-use Spherical Harmonics, based on
[`FastTransforms.jl`](https://github.com/JuliaApproximation/FastTransforms.jl)

* [![Documenter](https://img.shields.io/badge/docs-dev-blue.svg)](https://eschnett.github.io/FastSphericalHarmonics.jl/dev)
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

## Features

This package provides scalar, vector, and spin Spherical Harmonic
transforms for real and complex fields. The normalizations and
conventions are chosen to be convenient for real fields, and are
different from those usually used for complex spherical harmonics.
This is documented by
[FastTransforms.jl](https://juliaapproximation.github.io/FastTransforms.jl/stable/),
its underlying C implementation
[FastTransforms](https://mikaelslevinsky.github.io/FastTransforms/),
and in the references listed there (see bottom of the page there).

The
[documentation](https://eschnett.github.io/FastSphericalHarmonics.jl/dev)
lists the implemented functions as well as their Julia signatures.
Most functions come in two versions, one that mutates its arguments
and one that allocates its output. See also the test cases for more
examples.

## Examples

### Example 1: Transform, modify coefficients, and transform back

Make the example reproducible by choosing a particular random number
seed:
```Julia
julia> using Random

julia> Random.seed!(42);
```

Load the package and create a random scalar field on the sphere:
```Julia
julia> using FastSphericalHarmonics

julia> lmax = 4;

julia> F = randn(lmax+1, 2lmax+1)
5×9 Array{Float64,2}:
 -0.556027   -1.1449    1.08238   …   0.562669  -0.238284   1.81935
 -0.444383   -0.468606  0.187028      0.106869   1.01936   -0.36726
  0.0271553   0.156143  0.518149      0.569458   0.701771   0.756569
 -0.299484   -2.64199   1.49138       0.681085  -0.145211   0.0870168
  1.77786     1.00331   0.367563     -1.33913    0.642896  -0.851194
```

Transform to Spherical Harmonics:
```
julia> C = sph_transform(F)
5×9 Array{Float64,2}:
 -0.146293   -0.629688    0.223526     …  -0.216793  -0.777218    0.455651
  0.262666    0.0685201   0.000199464     -0.172973  -0.0955985  -0.421988
 -0.299492    0.258348   -0.090961         0.530321   0.156805    0.534255
  0.160477   -0.418951   -1.1076          -0.380814  -0.76739    -0.302279
  0.0249693  -0.668195    1.88872         -0.624655   0.382056    0.112519
```

Note that the coefficient array `C` contains `(lmax+1) * (2lmax+1) =
45` coefficients, more than the `(lmax+1)^2 = 25` coefficients we
expected. There are some higher modes (with `l > lmax`) present as
well. This makes the Spherical Harmonic transform invertible:
```Julia
julia> F′ = sph_evaluate(C)
5×9 Array{Float64,2}:
 -0.556027   -1.1449    1.08238   …   0.562669  -0.238284   1.81935
 -0.444383   -0.468606  0.187028      0.106869   1.01936   -0.36726
  0.0271553   0.156143  0.518149      0.569458   0.701771   0.756569
 -0.299484   -2.64199   1.49138       0.681085  -0.145211   0.0870168
  1.77786     1.00331   0.367563     -1.33913    0.642896  -0.851194
```
Evaluating the modes at the grid points gives the same values back.

Let's create a new coefficient array `C′` that contains only some of
the modes, and leaves all other modes set to zero:
```Julia
julia> C′ = zeros(size(C));

julia> for l in 0:lmax, m in -l:l
           C′[sph_mode(l,m)] = C[sph_mode(l,m)]
       end

julia> C′
5×9 Array{Float64,2}:
 -0.146293   -0.629688    0.223526     …  -0.216793  -0.777218  0.455651
  0.262666    0.0685201   0.000199464     -0.172973   0.0       0.0
 -0.299492    0.258348   -0.090961         0.0        0.0       0.0
  0.160477   -0.418951   -1.1076           0.0        0.0       0.0
  0.0249693   0.0         0.0              0.0        0.0       0.0
```

We can examine which of the coefficient array elements contains what
mode, and which are unused (or rather, used by the supernumerary
higher modes):
```Julia
julia> lm = similar(C, Any);

julia> for l in 0:lmax, m in -l:l
           lm[sph_mode(l,m)] = (l,m)
       end

julia> lm
5×9 Array{Any,2}:
 (0, 0)     (1, -1)     (1, 1)  …     (3, 3)     (4, -4)     (4, 4)
 (1, 0)     (2, -1)     (2, 1)        (4, 3)  #undef      #undef
 (2, 0)     (3, -1)     (3, 1)     #undef     #undef      #undef
 (3, 0)     (4, -1)     (4, 1)     #undef     #undef      #undef
 (4, 0)  #undef      #undef        #undef     #undef      #undef
```

Evaluating these modes gives us scalar field values `F″` that contain
only these modes, which are now different from the original values in
`F`:
```Julia
julia> F″ = sph_evaluate(C′)
5×9 Array{Float64,2}:
 -0.949163  -0.813295  -0.216166  …   0.547065   0.424954   -0.349691
 -0.445106   0.13906    0.151119     -0.152408   1.47035    -0.0645492
  0.185819  -0.525536   0.658583      0.725372   0.425367    0.407101
 -0.721984  -1.55685    1.17405       0.634849  -0.0154634   0.691279
  0.102395   0.142327   0.254127     -0.955378  -0.395419    0.0657883
```

### Example 2: Visualize some modes

Load the library and define the resolution of the images:
```Julia
julia> using FastSphericalHarmonics

julia> lmax = 100;

julia> Θ, Φ = sph_points(lmax+1);
```

Load [Makie](http://makie.juliaplots.org/stable/) to create the images:
```Julia
julia> using CairoMakie

julia> axis = (xticks=MultiplesTicks(4, π, "π"), yticks=MultiplesTicks(4, π, "π"));
```

Create images of some modes
```Julia
julia> for l in 4:4, m in 0:l
           C = zeros(lmax+1, 2lmax+1)
           C[sph_mode(4,2)] = 1
           F = sph_evaluate(C)
           save("mode$l$m.pdf", heatmap(Φ, Θ, F'; axis=axis))
       end
```

![](https://github.com/eschnett/FastSphericalHarmonics.jl/figures/mode40.pdf | width=100)
![](https://github.com/eschnett/FastSphericalHarmonics.jl/figures/mode41.pdf | width=100)
![](https://github.com/eschnett/FastSphericalHarmonics.jl/figures/mode42.pdf | width=100)
![](https://github.com/eschnett/FastSphericalHarmonics.jl/figures/mode43.pdf | width=100)
![](https://github.com/eschnett/FastSphericalHarmonics.jl/figures/mode44.pdf | width=100)
