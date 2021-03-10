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

julia> lmax = 40;

julia> Θ, Φ = sph_points(lmax+1);
```

Load [Makie](http://makie.juliaplots.org/stable/) to create the images:
```Julia
julia> using CairoMakie

julia> ticks = (xticks=MultiplesTicks(4, π, "π"), yticks=MultiplesTicks(4, π, "π"));
```

Create images of some modes:
```Julia
julia> for l in 4:4, m in 0:l
           C = zeros(lmax+1, 2lmax+1)
           C[sph_mode(l,m)] = 1
           F = sph_evaluate(C)

           scene, layout = layoutscene(; resolution=(400, 200))
           axis = layout[1, 1] = Axis(scene; xticks=MultiplesTicks(4, π, "π"), yticks=MultiplesTicks(4, π, "π"))
           heatmap!(axis, Φ, Θ, F')
           scale!(scene, 1, 1)
           save("mode$l$m.png", scene)
       end
```

The `l=0`, `m≥0` modes of the Spherical Harmonics look like this:
- `l=0`, `m=0`: ![l=4, m=0 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode40.png)
- `l=0`, `m=1`: ![l=4, m=1 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode41.png)
- `l=0`, `m=2`: ![l=4, m=2 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode42.png)
- `l=0`, `m=3`: ![l=4, m=3 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode43.png)
- `l=0`, `m=4`: ![l=4, m=4 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode44.png)

### Example 3: Vector Spherical Harmonics

The "usual" (scalar) spherical harmonics can only represent scalar
fields on the sphere. Vector fields require a different basis, and can
be represented via [Vector Spherical
Harmonics](https://en.wikipedia.org/wiki/Vector_spherical_harmonics).
Vector fields on the sphere have two components (`θ` and `ϕ`), both of
which have coordinate singularities at the poles. There are two
different sets of vector spherical harmonics, the *gradient* and the
*curl* harmonics. Each vector field on the sphere can be decomposed
into two parts, a gradient part and a curl part, and as their name
indicates, the gradient and curl harmonics can represent the
respective parts. The gradient harmonics are defined via the gradient
of the scalar spherical harmonics, and the curl harmonics via their
(2d) curl. 

The `θ` and `ϕ` components of the gradient and curl harmonics are not
independent. For example, the following four components are all
proportional to `(cos θ) (cos ϕ)`:
- `θ` component of`l=1`, `m=1` gradient harmonic 
- `ϕ` component of`l=2`, `m=1` gradient harmonic 
- `ϕ` component of`l=1`, `m=1` curl harmonic 
- `θ` component of`l=2`, `m=1` curl harmonic 
(see e.g. [here](https://en.wikipedia.org/wiki/Vector_spherical_harmonics)).

`FastSphericalHarmonics.jl` follows `FastTransforms.jl`'s choice to
transform the `θ` and `ϕ` components of a vector field separately,
each yielding its own set of coefficients. The respective harmonic
coefficients are stored in the coefficient arrays in different ways.

Setup:
```Julia
julia> using FastSphericalHarmonics

julia> lmax = 4;
```

`θ` components, calculated by setting the component index `v=1` when
calling `sphv_mode`:
```Julia
julia> Cθ = Array{Any}(undef, lmax+1, 2lmax+1);

julia> for l in 1:lmax, m in -l:l
           Cθ[sphv_mode(l,m,1)] = (l,m)
       end

julia> Cθ
5×9 Array{Any,2}:
    (1, 0)  #undef      #undef     …  #undef     #undef      #undef
    (2, 0)     (1, -1)     (1, 1)        (3, 3)     (4, -4)     (4, 4)
    (3, 0)     (2, -1)     (2, 1)        (4, 3)  #undef      #undef
    (4, 0)     (3, -1)     (3, 1)     #undef     #undef      #undef
 #undef        (4, -1)     (4, 1)     #undef     #undef      #undef
```

`ϕ` components, calculated by setting the component index `v=2` when
calling `sphv_mode`:
```Julia
julia> Cϕ = Array{Any}(undef, lmax+1, 2lmax+1);

julia> for l in 1:lmax, m in -l:l
           Cϕ[sphv_mode(l,m,2)] = (l,m)
       end

julia> Cϕ
5×9 Array{Any,2}:
    (1, 0)     (1, -1)     (1, 1)  …     (3, 3)     (4, -4)     (4, 4)
    (2, 0)     (2, -1)     (2, 1)        (4, 3)  #undef      #undef
    (3, 0)     (3, -1)     (3, 1)     #undef     #undef      #undef
    (4, 0)     (4, -1)     (4, 1)     #undef     #undef      #undef
 #undef     #undef      #undef        #undef     #undef      #undef
```

To represent a generic vector field on the sphere, one needs to
calculate both its gradient and curl components. However, it is
unlikely that one encounters such a field superposition in physics
since it doesn't have a well-defined parity, and one usually knows
whether a vector field is of gradient or of curl type. Transforming
the two components of the vector field yields two sets of
coefficients.

The scalar field `cos θ` is the `l=1`, `m=0` spherical harmonic; it is
smooth at the poles of the sphere -- the limit value at e.g. `θ = 0`
is independent of `ϕ`. The field `sin θ` cannot be represented as a
spherical harmonic since it has a kink at `θ = 0`. It is, however, the
`θ` component of the vector field `[sin θ, 0]ᵀ`, which is a gradient
vector spherical harmonic. Let us test this:

Prepare the two components of the vector field `Fθ` and `Fϕ`:
```Julia
julia> using FastSphericalHarmonics

julia> lmax = 4;

julia> Θ, Φ = sph_points(lmax+1);

julia> Fθ = [sin(θ) for θ in Θ, ϕ in Φ]
5×9 Matrix{Float64}:
 0.309017  0.309017  0.309017  0.309017  …  0.309017  0.309017  0.309017
 0.809017  0.809017  0.809017  0.809017     0.809017  0.809017  0.809017
 1.0       1.0       1.0       1.0          1.0       1.0       1.0
 0.809017  0.809017  0.809017  0.809017     0.809017  0.809017  0.809017
 0.309017  0.309017  0.309017  0.309017     0.309017  0.309017  0.309017

julia> Fϕ = [0.0 for θ in Θ, ϕ in Φ]
5×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```

Transform them:
```Julia
julia> Cθ = sphv_transform(Fθ)
5×9 Matrix{Float64}:
  2.89441      0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 -1.1907e-16   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  1.73238e-16  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 -1.06997e-16  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0          0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> Cϕ = sphv_transform(Fϕ)
5×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
Apart from round-off error, these coefficients contain only one
non-zero component, namely the `l=1`, `m=0` gradient harmonic
component (see also the output of `sphv_mode` for `Cθ` above to see
which coefficients are stored where.)

Let us now look at the vector field `[0, sin θ]ᵀ`:
```Julia
julia> using FastSphericalHarmonics

julia> lmax = 4;

julia> Θ, Φ = sph_points(lmax+1);

julia> Fθ = [0.0 for θ in Θ, ϕ in Φ]
5×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> Fϕ = [sin(θ) for θ in Θ, ϕ in Φ]
5×9 Matrix{Float64}:
 0.309017  0.309017  0.309017  0.309017  …  0.309017  0.309017  0.309017
 0.809017  0.809017  0.809017  0.809017     0.809017  0.809017  0.809017
 1.0       1.0       1.0       1.0          1.0       1.0       1.0
 0.809017  0.809017  0.809017  0.809017     0.809017  0.809017  0.809017
 0.309017  0.309017  0.309017  0.309017     0.309017  0.309017  0.309017
```

Transform them:
```Julia
julia> Cθ = sphv_transform(Fθ)
5×9 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> Cϕ = sphv_transform(Fϕ)
5×9 Matrix{Float64}:
  2.89441      0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 -1.1907e-16   0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  1.73238e-16  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 -1.06997e-16  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
  0.0          0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
In these coeffients, only the `l=1`, `m=0` curl harmonic coefficient
is non-zero (apart from round-off errors).
