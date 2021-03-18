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
Julia package to calculate spherical harmonics.

`FastSphericalHarmonics.jl` is a powerful, efficient, and well thought
out package. Unfortunately, its user interface is difficult to use for
beginners, and its documentation is very technical.
`FastSphericalHarmonics.jl` provides functions and documentation that
are easier to use. It would be worthwhile to fold
`FastSphericalHarmonics.jl` into `FastTransforms.jl` at some point.

## Features

This package provides scalar and spin spherical harmonic
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
5×9 Matrix{Float64}:
 -0.556027   -1.1449    1.08238   -0.886205  -1.05099    0.168341   0.36901      0.681085  -0.145211
 -0.444383   -0.468606  0.187028   0.684565   0.502079   0.284259  -0.00761298  -1.33913    0.642896
  0.0271553   0.156143  0.518149  -1.59058   -0.216248   0.569829   0.562669    -0.238284   1.81935
 -0.299484   -2.64199   1.49138    0.410653  -0.706424  -1.42206    0.106869     1.01936   -0.36726
  1.77786     1.00331   0.367563  -0.85635   -3.86593   -0.372401   0.569458     0.701771   0.756569
```

Transform to spherical harmonics:
```
julia> C = sph_transform(F)
5×9 Matrix{Float64}:
 -0.0971079  -0.480385    0.250188   -0.237836  -0.344185  -1.15475    -0.390183  -0.828615     -0.0736505
  0.173019    0.399055   -0.601621    0.235944   1.15235   -0.0698683  -0.222073   0.0436882    -1.03784
 -0.35473     0.28504    -0.0294237  -0.492649  -0.931321  -0.478315    0.684783  -0.000600234   0.577805
 -0.309621   -0.0516482  -0.99214    -0.318465  -0.243649   0.0434071  -0.452602  -0.338297      0.332204
  0.296832   -1.16363     1.65583     0.606001  -0.281404  -0.555203   -0.424356   0.21506       0.123637
```

Note that the coefficient array `C` contains `(lmax+1) * (2lmax+1) =
45` coefficients, more than the `(lmax+1)^2 = 25` coefficients we
expected. There are some higher modes (with `l > lmax`) present as
well. This makes the spherical harmonic transform invertible:
```Julia
julia> F′ = sph_evaluate(C)
5×9 Matrix{Float64}:
 -0.556027   -1.1449    1.08238   -0.886205  -1.05099    0.168341   0.36901      0.681085  -0.145211
 -0.444383   -0.468606  0.187028   0.684565   0.502079   0.284259  -0.00761298  -1.33913    0.642896
  0.0271553   0.156143  0.518149  -1.59058   -0.216248   0.569829   0.562669    -0.238284   1.81935
 -0.299484   -2.64199   1.49138    0.410653  -0.706424  -1.42206    0.106869     1.01936   -0.36726
  1.77786     1.00331   0.367563  -0.85635   -3.86593   -0.372401   0.569458     0.701771   0.756569
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
5×9 Matrix{Float64}:
 -0.0971079  -0.480385    0.250188   -0.237836  -0.344185  -1.15475    -0.390183  -0.828615  -0.0736505
  0.173019    0.399055   -0.601621    0.235944   1.15235   -0.0698683  -0.222073   0.0        0.0
 -0.35473     0.28504    -0.0294237  -0.492649  -0.931321   0.0         0.0        0.0        0.0
 -0.309621   -0.0516482  -0.99214     0.0        0.0        0.0         0.0        0.0        0.0
  0.296832    0.0         0.0         0.0        0.0        0.0         0.0        0.0        0.0
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
5×9 Matrix{Any}:
 (0, 0)     (1, -1)     (1, 1)     (2, -2)     (2, 2)     (3, -3)     (3, 3)     (4, -4)     (4, 4)
 (1, 0)     (2, -1)     (2, 1)     (3, -2)     (3, 2)     (4, -3)     (4, 3)  #undef      #undef
 (2, 0)     (3, -1)     (3, 1)     (4, -2)     (4, 2)  #undef      #undef     #undef      #undef
 (3, 0)     (4, -1)     (4, 1)  #undef      #undef     #undef      #undef     #undef      #undef
 (4, 0)  #undef      #undef     #undef      #undef     #undef      #undef     #undef      #undef
```

Evaluating these modes gives us scalar field values `F″` that contain
only these modes, which are now different from the original values in
`F`:
```Julia
julia> F″ = sph_evaluate(C′)
5×9 Matrix{Float64}:
 -1.09395   -0.814816  -0.0626794   0.497403   0.649439    0.440388   0.0488995  -0.36331   -0.783888
 -0.188812  -0.38777    0.535328   -0.220322   0.260065    0.224342  -0.309924   -0.529559   0.657753
  0.290636  -0.415552   0.643938   -1.07714    0.0812705   0.622275   0.631893   -0.554523   1.38539
 -1.23983   -1.27778    0.806158    0.115463  -1.27163    -1.4813     0.322699    0.908434   0.708826
  0.485627   0.410934   0.464254   -0.108253  -1.02271    -1.1944    -0.337374    0.589485   0.794295
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

The `l=0`, `m≥0` modes of the spherical harmonics look like this:
- `l=0`, `m=0`: ![l=4, m=0 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode40.png)
- `l=0`, `m=1`: ![l=4, m=1 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode41.png)
- `l=0`, `m=2`: ![l=4, m=2 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode42.png)
- `l=0`, `m=3`: ![l=4, m=3 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode43.png)
- `l=0`, `m=4`: ![l=4, m=4 mode](https://github.com/eschnett/FastSphericalHarmonics.jl/blob/main/figures/mode44.png)

### Example 3: Laplacian of Spherical Harmonics

We use both scalar and spin-weighted spherical harmonics to calculate
derivatives. Spin-weighted spherical harmonics with spin-weight 0 are
the same as scalar spherical harmonics, except with a different
normalization.

Some preliminaries:
```Julia
julia> using FastSphericalHarmonics

julia> chop(x) = abs2(x) < 100eps(x) ? zero(x) : x;

julia> chop(x::Complex) = Complex(chop(real(x)), chop(imag(x)));
```

Choose a function (here `z + 2x` with `z = cos θ` and `x = sin θ cos
ϕ`:
```Julia
julia> lmax = 4;

julia> Θ, Φ = sph_points(lmax+1);

julia> F = [cos(θ) + sin(θ)*cos(ϕ) for θ in Θ, ϕ in Φ]
5×9 Matrix{Float64}:
  1.26007    1.18778     1.00472   …   0.796548   1.00472    1.18778
  1.3968     1.20753     0.72827       0.183277   0.72827    1.20753
  1.0        0.766044    0.173648     -0.5        0.173648   0.766044
  0.221232   0.0319577  -0.447301     -0.992294  -0.447301   0.0319577
 -0.64204   -0.714336   -0.897396     -1.10557   -0.897396  -0.714336
```

We transform to scalar spherical harmonics, calculate the Laplacian
(which is efficient in spectral space), and convert back to point
values:
```Julia
julia> C = sph_transform(F); chop.(C)
5×9 Matrix{Float64}:
 0.0      0.0  2.04665  0.0  0.0  0.0  0.0  0.0  0.0
 2.04665  0.0  0.0      0.0  0.0  0.0  0.0  0.0  0.0
 0.0      0.0  0.0      0.0  0.0  0.0  0.0  0.0  0.0
 0.0      0.0  0.0      0.0  0.0  0.0  0.0  0.0  0.0
 0.0      0.0  0.0      0.0  0.0  0.0  0.0  0.0  0.0

julia> ΔC = sph_laplace(C); chop.(ΔC)
5×9 Matrix{Float64}:
  0.0      0.0  -4.09331  0.0  0.0  0.0  0.0  0.0  0.0
 -4.09331  0.0   0.0      0.0  0.0  0.0  0.0  0.0  0.0
  0.0      0.0   0.0      0.0  0.0  0.0  0.0  0.0  0.0
  0.0      0.0   0.0      0.0  0.0  0.0  0.0  0.0  0.0
  0.0      0.0   0.0      0.0  0.0  0.0  0.0  0.0  0.0

julia> ΔF = sph_evaluate(ΔC)
5×9 Matrix{Float64}:
 -2.52015   -2.37555    -2.00943   …  -1.5931    -2.00943   -2.37555
 -2.7936    -2.41506    -1.45654      -0.366554  -1.45654   -2.41506
 -2.0       -1.53209    -0.347296      1.0       -0.347296  -1.53209
 -0.442463  -0.0639154   0.894602      1.98459    0.894602  -0.0639154
  1.28408    1.42867     1.79479       2.21113    1.79479    1.42867
```

Since we chose `F` to consist of two `l=1` modes, we know its
Laplacian: `ΔF = -l(l+1) F`:
```Julia
julia> ΔF ≈ -2F
true
```

Now let's do the same calculation with spin-weighted spherical
harmonics. These are defined for complex functions, so we first
convert the real array to a complex array, and then transform to
spin-weighted spherical harmonics (of spin-weight `0`).
```Julia
julia> F⁰ = Complex.(F);

julia> C⁰ = spinsph_transform(F⁰, 0); chop.(C⁰)
5×9 Matrix{ComplexF64}:
     0.0+0.0im  1.4472+0.0im  …  0.0+0.0im  0.0+0.0im  0.0+0.0im
 2.04665+0.0im     0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im     0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im     0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im     0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
```
Due to the way in which `FastTransforms.jl` defines spherical
harmonics, and because we started with a real function `F`, the
imaginary part of `C⁰` is actually zero.

We test that evaluating the spin-weighted spherical harmonics gives us
back the original function `F`:
```Julia
julia> spinsph_evaluate(C⁰, 0) ≈ F
true
```

We then apply the ð (["eth"](https://en.wikipedia.org/wiki/Eth))
operator, which is the gradient when applied to a spin-0 function,
yielding a spin-1 function.
```Julia
julia> ðC¹ = spinsph_eth(C⁰, 0); chop.(ðC¹)
5×9 Matrix{ComplexF64}:
     0.0+0.0im  -2.04665+0.0im  …  0.0+0.0im  0.0+0.0im  0.0+0.0im
 2.89441+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
     0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
```

We can evaluate this spin-1 function on the points on the sphere and
thus read off the gradient of `F`:
```Julia
julia> ðF¹ = spinsph_evaluate(ðC¹, 1);

julia> ∂θF = real(ðF¹)
5×9 Matrix{Float64}:
  0.657164      0.28605    0.0885849  …   1.15716    1.22574     1.02828
  1.06331       0.6922     0.494734       1.56331    1.63189     1.43443
  1.37051e-16  -0.371114  -0.568579       0.5        0.568579    0.371114
 -1.06331      -1.43443   -1.63189       -0.563314  -0.494734   -0.6922
 -0.657164     -1.02828   -1.22574       -0.157164  -0.0885849  -0.28605

julia> cscθ∂ϕF = imag(ðF¹)
5×9 Matrix{Float64}:
  8.66754e-16  0.642788  0.984808  …  -0.866025  -0.984808  -0.642788
 -5.4187e-17   0.642788  0.984808     -0.866025  -0.984808  -0.642788
 -2.88227e-16  0.642788  0.984808     -0.866025  -0.984808  -0.642788
  1.87818e-16  0.642788  0.984808     -0.866025  -0.984808  -0.642788
 -4.87207e-16  0.642788  0.984808     -0.866025  -0.984808  -0.642788
```
We thus have `julia> ∇F = (∂θF, sin(θ) * cscθ∂ϕF)`. However, since
`sin(θ)` has a coordinate singularity at the pole, it's best not to
actually perform this multiplication.

Next we apply the ð̄ ("eth-bar") operator, which is the divergence when
applied to a spin-1 function, yielding a spin-0 function again, which
we evaluate on the grid points.
```Julia
julia> ð̄ðC⁰ = spinsph_ethbar(ðC¹, 1); chop.(ð̄ðC⁰)
5×9 Matrix{ComplexF64}:
      0.0+0.0im  -2.89441+0.0im  …  0.0+0.0im  0.0+0.0im  0.0+0.0im
 -4.09331+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
      0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
      0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im
      0.0+0.0im       0.0+0.0im     0.0+0.0im  0.0+0.0im  0.0+0.0im

julia> ð̄ðF⁰ = spinsph_evaluate(ð̄ðC⁰, 0); chop.(ð̄ðF⁰)
5×9 Matrix{ComplexF64}:
  -2.52015+0.0im    -2.37555+0.0im  …   -2.00943+0.0im    -2.37555+0.0im
   -2.7936+0.0im    -2.41506+0.0im      -1.45654+0.0im    -2.41506+0.0im
      -2.0+0.0im    -1.53209+0.0im     -0.347296+0.0im    -1.53209+0.0im
 -0.442463+0.0im  -0.0639154+0.0im      0.894602+0.0im  -0.0639154+0.0im
   1.28408+0.0im     1.42867+0.0im       1.79479+0.0im     1.42867+0.0im
```

This function `ð̄ðF⁰` is the Laplacian of our original function `F`
above. It is complex, but since we started with a real function `F`,
`ð̄ðF⁰` has a zero imaginay part (up to round-off):
```Julia
julia> maximum(abs.(imag.(ð̄ðF⁰)))
1.1102230246251565e-16
```

Of course, both ways of evaluating the Laplacian give the same result:
```Julia
julia> real.(ð̄ðF⁰) ≈ ΔF
true
```

We can also apply the ð̄ operator first, and then the ð operator:
```Julia
julia> ð̄C⁻¹ = spinsph_ethbar(C⁰, 0);

julia> ðð̄C⁰ = spinsph_eth(ð̄C⁻¹, -1);

julia> ðð̄F⁰ = spinsph_evaluate(ðð̄C⁰, 0);

julia> ðð̄F⁰ ≈ ΔF
true
```
