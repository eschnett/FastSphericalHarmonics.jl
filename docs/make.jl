# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using FastSphericalHarmonics

makedocs(; sitename="FastSphericalHarmonics", format=Documenter.HTML(),
         modules=[FastSphericalHarmonics])

deploydocs(; repo="github.com/eschnett/FastSphericalHarmonics.jl.git",
           devbranch="main", push_preview=true)
