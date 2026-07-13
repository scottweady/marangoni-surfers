
# 2D boundary integral equation tools: quadrature, boundary representation,
# kernel matrix assembly, layer potentials, and close evaluation.

using FastGaussQuadrature
using FFTW
using FMM3D
using LinearAlgebra
using SpecialFunctions

let dir = @__DIR__
    include(joinpath(dir, "boundary.jl"))
    include(joinpath(dir, "quadrature.jl"))
    include(joinpath(dir, "matrices.jl"))
    include(joinpath(dir, "operators.jl"))
    include(joinpath(dir, "solvers.jl"))
    include(joinpath(dir, "closeeval.jl"))
end
