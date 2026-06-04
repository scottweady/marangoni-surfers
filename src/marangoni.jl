
using FFTW
using FastGaussQuadrature
using FMM3D
using LinearAlgebra
using ProjectedSphericalHarmonics

using Printf
using Plots

dir = @__DIR__

include(joinpath(dir, "boundary.jl"))
include(joinpath(dir, "quadrature.jl"))
include(joinpath(dir, "matrices.jl"))
include(joinpath(dir, "operators.jl"))
include(joinpath(dir, "solvers.jl"))
include(joinpath(dir, "simulation.jl"))
include(joinpath(dir, "plotting.jl"))