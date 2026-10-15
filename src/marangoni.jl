
dir = @__DIR__

include(joinpath(dir, "script_params.jl"))
include(joinpath(dir, "bie", "bie.jl"))

using Printf
using ProjectedSphericalHarmonics

include(joinpath(dir, "solvers.jl"))
include(joinpath(dir, "simulation.jl"))

using Plots

include(joinpath(dir, "plotting.jl"))