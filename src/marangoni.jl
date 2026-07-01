
using Printf
using Plots

dir = @__DIR__

include(joinpath(dir, "bie", "bie.jl"))
include(joinpath(dir, "marangoni_solvers.jl"))
include(joinpath(dir, "simulation.jl"))
include(joinpath(dir, "plotting.jl"))