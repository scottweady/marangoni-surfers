
using Printf
using Plots
using ProjectedSphericalHarmonics

include("bim/boundary.jl")
include("bim/operators.jl")
include("bim/solvers.jl")
include("simulation/plotting.jl")
include("simulation/simulate.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

M_disk = 16
Nₛ    = 128
α     = π/2
g(θ)  = cos(2θ)
D     = disk(M_disk)

# ── Run ───────────────────────────────────────────────────────────────────────

Nparticles   = 2
centers      = [0.0 - 3.0im, 0.0 + 3.0im]
orientations = [0.0, 0.0]

println("Running simulation...")
traj_pos, traj_ang, traj_t = simulate(centers, orientations; t_end=1.0e8, tol=0.05, plot_mode=:concentration)

readline()
