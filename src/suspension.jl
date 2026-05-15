
using Printf
using Plots
using ProjectedSphericalHarmonics

include("bim/boundary.jl")
include("bim/operators.jl")
include("bim/solvers.jl")
include("simulation/plotting.jl")
include("simulation/simulate.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Nₛ    = 128
M     = 16
α     = π/2
a, b  = 1.0, 1.0
g(θ)  = cos(2θ)

Ω   = ellipse(a, b, M)
γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = arc(γ, γs, γss, Nₛ; dom=[α, 2π-α])
Γ_N = arc(γ, γs, γss, Nₛ; dom=[-α, α])

# ── Run ───────────────────────────────────────────────────────────────────────

Nparticles   = 2
centers      = [0.0 - 2im, 0.0 + 2im]
orientations = [0.0, 0.0]

println("Running simulation...")
traj_pos, traj_ang, traj_t = simulate(centers, orientations, Ω, Γ_D, Γ_N, g; t_end=1.0e8, plot_mode=:distance, plot_every=0.01, output_file="distances.dat")

readline()
