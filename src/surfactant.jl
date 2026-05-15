
using Plots

include("bim/boundary.jl")
include("bim/operators.jl")
include("bim/solvers.jl")
include("simulation/plotting.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Nₛ   = 256
α    = π/2
a, b = 2.0, 1.0
g(θ) = cos(2θ)

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = arc(γ, γs, γss, Nₛ; dom=[α, 2π-α])
Γ_N = arc(γ, γs, γss, Nₛ; dom=[-α, α])

# ── Solve ─────────────────────────────────────────────────────────────────────
N_body   = 13
r0, dr   = 8.0, 4.0          # inner radius, radial growth per turn
θ_spiral = range(0, 4π, length=N_body+1)[1:end-1]
r        = r0 .+ dr .* θ_spiral ./ (2π)
centers  = r .* exp.(im .* θ_spiral)
# tangent angle of the spiral at each point
drdθ         = dr / (2π)
orientations = atan.(drdθ .* sin.(θ_spiral) .+ r .* cos.(θ_spiral),
                     drdθ .* cos.(θ_spiral) .- r .* sin.(θ_spiral))

σs, boundaries, ϕinf = surfactant_solve(centers, orientations, Γ_D, Γ_N, g)

# ── Plot ──────────────────────────────────────────────────────────────────────

p = plot_concentration(centers, σs, boundaries, ϕinf)
display(p)
readline()