
include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

a, b = 1.0, 1.0
α = π/2
g(θ) = -cos(θ)

# ── Geometry ──────────────────────────────────────────────────────────────────

Ns   = 32

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

# ── Configuration ──────────────────────────────────────────────────────────────

Nparticles   = 2
r0, dr   = 8.0, 4.0 # inner radius, radial growth per turn
θ_spiral = range(0, 4π, length=Nparticles+1)[1:end-1]
r = r0 .+ dr .* θ_spiral ./ (2π)
xcm  = r .* exp.(im .* θ_spiral)

# tangent angle of the spiral at each point
drdθ = dr / (2π)
θcm = atan.(drdθ .* sin.(θ_spiral) .+ r .* cos.(θ_spiral),
            drdθ .* cos.(θ_spiral) .- r .* sin.(θ_spiral))

# ── Solve ──────────────────────────────────────────────────────────────────────

σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

p = plot_concentration(xcm, σ, ϕinf, Γ)
display(p)
readline()