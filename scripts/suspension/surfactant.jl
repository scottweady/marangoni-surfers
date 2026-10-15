# Surfactant concentration for a static configuration of surfers placed along a spiral,
# each oriented tangent to it.

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

# -- Geometry ------------------------------------------------------------------

Ns   = 64

Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)

# -- Configuration -------------------------------------------------------------

Nparticles   = 2
r0, dr   = 8.0, 4.0 # inner radius, radial growth per turn
θ_spiral = range(0, 4π, length=Nparticles+1)[1:end-1]
r = r0 .+ dr .* θ_spiral ./ (2π)
xcm  = r .* exp.(im .* θ_spiral)

# tangent angle of the spiral at each point
drdθ = dr / (2π)
θcm = atan.(drdθ .* sin.(θ_spiral) .+ r .* cos.(θ_spiral),
            drdθ .* cos.(θ_spiral) .- r .* sin.(θ_spiral))

# -- Solve ---------------------------------------------------------------------

σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

p = plot_concentration(xcm, σ, ϕinf, Γ)
display(p)
readline()