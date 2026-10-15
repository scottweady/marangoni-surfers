# Single elliptical surfer: solves for the rigid-body velocity, prints the multipole
# moments, and plots the surfactant concentration.

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

aspect_ratio = 2.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

# -- Geometry ------------------------------------------------------------------

Ns = 64
M  = 16

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Solve ---------------------------------------------------------------------

U, σ, Γ, ϕinf = marangoni_solve(Ω, Γ_D, Γ_N, g)

# -- Multipoles ----------------------------------------------------------------

multipoles = marangoni_multipoles(Ω, Γ_D, Γ_N, g)

@printf("velocity            = (%.6f, %.6f, %.6f)\n", multipoles.velocity...)
@printf("force monopole      = (%.6f, %.6f)\n", multipoles.force_monopole...)
@printf("force dipole        = [%.6f %.6f; %.6f %.6f]\n", multipoles.force_dipole...)
@printf("surfactant monopole = (%.6f, %.6f)\n", multipoles.surf_monopole...)
@printf("surfactant dipole   = [%.6f %.6f; %.6f %.6f]\n", multipoles.surf_dipole...)

# -- Plot concentration --------------------------------------------------------

p = plot_concentration([0.0im], σ, ϕinf, Γ)
display(p)
readline()
