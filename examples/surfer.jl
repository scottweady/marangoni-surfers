
using ProjectedSphericalHarmonics
using Printf

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

a, b = 1.0, 1.0
α = π/2
g(θ) = -cos(θ)

# ── Geometry ──────────────────────────────────────────────────────────────────

Ns = 256
M = 32

γ(s)   =  a * cos.(s) + im * b * sin.(s)
γs(s)  = -a * sin.(s) + im * b * cos.(s)
γss(s) = -a * cos.(s) - im * b * sin.(s)

Ω   = ellipse(a, b, M)
Γ_D = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α, α])
Γ = [Γ_D, Γ_N]

# ── Solve ────────────────────────────────────────────────────────────────────
multipoles = marangoni_multipoles(Ω, Γ, g)

@printf("velocity            = (%.6f, %.6f, %.6f)\n", multipoles.velocity...)
@printf("force monopole      = (%.6f, %.6f)\n", multipoles.force_monopole...)
@printf("force dipole        = [%.6f %.6f; %.6f %.6f]\n", multipoles.force_dipole...)
@printf("surfactant monopole = (%.6f, %.6f)\n", multipoles.surf_monopole...)
@printf("surfactant dipole   = [%.6f %.6f; %.6f %.6f]\n", multipoles.surf_dipole...)
