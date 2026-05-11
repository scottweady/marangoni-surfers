
using Printf
using ProjectedSphericalHarmonics

include("bim/boundary.jl")
include("bim/operators.jl")
include("bim/solvers.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

α    = π/2
g(θ) = cos(2θ)
a, b = 1.0, 1.0

# ── Run ───────────────────────────────────────────────────────────────────────

Ns = 256
Nq = 64
U, ω, α_D, α_N, Γ_D, Γ_N, Ω = marangoni_mobility(Ns, Nq, a, b, α, g)
mp = marangoni_multipoles(U, ω, α_D, α_N, Γ_D, Γ_N, Ω)

@printf("U                   = (%.6f, %.6f)\n", U[1], U[2])
@printf("force monopole      = (%.6f, %.6f)\n", mp.force_monopole...)
@printf("force dipole        = [%.6f %.6f; %.6f %.6f]\n", mp.force_dipole...)
@printf("surfactant monopole = (%.6f, %.6f)\n", mp.surf_monopole...)
@printf("surfactant dipole   = [%.6f %.6f; %.6f %.6f]\n", mp.surf_dipole...)
