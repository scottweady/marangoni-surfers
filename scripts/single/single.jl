

include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

a = sqrt(1.0/16.0)
b = 1.0/a
α = π/2
g(θ) = -cos(θ)

# ── Geometry ──────────────────────────────────────────────────────────────────

Ns = 1024
M = 32

γ(s)   =  a * cos(s) + im * b * sin(s)
γs(s)  = -a * sin(s) + im * b * cos(s)
γss(s) = -a * cos(s) - im * b * sin(s)

Ω   = ellipse(a, b, M)
Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

# ── Solve ─────────────────────────────────────────────────────────────────────

U, σ, Γ, ϕinf = marangoni_solve(Ω, Γ_D, Γ_N, g)

# ── Multipoles ────────────────────────────────────────────────────────────────

multipoles = marangoni_multipoles(Ω, Γ_D, Γ_N, g)

@printf("velocity            = (%.6f, %.6f, %.6f)\n", multipoles.velocity...)
@printf("force monopole      = (%.6f, %.6f)\n", multipoles.force_monopole...)
@printf("force dipole        = [%.6f %.6f; %.6f %.6f]\n", multipoles.force_dipole...)
@printf("surfactant monopole = (%.6f, %.6f)\n", multipoles.surf_monopole...)
@printf("surfactant dipole   = [%.6f %.6f; %.6f %.6f]\n", multipoles.surf_dipole...)

# ── Plot concentration ────────────────────────────────────────────────────────

p = plot_concentration([0.0im], σ, ϕinf, Γ)
display(p)
readline()
