
# Ns = 2048 ... Ux = -0.002015491351649, Uy = 0.149827388413482

"""
Convergence study: two unit circles, boundary-to-boundary gap = 2e-2.

COM distance d_cm = 2 + 2e-2.  Both particles oriented at θ = π/2.
Sweeps arc quadrature points N and records:
  - x-velocity of particle 1 (proxy for Marangoni force)
  - surfactant layer densities σ on all four arcs

Produces two plots:
  (1) Ux convergence vs N
  (2) σ profiles on particle 1 arcs for each N
"""

using ProjectedSphericalHarmonics
using Printf
using Plots

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

# ── Parameters ─────────────────────────────────────────────────────────────────

α_arc = π/2
g(θ)  = -cos(θ)
a, b  = 1.0, 1.0
gap   = 1e-2
d_cm  = 2.0 + gap          # COM distance; boundary-to-boundary gap = 2e-2

M     = 32                  # spectral degree for Stokes domain (fixed)
Ns_vals = [16, 32, 64, 128, 256, 512, 1024, 2048]
xcm = [-d_cm/2 + 0im, d_cm/2 + 0im]
θcm = [π/2, π/2]

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Ω = ellipse(a, b, M)

# ── Sweep ──────────────────────────────────────────────────────────────────────

Ux_vals, Uy_vals = Float64[], Float64[]
σ_runs  = []   # σ_runs[k] = σ vector (4 arcs) for Ns_vals[k]
Γ_runs  = []

for Ns in Ns_vals
    @printf("Ns = %d ... ", Ns)
    Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
    Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

    U, σ, Γ, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g, eps=1e-12)

    push!(Ux_vals, U[1][1])
    push!(Uy_vals, U[1][2])
    push!(σ_runs, σ)
    push!(Γ_runs, Γ)
    @printf("Ux = %.15f, Uy = %.15f\n", U[1][1], U[1][2])
end

# ── Plot 1: Ux convergence ─────────────────────────────────────────────────────

Ux_ref, Uy_ref = Ux_vals[end], Uy_vals[end]
err = sqrt.(abs.(Ux_vals[1:end-1] .- Ux_ref).^2 .+ abs.(Uy_vals[1:end-1] .- Uy_ref).^2)

p1 = plot(Ns_vals[1:end-1], err,
          marker=:circle, markersize=5,
          xscale=:log2, yscale=:log10,
          xlabel="N", ylabel="|U − U_ref|",
          label="error vs N=$(Ns_vals[end])",
          title="Convergence: U of particle 1  (gap = $gap)",
          legend=:topright)

# ── Plot 2: σ profiles on particle 1, D and N arcs separately ─────────────────

colors = palette(:viridis, length(Ns_vals))
p2 = plot(title="σ_D  (particle 1)", legend=:outertopright)
p3 = plot(title="σ_N  (particle 1)", legend=:outertopright)

for (k, Ns) in enumerate(Ns_vals)
    σk = σ_runs[k]
    Γk = Γ_runs[k]
    # Arcs: [D1, N1, D2, N2]  (order from surfactant_solve)
    lbl = "N=$Ns"
    plot!(p2, Γk[1].s, σk[1], color=colors[k], label=lbl)
    plot!(p3, Γk[2].s, σk[2], color=colors[k], label=lbl)
end
xlabel!(p2, "s");  ylabel!(p2, "σ")
xlabel!(p3, "s");  ylabel!(p3, "σ")

p = plot(p1, p2, p3, layout=(1,3), size=(1200, 380))
display(p)
readline()

# 64, 0.00627231
# 128, 0.00627299
# 256, 0.00626794
# 512, 0.00627012
# 1024, 0.00627021
# 2048, 0.00627023
# 4096, 0.00627001