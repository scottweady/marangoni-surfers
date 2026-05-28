"""
Convergence study: single unit circle, oriented at θ = π/2.

Sweeps arc quadrature points Ns and records:
  - translational velocity U = (Ux, Uy)
  - angular velocity ω
  - surfactant layer densities σ on the two arcs

Produces two plots:
  (1) |U − U_ref| convergence vs Ns
  (2) σ profiles on D and N arcs for each Ns
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

M       = 16
Ns_vals = [16, 32, 64, 128, 256, 512, 1024, 2048]

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Ω = ellipse(a, b, M)

# ── Sweep ──────────────────────────────────────────────────────────────────────

Ux_vals = Float64[]
Uy_vals = Float64[]
ω_vals  = Float64[]
σ_runs  = []
Γ_runs  = []

for Ns in Ns_vals
    @printf("Ns = %4d ... ", Ns)
    Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
    Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

    U, ω, α_D, α_N, Γ_D, Γ_N, _ = marangoni_mobility(Ω, Γ_D, Γ_N, g)

    push!(Ux_vals, U[1])
    push!(Uy_vals, U[2])
    push!(ω_vals,  ω)
    push!(σ_runs, (α_D, α_N))
    push!(Γ_runs, (Γ_D, Γ_N))
    @printf("Ux = %.8f, Uy = %.8f, ω = %.8f\n", U[1], U[2], ω)
end

# ── Plot 1: convergence of |U − U_ref| ────────────────────────────────────────

Ux_ref = Ux_vals[end]
Uy_ref = Uy_vals[end]
err = sqrt.(abs.(Ux_vals[1:end-1] .- Ux_ref).^2 .+ abs.(Uy_vals[1:end-1] .- Uy_ref).^2)

p1 = plot(Ns_vals[1:end-1], err,
          marker=:circle, markersize=5,
          xscale=:log2, yscale=:log10,
          xlabel="Ns", ylabel="|U − U_ref|",
          label="error vs Ns=$(Ns_vals[end])",
          title="Convergence: single circle  (α = π/2)",
          legend=:topright)

# ── Plot 2: σ profiles on D and N arcs ────────────────────────────────────────

colors = palette(:viridis, length(Ns_vals))
p2 = plot(title="σ_D", legend=:outertopright)
p3 = plot(title="σ_N", legend=:outertopright)

for (k, Ns) in enumerate(Ns_vals)
    α_D, α_N = σ_runs[k]
    Γ_D, Γ_N = Γ_runs[k]
    lbl = "Ns=$Ns"
    plot!(p2, Γ_D.s, α_D, color=colors[k], label=lbl)
    plot!(p3, Γ_N.s, α_N, color=colors[k], label=lbl)
end
xlabel!(p2, "s"); ylabel!(p2, "σ")
xlabel!(p3, "s"); ylabel!(p3, "σ")

p = plot(p1, p2, p3, layout=(1,3), size=(1200, 380))
display(p)
readline()
