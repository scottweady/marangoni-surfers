
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

M  = 64    # spectral degree for Stokes domain
Ns = 1024   # quadrature nodes per arc

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Ω   = ellipse(a, b, M)
Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Sweep over separation distance ─────────────────────────────────────────────

gaps     = exp10.(range(log10(0.01), log10(1.0), length=10))
Ux_vals  = Float64[]

for gap in gaps
    d_cm = 2.0 + gap   # center-to-center distance; radius = 1 so gap = d_cm - 2
    xcm  = [-d_cm/2 + 0im, d_cm/2 + 0im]
    θcm  = [π/2, π/2]

    @printf("gap = %.4f ... ", gap)
    U, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=1e-10)
    Ux = U[1][1]
    push!(Ux_vals, Ux)
    @printf("Ux = %.10f\n", Ux)
end

# ── Plot ───────────────────────────────────────────────────────────────────────

p = plot(gaps, Ux_vals,
         marker=:circle, markersize=4,
         xscale=:log10,
         xlabel="gap (boundary-to-boundary)",
         ylabel="Ux (particle 1)",
         title="Marangoni pair: x-velocity vs separation (Ns=$Ns)",
         legend=false)

display(p)
readline()
