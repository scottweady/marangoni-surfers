

include("../../src/marangoni.jl")

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
Γ_D = openarc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Sweep over separation distance ─────────────────────────────────────────────

gaps    = exp.(range(log(0.01), log(10.0), length=20))
ucm_vals = Float64[]
ω_vals  = Float64[]

for gap in gaps
    d_cm = 2.0 + gap   # center-to-center distance; radius = 1 so gap = d_cm - 2
    xcm  = [-d_cm/2 + 0im, d_cm/2 + 0im]
    θcm  = [π/2, π/2]

    @printf("gap = %.4f ... ", gap)
    U, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
    ucm = U[1][1]
    ω  = U[1][3]
    push!(ucm_vals, ucm)
    push!(ω_vals,  ω)
    @printf("ucm = %.10f  ω = %.10f\n", ucm, ω)
end

# ── Write results to file ──────────────────────────────────────────────────────

outfile = joinpath(@__DIR__, "pair_velocity_vs_separation.dat")
open(outfile, "w") do io
    println(io, "# gap ucm omega")
    for (g, u, w) in zip(gaps, ucm_vals, ω_vals)
        @printf(io, "%.15e %.15e %.15e\n", g, u, w)
    end
end
println("Results written to $outfile")

# ── Plot ───────────────────────────────────────────────────────────────────────

p = plot(gaps, ucm_vals,
         marker=:circle, markersize=4,
         xscale=:log10,
         xlabel="gap (boundary-to-boundary)",
         ylabel="ucm (particle 1)",
         title="Marangoni pair: x-velocity vs separation (Ns=$Ns)",
         legend=false)

display(p)
readline()
