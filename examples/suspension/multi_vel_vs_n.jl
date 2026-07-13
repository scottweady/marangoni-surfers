

include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Ns   = 1024
M    = 32
α    = π/2
a, b = 1.0, 1.0
g(θ) = -cos(θ)
d    = 20.0   # center-to-center spacing between adjacent particles
θ0   = 0.0   # orientation angle of all particles

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Ω   = ellipse(a, b, M)
Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

# ── Compute single-particle speed for normalization ───────────────────────────

U0_vec, _, _, _ = marangoni_solve([0.0+0im], [θ0], Ω, Γ_D, Γ_N, g)
U0 = hypot(U0_vec[1][1], U0_vec[1][2])
@printf("Isolated particle speed U0 = %.6f\n", U0)

# ── Sweep over odd N ──────────────────────────────────────────────────────────

Nvals    = 1:2:13
U_center = Float64[]

for N in Nvals
    ic  = (N + 1) ÷ 2
    xcm = ComplexF64[d * (k - (N + 1) / 2) for k in 1:N]
    θcm = fill(θ0, N)

    U, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
    uc = hypot(U[ic][1], U[ic][2])
    push!(U_center, uc)
    @printf("N = %2d:  |U| = %.6f  (|U|/U0 = %.4f)\n", N, uc, uc / U0)
end

# ── Write results ─────────────────────────────────────────────────────────────

outfile = joinpath(@__DIR__, "center_speed_vs_N.dat")
open(outfile, "w") do io
    @printf(io, "# N  |U|/U0  (d=%.1f  theta0=%.4f  U0=%.10f)\n", d, θ0, U0)
    for (N, uc) in zip(Nvals, U_center)
        @printf(io, "%d  %.15e\n", N, uc / U0)
    end
end
println("Written to $outfile")

# ── Plot ──────────────────────────────────────────────────────────────────────

p = plot(collect(Nvals), U_center ./ U0;
         xlabel  = "N",
         ylabel  = "|U| / U₀",
         title   = "Center particle speed (d = $d, θ = $(round(θ0/π, digits=2))π)",
         label   = "|U| / U₀",
         marker  = :circle,
         lw      = 2,
         ylims   = (0, nothing))
hline!(p, [1.0]; ls=:dash, color=:gray, label="isolated")

display(p)
readline()
