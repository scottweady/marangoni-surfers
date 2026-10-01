

include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

a, b  = 1.0, 1.0
α = π/2
g(θ)  = -cos(θ)


# ── Geometry ──────────────────────────────────────────────────────────────────

Ns = 1024

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

separations = range(2.0, 100.0, length=100)[2:end]   # skip d=0 to avoid singularity

# ── Sweep: θ = 0 ──────────────────────────────────────────────────────────────

θcm_0 = [0.0, 0.0]
integrals_0 = map(separations) do d
  xcm = [-d/2 + 0im, d/2 + 0im]
  σ, _, Γ = surfactant_solve(xcm, θcm_0, Γ_D, Γ_N, g)
  I1 = integral(σ[1], Γ[1]) + integral(σ[2], Γ[2])
  I2 = integral(σ[3], Γ[3]) + integral(σ[4], Γ[4])
  (I1, I2)
end

# ── Sweep: θ = π/2 ────────────────────────────────────────────────────────────

θcm_pi2 = [π/2, π/2]
integrals_pi2 = map(separations) do d
  xcm = [-d/2 + 0im, d/2 + 0im]
  σ, _, Γ = surfactant_solve(xcm, θcm_pi2, Γ_D, Γ_N, g)
  I1 = integral(σ[1], Γ[1]) + integral(σ[2], Γ[2])
  I2 = integral(σ[3], Γ[3]) + integral(σ[4], Γ[4])
  (I1, I2)
end

# ── Write .dat file ───────────────────────────────────────────────────────────
# Columns: d  I1_th0  I2_th0  I1_thpi2  I2_thpi2

out_dir = joinpath(@__DIR__, "output")
mkpath(out_dir)

open(joinpath(out_dir, "sigma_vs_separation.dat"), "w") do io
    for (i, d) in enumerate(separations)
        I1_0,   I2_0   = integrals_0[i]
        I1_pi2, I2_pi2 = integrals_pi2[i]
        @printf(io, "%.6f  %.10f  %.10f  %.10f  %.10f\n",
                d, I1_0, I2_0, I1_pi2, I2_pi2)
    end
end

println("\nWrote sigma_vs_separation.dat")
