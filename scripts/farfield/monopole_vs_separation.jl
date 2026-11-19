# Total surfactant density ∫ σ dΓ on each particle of a pair versus center-to-center
# separation, with both particles oriented along or perpendicular to the line of centers.
#
# Output: output/monopole_vs_separation.dat
#   d  I1  I2 (aligned)  I1  I2 (perpendicular)

using MarangoniSurfers
using Printf

# -- Parameters ----------------------------------------------------------------

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

# -- Geometry ------------------------------------------------------------------

Ns = 64

Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)

separations = 10 .^ range(log10(3.0), 6.0, length=200)

# -- Sweep: aligned ------------------------------------------------------------

θcm_0 = [0.0, 0.0]
integrals_0 = map(separations) do d
    xcm = [-d/2 + 0im, d/2 + 0im]
    σ, _, Γ = surfactant_solve(xcm, θcm_0, Γ_D, Γ_N, g)
    I1 = integral(σ[1], Γ[1]) + integral(σ[2], Γ[2])
    I2 = integral(σ[3], Γ[3]) + integral(σ[4], Γ[4])
    (I1, I2)
end

# -- Sweep: perpendicular ------------------------------------------------------

θcm_pi2 = [π/2, π/2]
integrals_pi2 = map(separations) do d
    xcm = [-d/2 + 0im, d/2 + 0im]
    σ, _, Γ = surfactant_solve(xcm, θcm_pi2, Γ_D, Γ_N, g)
    I1 = integral(σ[1], Γ[1]) + integral(σ[2], Γ[2])
    I2 = integral(σ[3], Γ[3]) + integral(σ[4], Γ[4])
    (I1, I2)
end

# -- Write results -------------------------------------------------------------

outdir = script_param(:outdir, joinpath(@__DIR__, "output"))
mkpath(outdir)
outfile = joinpath(outdir, "monopole_vs_separation.dat")

open(outfile, "w") do io
    for (i, d) in enumerate(separations)
        I1_0,   I2_0   = integrals_0[i]
        I1_pi2, I2_pi2 = integrals_pi2[i]
        @printf(io, "%.10e  %.10e  %.10e  %.10e  %.10e\n",
                d, I1_0, I2_0, I1_pi2, I2_pi2)
    end
end

println("Wrote ", relpath(outfile, @__DIR__))
