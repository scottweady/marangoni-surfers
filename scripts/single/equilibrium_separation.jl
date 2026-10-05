# Equilibrium separation of an anti-parallel pair, predicted by the far-field
# multipole expansion:
#
#     r = β / (4π uˢ) * (1 - λ/2),     λ = (a² - b²) / (a² + b²)
#
# where uˢ is the single-particle propulsion speed and β is the surfactant
# dipole strength, ∫_Γ x σ dΓ = -β, both computed from the single-particle
# boundary integral solve. The sweep is over aspect ratio a/b at fixed area
# (ab = 1), consistent with the length scale ℓ_c = (|Ω|/π)^(1/2).
#
# In this configuration the pair is offset perpendicular to the propulsion
# axis, so the particles touch at r = 2b; r < 2b is unphysical.

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀ = π/2
g(s) = -cos(s)

aspect_ratios = exp.(range(log(1.0/16.0), log(16.0), length=25))   # log-spaced, includes 1

# -- Geometry ------------------------------------------------------------------

Ns = 64
M = 16

# -- Sweep ---------------------------------------------------------------------

uS_vals = Float64[]
β_vals  = Float64[]
λ_vals  = Float64[]
r_vals  = Float64[]

for aspect_ratio in aspect_ratios

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

    U, σ, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
    σ_D, σ_N = σ[1], σ[2]

    uS = U[1]

    # Surfactant dipole: ∫_Γ x σ dΓ = -β, with x the physical coordinate
    β = -(integral(real.(Γ_D.γ) .* σ_D, Γ_D) + integral(real.(Γ_N.γ) .* σ_N, Γ_N))

    λ = (a^2 - b^2) / (a^2 + b^2)
    r = β / (4π * uS) * (1 - λ/2)

    push!(uS_vals, uS)
    push!(β_vals,  β)
    push!(λ_vals,  λ)
    push!(r_vals,  r)

    @printf("a/b = %8.4f:  uS = %.6f  β = %.6f  λ = %+.6f  r = %.6f  (2b = %.4f)\n",
            aspect_ratio, uS, β, λ, r, 2b)

end

# -- Write results -------------------------------------------------------------

outfile = joinpath(@__DIR__, "..", "..", "data", "single", "equilibrium_separation.dat")
open(outfile, "w") do io
    println(io, "# a/b  a  b  uS  beta  lambda  r  r_contact")
    for (i, aspect_ratio) in enumerate(aspect_ratios)
        a = sqrt(aspect_ratio)
        b = 1.0 / a
        @printf(io, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
                aspect_ratio, a, b, uS_vals[i], β_vals[i], λ_vals[i], r_vals[i], 2b)
    end
end
println("Results written to $outfile")

# -- Plot ----------------------------------------------------------------------

p = plot(aspect_ratios, r_vals;
         xlabel = "aspect ratio, a/b",
         ylabel = "equilibrium separation, r",
         label  = "β/(4π uˢ)(1 - λ/2)",
         xscale = :log10,
         marker = :circle,
         legend = :topright)

plot!(p, aspect_ratios, 2 ./ sqrt.(aspect_ratios); label = "contact, 2b", linestyle = :dash, color = :black)

display(p)
readline()
