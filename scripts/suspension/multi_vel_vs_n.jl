# Speed of the center particle in a line of equally spaced surfers with a common
# orientation, normalized by the isolated-particle speed, sweeping the (odd) number of
# particles over `Nvals` and the spacing over `dvals`.
#
# Output: center_speed_vs_N.dat
#   N  followed by |U|/U0 for each spacing

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

Ns   = 128
M    = 16
s₀   = π/2
g(s) = -cos(s)

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

dvals = [5.0, 10.0, 15.0, 20.0]   # center-to-center spacings between adjacent particles
θ0   = π/2   # orientation angle of all particles

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Compute single-particle speed for normalization ---------------------------

U0_vec, _, _, _ = marangoni_solve([0.0+0im], [θ0], Ω, Γ_D, Γ_N, g)
U0 = hypot(U0_vec[1][1], U0_vec[1][2])
@printf("Isolated particle speed U0 = %.6f\n", U0)

# -- Sweep over odd N and spacing d --------------------------------------------

Nvals    = 1:2:13
U_center = zeros(length(Nvals), length(dvals))

for (j, d) in enumerate(dvals), (i, N) in enumerate(Nvals)
    ic  = (N + 1) ÷ 2
    xcm = ComplexF64[d * (k - (N + 1) / 2) for k in 1:N]
    θcm = fill(θ0, N)

    U, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
    uc = hypot(U[ic][1], U[ic][2])
    U_center[i, j] = uc
    @printf("d = %4.1f  N = %2d:  |U| = %.6f  (|U|/U0 = %.4f)\n", d, N, uc, uc / U0)
end

# -- Write results -------------------------------------------------------------

outfile = joinpath(@__DIR__, "center_speed_vs_N.dat")
open(outfile, "w") do io
    @printf(io, "# N  |U|/U0 for d = %s  (theta0=%.4f  U0=%.10f)\n", join(dvals, " "), θ0, U0)
    for (i, N) in enumerate(Nvals)
        @printf(io, "%d  %s\n", N, join([@sprintf("%.15e", uc / U0) for uc in U_center[i, :]], "  "))
    end
end
println("Written to $outfile")

# -- Plot ----------------------------------------------------------------------

p = plot(collect(Nvals), U_center ./ U0;
         xlabel  = "N",
         ylabel  = "|U| / U₀",
         title   = "Center particle speed (θ = $(round(θ0/π, digits=2))π)",
         label   = permutedims(["d = $d" for d in dvals]),
         marker  = :circle,
         lw      = 2,
         ylims   = (0, nothing))
hline!(p, [1.0]; ls=:dash, color=:gray, label="isolated")

display(p)
readline()
