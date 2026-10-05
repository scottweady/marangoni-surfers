# Speed of the center particle in a line of equally spaced surfers with a common
# orientation, normalized by the isolated-particle speed, sweeping the (odd) number of
# particles over `Nvals` and the spacing over `dvals`.
#
# Output: output/center_speed_vs_N.dat
#   N  followed by |u_c|/U for each spacing

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

Ns   = 128
M    = 32
s₀   = π/2
g(s) = -cos(s)

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

dvals = [5.0, 10.0, 15.0, 20.0]   # center-to-center spacings between adjacent particles
θ0   = π/2   # orientation angle of all particles

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Compute single-particle speed for normalization ---------------------------

U_single, = marangoni_solve([0.0+0im], [θ0], Ω, Γ_D, Γ_N, g)
U = hypot(U_single[1][1], U_single[1][2])
@printf("Isolated particle speed U = %.6f\n", U)

# -- Sweep over odd N and spacing d --------------------------------------------

Nvals    = 1:2:13
U_center = zeros(length(Nvals), length(dvals))

for (j, d) in enumerate(dvals), (i, N) in enumerate(Nvals)
    ic  = (N + 1) ÷ 2
    xcm = ComplexF64[d * (k - (N + 1) / 2) for k in 1:N]
    θcm = fill(θ0, N)

    Ucm, = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
    uc = hypot(Ucm[ic][1], Ucm[ic][2])
    U_center[i, j] = uc
    @printf("d = %4.1f  N = %2d:  |u_c| = %.6f  (|u_c|/U = %.4f)\n", d, N, uc, uc / U)
end

# -- Write results -------------------------------------------------------------

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
outfile = joinpath(outdir, "center_speed_vs_N.dat")
open(outfile, "w") do io
    @printf(io, "# N  |u_c|/U for d = %s  (theta0=%.4f  U=%.10f)\n", join(dvals, " "), θ0, U)
    for (i, N) in enumerate(Nvals)
        @printf(io, "%d  %s\n", N, join([@sprintf("%.15e", uc / U) for uc in U_center[i, :]], "  "))
    end
end
println("Written to $outfile")

# -- Plot ----------------------------------------------------------------------

p = plot(collect(Nvals), U_center ./ U;
         xlabel  = "N",
         ylabel  = "|u_c| / U",
         title   = "Center particle speed (θ = $(round(θ0/π, digits=2))π)",
         label   = permutedims(["d = $d" for d in dvals]),
         marker  = :circle,
         lw      = 2,
         ylims   = (0, nothing))
hline!(p, [1.0]; ls=:dash, color=:gray, label="isolated")

display(p)
readline()
