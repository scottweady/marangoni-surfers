
# Locate the anti-parallel orbit directly from the full boundary integral solve,
# without the multipole expansion.
#
# Symmetric family: particle 1 at +(r/2) e_x with orientation p1 at angle chi
# from e_x, particle 2 at -(r/2) e_x with p2 = -p1. The pi-rotation x -> -x maps
# particle 1 onto particle 2 for any chi, which forces u2 = -u1 and w2 = w1.
#
# Translation and rotation invariance leave the two-dimensional shape space
# (r, chi). A steady orbit of constant separation is a FIXED POINT of that shape
# space, so both shape coordinates must be stationary:
#
#     rdot   = 2 * ux1                    = 0
#     chidot = w1 - Om_r = w1 - 2*uy1/r   = 0
#
# We write these as
#
#     F1 = ux1                = 0
#     F2 = uy1 - w1 * r / 2   = 0        (chidot = -2*F2/r)
#
# The paper's ansatz (5.2) fixes chi = pi/2. That is one slice of the shape
# space, and on it F1 and F2 generally do not vanish at the same r, so chi must
# be solved for rather than assumed. This script scans the (r, chi) plane,
# extracts the F1 = 0 and F2 = 0 nullclines, and reports their intersection --
# the rigidly rotating orbit, if one exists.
#
# The multipole model satisfies F1 identically for every (r, chi = pi/2), so
# within it the ansatz collapses to F2 = 0 with root r* = beta/(4 pi uS)(1-lam/2).

include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

α = π/2
g(θ) = -cos(θ)

ratios = [1/16, 1.0, 16.0]                             # a/b, fixed area ab = 1
chis   = range(-π, 0, length=8)                         # angle from e_r to p1
gaps   = exp.(range(log(5.0), log(20.0), length=8))   # edge-to-edge separation

Ns = 128
M  = 16

# ── Sweep ─────────────────────────────────────────────────────────────────────

rows = NTuple{11,Float64}[]

for ratio in ratios

    a = sqrt(ratio)
    b = 1.0 / sqrt(ratio)

    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    Ω   = ellipse(a, b, M)
    Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
    Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

    # Single-particle quantities and the far-field prediction
    U1, σ1, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
    uS = U1[1]
    β  = -(integral(real.(Γ_D.γ) .* σ1[1], Γ_D) + integral(real.(Γ_N.γ) .* σ1[2], Γ_N))
    λ  = (a^2 - b^2) / (a^2 + b^2)
    rstar = β / (4π * uS) * (1 - λ/2)

    @printf("\n=== a/b = %.4g   uS = %.6f   beta = %.6f   lam = %+.4f   r* = %.4f   2b = %.4f\n",
            ratio, uS, β, λ, rstar, 2b)

    # Nullclines: for each chi, the r at which F1 (resp. F2) changes sign
    null1 = Tuple{Float64,Float64}[]
    null2 = Tuple{Float64,Float64}[]

    for χ in chis

        F1s = Float64[]
        F2s = Float64[]
        rs  = Float64[]

        for gap in gaps

            r = 2b + gap

            xcm = [r/2 + 0im, -r/2 + 0im]
            θcm = [χ, χ + π]

            U, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)

            ux1, uy1, ω1 = U[1]
            ux2, uy2, ω2 = U[2]

            F1 = ux1
            F2 = uy1 - ω1 * r / 2

            # symmetry check: u2 = -u1 and w2 = w1
            sym = max(abs(ux2 + ux1), abs(uy2 + uy1), abs(ω2 - ω1))

            push!(rs, r); push!(F1s, F1); push!(F2s, F2)
            push!(rows, (ratio, a, b, χ, r, gap, ux1, uy1, ω1, F1, F2))

            @printf("  chi/pi = %.3f  r = %8.4f  ux1 = %11.4e  uy1 = %11.4e  w1 = %11.4e  F2 = %11.4e  sym = %.0e\n",
                    χ/π, r, ux1, uy1, ω1, F2, sym)

        end

        for (F, null) in ((F1s, null1), (F2s, null2))
            ic = findfirst(k -> F[k] * F[k+1] < 0, 1:length(F)-1)
            ic === nothing && continue
            r0, r1 = rs[ic], rs[ic+1]
            f0, f1 = F[ic], F[ic+1]
            push!(null, (χ, r0 - f0 * (r1 - r0) / (f1 - f0)))
        end

    end

    @printf("\n  F1 = 0 nullcline (chi/pi, r): %s\n",
            isempty(null1) ? "none" : join([@sprintf("(%.3f, %.3f)", c/π, r) for (c,r) in null1], "  "))
    @printf("  F2 = 0 nullcline (chi/pi, r): %s\n",
            isempty(null2) ? "none" : join([@sprintf("(%.3f, %.3f)", c/π, r) for (c,r) in null2], "  "))

end

# ── Write results ─────────────────────────────────────────────────────────────

outfile = joinpath(@__DIR__, "antiparallel_equilibrium.dat")
open(outfile, "w") do io
    println(io, "# a/b  a  b  chi  r  gap  ux1  uy1  omega1  F1  F2")
    for row in rows
        @printf(io, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", row...)
    end
end
println("\nResults written to $outfile")

# ── Plot: nullclines in the (chi, r) shape space ──────────────────────────────

plots = map(ratios) do ratio
    sel = [row for row in rows if row[1] == ratio]
    p = plot(xlabel = "χ/π", ylabel = "separation, r",
             title = @sprintf("a/b = %.4g", ratio), legend = :topright)
    for (col, lab, c) in ((10, "ṙ = 0", :crimson), (11, "χ̇ = 0", :navy))
        pts = Tuple{Float64,Float64}[]
        for χ in chis
            sub = [row for row in sel if row[4] == χ]
            F   = [row[col] for row in sub]
            rs  = [row[5]   for row in sub]
            ic  = findfirst(k -> F[k] * F[k+1] < 0, 1:length(F)-1)
            ic === nothing && continue
            push!(pts, (χ/π, rs[ic] - F[ic] * (rs[ic+1] - rs[ic]) / (F[ic+1] - F[ic])))
        end
        isempty(pts) || plot!(p, first.(pts), last.(pts);
                              marker = :circle, markersize = 3, color = c, label = lab)
    end
    p
end

display(plot(plots...; layout = (1, length(ratios)), size = (1100, 380)))
readline()
