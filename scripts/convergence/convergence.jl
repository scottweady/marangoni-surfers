# Convergence of the single-particle swimming speed against the exact solution for an
# elliptical surfer, sweeping aspect ratio over `aspect_ratios` at fixed area. The error is
# measured separately in the Stokes expansion order (`Mspan`, with the boundary resolved
# at `Ns_ref`) and in the boundary resolution (`Nspan`, with expansion order `M_ref`).
#
# Output: output/convergence_<expansion|boundary>.dat
#   Header: Ns_ref, M_ref, aspect ratios, exact speeds
#   Rows:   M (or N)  followed by the error for each aspect ratio

include("../../src/marangoni.jl")

using SpecialFunctions: ellipk, ellipe

# -- Parameters ----------------------------------------------------------------

s₀ = π/2
g(s) = -cos(s)

aspect_ratios = [0.25, 0.5, 1.0, 2.0, 4.0]

Ns_ref = 2048
M_ref  = 32

Mspan = 2:2:32
Nspan = Int.(2.0.^(2:10)) # boundary resolution sweep

# -- Exact solution ------------------------------------------------------------

# Resistance formulation: U = Δm/(2πab) ∫ x σ ds. The Joukowski map takes the disk
# surfactant problem to the ellipse one (g is written in the ellipse parameter s), scaling
# the dipole moment ∫ x σ ds by (a+b)/2. The disk value -π(1 + cos s₀)²/2 follows from
# Lauga & Davis (2012, eqs. 3.39-3.40), U = π(1 + cos s₀)²/32 with RT = μ = a = 1.
function u_exact(a, b)
    Δm = if abs(a - b) < 1e-8 * (a + b)
        -π * min(a, b) / 8
    elseif a > b
        ε² = 1 - (b/a)^2
        b/(2ε²) * (ellipk(ε²) - ellipe(ε²)) - b/2 * ellipk(ε²)
    else
        ε² = 1 - (a/b)^2
        a/(2ε²) * (ellipe(ε²) - (1 - ε²) * ellipk(ε²)) - a/2 * ellipk(ε²)
    end
    return -(a + b) * (1 + cos(s₀))^2 * Δm / (8 * a * b)
end

# -- Main loop -----------------------------------------------------------------

err_M = zeros(length(Mspan), length(aspect_ratios))
err_N = zeros(length(Nspan), length(aspect_ratios))
Uref  = zeros(length(aspect_ratios))

for (j, aspect_ratio) in enumerate(aspect_ratios)

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    u_ref = u_exact(a, b)
    Uref[j] = u_ref

    @printf("\na/b = %.4f (a = %.6f, b = %.6f):  u* = %.15e\n", aspect_ratio, a, b, u_ref)

    # (i) Expansion order M, boundary fixed at a well-resolved Ns_ref
    Γ_D_ref, Γ_N_ref = surfer_boundary(a, b, s₀, Ns_ref)
    println("  expansion order M (N = $Ns_ref):")
    for (i, M) in enumerate(Mspan)
        Ω = ellipse(a, b, M)
        U, _, _, _ = marangoni_solve(Ω, Γ_D_ref, Γ_N_ref, g)
        err_M[i, j] = abs(U[1] - u_ref)
        @printf("    M = %3d:  u = %+.15e  err = %.6e\n", M, U[1], err_M[i, j])
    end

    # (ii) Boundary resolution N, expansion fixed at M_ref
    Ω_ref = ellipse(a, b, M_ref)
    println("  boundary resolution N (M = $M_ref):")
    for (i, Ns) in enumerate(Nspan)
        Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)
        U, _, _, _ = marangoni_solve(Ω_ref, Γ_D, Γ_N, g)
        err_N[i, j] = abs(U[1] - u_ref)
        @printf("    N = %3d:  u = %+.15e  err = %.6e\n", Ns, U[1], err_N[i, j])
    end

end

# -- Write results -------------------------------------------------------------

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

function write_errors(filename, label, span, err)
    outfile = joinpath(outdir, filename)
    open(outfile, "w") do io
        println(io, "# single-particle velocity error vs exact solution")
        @printf(io, "# Ns_ref = %d, M_ref = %d\n", Ns_ref, M_ref)
        println(io, "# aspect ratios a/b: ", join(aspect_ratios, " "))
        println(io, "# exact u*:          ", join(Uref, " "))
        println(io, "# columns: $label  err(a/b = ", join(aspect_ratios, ")  err(a/b = "), ")")
        for (i, n) in enumerate(span)
            @printf(io, "%d", n)
            for j in eachindex(aspect_ratios)
                @printf(io, " %.15e", err[i, j])
            end
            println(io)
        end
    end
    println("Results written to $outfile")
end

write_errors("convergence_expansion.dat", "M", Mspan, err_M)
write_errors("convergence_boundary.dat", "N", Nspan, err_N)

# -- Plot results --------------------------------------------------------------

labels = reshape(["a/b = $(r)" for r in aspect_ratios], 1, :)

p1 = plot(collect(Mspan), err_M;
          xlabel = "expansion order, M",
          ylabel = "error, |u_cm - u*|",
          yscale = :log10,
          marker = :circle,
          label  = labels,
          legend = :topright)

p2 = plot(collect(Nspan), err_N;
          xlabel = "boundary resolution, N",
          ylabel = "error, |u_cm - u*|",
          yscale = :log10,
          marker = :circle,
          label  = labels,
          legend = :topright)

p = plot(p1, p2; layout=(1, 2), size=(900, 380))
display(p)
readline()
