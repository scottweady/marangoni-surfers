
include("../../src/marangoni.jl")

using SpecialFunctions: ellipk, ellipe

##### Parameters #####

α = π/2
g(θ) = exp.(-cos(θ))

ratios = [0.25, 0.5, 1.0, 2.0, 4.0]  # a/b, at fixed area ab = 1

Ns_ref = 2048
M_ref  = 32

Mspan = 2:2:32
Nspan = Int.(2.0.^(2:10)) # boundary resolution sweep

##### Solution from resistance formulation #####
function ellipse_mobility(a, b)
    Δm = if abs(a - b) < 1e-8 * (a + b)
        -π * min(a, b) / 8
    elseif a > b
        ε² = 1 - (b/a)^2
        b/(2ε²) * (ellipk(ε²) - ellipe(ε²)) - b/2 * ellipk(ε²)
    else
        ε² = 1 - (a/b)^2
        a/(2ε²) * (ellipe(ε²) - (1 - ε²) * ellipk(ε²)) - a/2 * ellipk(ε²)
    end
    return Δm / (2π * a * b)
end

semi_analytical(σ, Γ, a, b) =
    ellipse_mobility(a, b) * (integral(real.(Γ[1].γ) .* σ[1], Γ[1]) +
                              integral(real.(Γ[2].γ) .* σ[2], Γ[2]))

##### Main loop #####

err_M = zeros(length(Mspan), length(ratios))
err_N = zeros(length(Nspan), length(ratios))
Uref  = zeros(length(ratios))

for (j, ratio) in enumerate(ratios)

    a = sqrt(ratio)
    b = 1.0 / sqrt(ratio)

    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    # Reference from a highly resolved surfactant solve
    Γ_D_ref = openarc(γ, γs, γss, Ns_ref; dom=[α, 2π-α], quad=:open_weight)
    Γ_N_ref = openarc(γ, γs, γss, Ns_ref; dom=[-α, α], quad=:open)

    σ_ref, _, Γ_ref = surfactant_solve([0.0im], [0.0], Γ_D_ref, Γ_N_ref, g)
    u_ref = semi_analytical(σ_ref, Γ_ref, a, b)
    Uref[j] = u_ref

    @printf("\na/b = %.4f (a = %.6f, b = %.6f):  u* = %.15e\n", ratio, a, b, u_ref)

    # (i) Expansion order M, boundary fixed at Ns_ref
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
        Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α], quad=:open_weight)
        Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α], quad=:open)
        U, _, _, _ = marangoni_solve(Ω_ref, Γ_D, Γ_N, g)
        err_N[i, j] = abs(U[1] - u_ref)
        @printf("    N = %3d:  u = %+.15e  err = %.6e\n", Ns, U[1], err_N[i, j])
    end

end

##### Write results #####

outdir = joinpath(@__DIR__, "../../data/convergence")
mkpath(outdir)

function write_errors(filename, label, span, err)
    outfile = joinpath(outdir, filename)
    open(outfile, "w") do io
        println(io, "# single-particle velocity error vs semi-analytical solution")
        @printf(io, "# Ns_ref = %d, M_ref = %d\n", Ns_ref, M_ref)
        println(io, "# aspect ratios a/b: ", join(ratios, " "))
        println(io, "# reference u*:      ", join(Uref, " "))
        println(io, "# columns: $label  err(a/b = ", join(ratios, ")  err(a/b = "), ")")
        for (i, n) in enumerate(span)
            @printf(io, "%d", n)
            for j in eachindex(ratios)
                @printf(io, " %.15e", err[i, j])
            end
            println(io)
        end
    end
    println("Results written to $outfile")
end

write_errors("convergence_expansion.dat", "M", Mspan, err_M)
write_errors("convergence_boundary.dat", "N", Nspan, err_N)

##### Plot results #####

labels = reshape(["a/b = $(ratio)" for ratio in ratios], 1, :)

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
