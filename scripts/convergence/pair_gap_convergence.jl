# Self-convergence of the rigid-body velocities of an interacting pair of elliptical
# surfers as the boundaries approach contact. The pair sits side by side with the body
# x axes (semi-axis a) perpendicular to the line of centers, parallel or antiparallel
# (`orient`), as in scripts/pair/pair.jl, so the edge gap is the center separation minus
# 2b. Each configuration
# is an (aspect ratio, edge gap) pair at fixed area, taken over `aspect_ratios` × `gaps`.
# The error is the largest difference over both particles and all components (ux, uy, ω)
# from a reference solve at (`M_ref`, `N_ref`), measured separately in the Stokes
# expansion order (`Mspan`, with the boundary at `N_ref`) and in the boundary resolution
# (`Nspan`, with the expansion at `M_ref`).
#
# Output: output/pair_gap/pair_gap_<expansion|boundary>.dat
#   Header: N_ref, M_ref, tol, orientation, aspect ratios, edge gaps
#   Rows:   M (or N)  followed by the error for each configuration

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀   = π/2
g(s) = -cos(s)

aspect_ratios = script_param(:aspect_ratios, [0.25, 1.0, 4.0])
gaps          = script_param(:gaps, [1.0, 0.1, 0.01])
orient        = script_param(:orient, :antiparallel)   # or :parallel
show_plots    = script_param(:plot, true)              # display the convergence plots at the end

N_ref = script_param(:N_ref, 1024)
M_ref = script_param(:M_ref, 32)
tol   = script_param(:tol, 1e-12)

Mspan = script_param(:Mspan, 5 .* (1:5))
Nspan = script_param(:Nspan, [])

configs = vec([(r, d) for d in gaps, r in aspect_ratios])

# -- Helpers -------------------------------------------------------------------

"Centers and orientations of a side-by-side pair with edge gap `gap` and body y semi-axis `b`."
function placement(b, gap)
    Δx  = 2b + gap
    xcm = [-Δx/2 + 0im, Δx/2 + 0im]
    θcm = orient == :parallel ? [π/2, π/2] : [-π/2, π/2]
    return xcm, θcm
end

"Rigid-body velocities `[ux uy ω]` of both particles, one row per particle."
function velocities(Ω, Γ_D, Γ_N, xcm, θcm)
    U, = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol,
                         self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))
    return permutedims(hcat([collect(u) for u in U]...))
end

# -- Main loop -----------------------------------------------------------------

err_M = zeros(length(Mspan), length(configs))
err_N = zeros(length(Nspan), length(configs))

for (j, (aspect_ratio, gap)) in enumerate(configs)

    a = sqrt(aspect_ratio)
    b = 1.0 / a
    xcm, θcm = placement(b, gap)

    Γ_D_ref, Γ_N_ref = surfer_boundary(a, b, s₀, N_ref)
    Ω_ref = ellipse(a, b, M_ref)
    U_ref = velocities(Ω_ref, Γ_D_ref, Γ_N_ref, xcm, θcm)

    @printf("\na/b = %g, gap = %.1e (measured %.3e), %s\n", aspect_ratio, gap,
            min_gap(xcm, θcm, a, b), orient)
    for i in 1:2
        @printf("  reference U%d = (%.12f, %.12f, %.12f)\n", i, U_ref[i, :]...)
    end
    flush(stdout)

    # (i) Expansion order M, boundary fixed at N_ref
    @printf("  expansion order (Ns = %d)\n%8s%14s\n", N_ref, "M", "error")
    for (i, M) in enumerate(Mspan)
        U = velocities(ellipse(a, b, M), Γ_D_ref, Γ_N_ref, xcm, θcm)
        err_M[i, j] = maximum(abs, U .- U_ref)
        @printf("%8d%14.2e\n", M, err_M[i, j])
        flush(stdout)
    end

    # (ii) Boundary resolution N, expansion fixed at M_ref
    @printf("  boundary resolution (M = %d)\n%8s%14s\n", M_ref, "Ns", "error")
    for (i, Ns) in enumerate(Nspan)
        Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)
        U = velocities(Ω_ref, Γ_D, Γ_N, xcm, θcm)
        err_N[i, j] = maximum(abs, U .- U_ref)
        @printf("%8d%14.2e\n", Ns, err_N[i, j])
        flush(stdout)
    end

end

# -- Write results -------------------------------------------------------------

outdir = script_param(:outdir, joinpath(@__DIR__, "output", "pair_gap"))
mkpath(outdir)

function write_errors(filename, label, span, err)
    outfile = joinpath(outdir, filename)
    open(outfile, "w") do io
        println(io, "# pair velocity error vs reference solution, side-by-side pair")
        @printf(io, "# N_ref = %d, M_ref = %d, tol = %.1e, orientation = %s\n", N_ref, M_ref, tol, orient)
        println(io, "# aspect ratios a/b: ", join(first.(configs), " "))
        println(io, "# edge gaps:         ", join(last.(configs), " "))
        println(io, "# columns: $label  err(config 1)  err(config 2) ...")
        for (i, n) in enumerate(span)
            @printf(io, "%d", n)
            for j in eachindex(configs)
                @printf(io, " %.15e", err[i, j])
            end
            println(io)
        end
    end
    println("Wrote ", relpath(outfile, @__DIR__))
end

write_errors("pair_gap_expansion.dat", "M", Mspan, err_M)
write_errors("pair_gap_boundary.dat", "N", Nspan, err_N)

# -- Plot results --------------------------------------------------------------

show_plots || exit()

panels = map(aspect_ratios) do r
    cols = findall(c -> first(c) == r, configs)
    labels = reshape(["gap = $(last(configs[j]))" for j in cols], 1, :)
    pM = plot(collect(Mspan), err_M[:, cols]; yscale=:log10, marker=:circle, label=labels,
              xlabel="expansion order, M", ylabel="error, |U - U_ref|", title="a/b = $r")
    pN = plot(collect(Nspan), err_N[:, cols]; yscale=:log10, marker=:circle, label=labels,
              xlabel="boundary resolution, N", ylabel="error, |U - U_ref|", title="a/b = $r")
    (pM, pN)
end

p = plot(first.(panels)..., last.(panels)...; layout=(2, length(aspect_ratios)),
         size=(400 * length(aspect_ratios), 700))
display(p)
readline()
