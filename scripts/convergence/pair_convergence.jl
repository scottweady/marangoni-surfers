# Self-convergence of the rigid-body velocities of an interacting pair of elliptical
# surfers, at the closest pair configurations of the trajectory runs. Each configuration is
# an (aspect ratio, center separation) pair at fixed area, with the pair oriented parallel
# or antiparallel (`orient`). The error is the largest difference over both particles and
# all components (ux, uy, ω) from a reference solve at (`M_ref`, `N_ref`), measured
# separately in the Stokes expansion order (`Mspan`, with the boundary at `N_ref`) and in
# the boundary resolution (`Nspan`, with the expansion at `M_ref`). The FMM/GMRES
# tolerance `tol` is set well below the errors being measured.
#
# Output: output/pair_convergence_<expansion|boundary>.dat
#   Header: N_ref, M_ref, tol, orientation, aspect ratios, separations, edge gaps
#   Rows:   M (or N)  followed by the error for each configuration

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀   = π/2
g(s) = -cos(s)

configs    = script_param(:configs, [(4.0, 1.1)])
orient     = script_param(:orient, :antiparallel)   # or :parallel
show_plots = script_param(:plot, true)              # display the convergence plots at the end

N_ref = script_param(:N_ref, 128)
M_ref = script_param(:M_ref, 32)
tol   = script_param(:tol, 1e-13)

Mspan = script_param(:Mspan, 4:4:24)
Nspan = script_param(:Nspan, 20:20:100)

# -- Helpers -------------------------------------------------------------------

"Initial centers and orientations of the pair, matching scripts/pair/pair.jl."
function placement(Δx)
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

"Minimum distance between the boundaries of the two placed particles."
function edge_gap(a, b, xcm, θcm)
    z = [a*cos(s) + im*b*sin(s) for s in range(0, 2π; length=4001)]
    z1 = exp(im*θcm[1]) .* z .+ xcm[1]
    z2 = exp(im*θcm[2]) .* z .+ xcm[2]
    return sqrt(minimum(abs2.(z1 .- permutedims(z2))))
end

# -- Main loop -----------------------------------------------------------------

err_M = zeros(length(Mspan), length(configs))
err_N = zeros(length(Nspan), length(configs))
gaps  = zeros(length(configs))

for (j, (aspect_ratio, Δx)) in enumerate(configs)

    a = sqrt(aspect_ratio)
    b = 1.0 / a
    xcm, θcm = placement(Δx)
    gaps[j] = edge_gap(a, b, xcm, θcm)

    Γ_D_ref, Γ_N_ref = surfer_boundary(a, b, s₀, N_ref)
    Ω_ref = ellipse(a, b, M_ref)
    U_ref = velocities(Ω_ref, Γ_D_ref, Γ_N_ref, xcm, θcm)

    @printf("\na/b = %g, Δx = %g (gap %.3f), %s\n", aspect_ratio, Δx, gaps[j], orient)
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

outdir = script_param(:outdir, joinpath(@__DIR__, "output"))
mkpath(outdir)

function write_errors(filename, label, span, err)
    outfile = joinpath(outdir, filename)
    open(outfile, "w") do io
        println(io, "# pair velocity error vs reference solution")
        @printf(io, "# N_ref = %d, M_ref = %d, tol = %.1e, orientation = %s\n", N_ref, M_ref, tol, orient)
        println(io, "# aspect ratios a/b: ", join(first.(configs), " "))
        println(io, "# separations Δx:    ", join(last.(configs), " "))
        println(io, "# edge gaps:         ", join(round.(gaps; digits=6), " "))
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

write_errors("pair_convergence_expansion.dat", "M", Mspan, err_M)
write_errors("pair_convergence_boundary.dat", "N", Nspan, err_N)

# -- Plot results --------------------------------------------------------------

show_plots || exit()

labels = reshape(["a/b = $(r), Δx = $(d)" for (r, d) in configs], 1, :)

p1 = plot(collect(Mspan), err_M;
          xlabel = "expansion order, M",
          ylabel = "error, |U - U_ref|",
          yscale = :log10,
          marker = :circle,
          label  = labels,
          legend = :topright)

p2 = plot(collect(Nspan), err_N;
          xlabel = "boundary resolution, N",
          ylabel = "error, |U - U_ref|",
          yscale = :log10,
          marker = :circle,
          label  = labels,
          legend = :topright)

p = plot(p1, p2; layout=(1, 2), size=(900, 380))
display(p)
readline()
