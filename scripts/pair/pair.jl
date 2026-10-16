# Trajectories of two interacting elliptical surfers, integrated with midpoint RK2.
#
# Sweeps aspect ratio over `aspect_ratios` at fixed area and initial separation over
# `Δx_range`, with the pair parallel or antiparallel (`orient`).
#
# Output: output/<orient>/a<a>_b<b>_d<Δx>.dat, one row per time step
#   t  x1  y1  θ1  x2  y2  θ2  gap
# where gap is the minimum boundary-to-boundary distance.

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀   = π/2
g(s) = -cos(s)

aspect_ratios = script_param(:aspect_ratios, [4.0])
Δx_range      = script_param(:Δx_range, [4.0])
orient        = script_param(:orient, :parallel)   # or :parallel
show_plots    = script_param(:plot, true)              # display the pair at each step

# -- Discretization ------------------------------------------------------------

Ns = 64    # arc quadrature points
M  = 16    # Stokes expansion order

t_max          = script_param(:t_max, 500.0)
dt_min, dt_max = 1e-3, 5.0
safety         = 0.5    # fraction of the gap that may close in one step
stop_dist      = 1e-3   # stop when the boundary gap falls below this
max_separation = Inf    # stop when the center-to-center distance exceeds this

outdir = script_param(:outdir, joinpath(@__DIR__, "output", string(orient)))
mkpath(outdir)

# -- Helpers -------------------------------------------------------------------

"Place body-frame points `z` at center `xc` with orientation `θ`."
place(z, xc, θ) = exp(im * θ) .* z .+ xc

"Append one trajectory row `t x1 y1 θ1 x2 y2 θ2 gap` to `io`."
function write_row(io, t, xcm, θcm, d)
    @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
            t, real(xcm[1]), imag(xcm[1]), θcm[1], real(xcm[2]), imag(xcm[2]), θcm[2], d)
    flush(io)
end

"Plot both particles with their Dirichlet arcs highlighted."
function plot_pair(xcm, θcm, t, z_outline, z_arc)
    plt = plot(aspect_ratio=:equal, legend=false, xlabel="x", ylabel="y",
               title=@sprintf("t = %.2f", t), xlims=(-20, 20), ylims=(-20, 20))
    for i in 1:2
        zo = place(z_outline, xcm[i], θcm[i])
        za = place(z_arc, xcm[i], θcm[i])
        plot!(plt, Shape(real.(zo), imag.(zo)), fillcolor=:lightgray, linecolor=:black, linewidth=1.5)
        plot!(plt, real.(za), imag.(za), linecolor=:crimson, linewidth=3)
    end
    return plt
end

"Plot the boundary gap history on a log scale."
function plot_gap(ts, ds)
    return plot(ts, ds, legend=false, xlabel="t", ylabel="gap",
                linecolor=:black, linewidth=1.5)
end

# -- Main loop: aspect ratio and separation ------------------------------------

for aspect_ratio in aspect_ratios

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    γ(s) = a*cos(s) + im*b*sin(s)
    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)
    self_op = marangoni_self_operator(Ω, Γ_D, Γ_N)
    velocity(x, θ) = rigid_velocities(first(marangoni_solve(x, θ, Ω, Γ_D, Γ_N, g; self_op=self_op)))

    z_outline = γ.(range(0, 2π; length=201))
    z_arc     = γ.(range(s₀, 2π - s₀; length=101))

    for Δx in Δx_range

        @printf("a/b = %g, Δx = %g\n", aspect_ratio, Δx)

        xcm = [-Δx/2 + 0im, Δx/2 + 0im]
        θcm = orient == :parallel ? [π/2, π/2] : [-π/2, π/2]

        outfile = joinpath(outdir, @sprintf("a%.4f_b%.4f_d%.1f.dat", a, b, Δx))
        io = open(outfile, "w")

        t, step = 0.0, 0
        d = min_gap(xcm, θcm, a, b)
        write_row(io, t, xcm, θcm, d)
        ts, ds = [t], [d]

        max_steps = ceil(Int, t_max / dt_max) + 1000
        while t < t_max && step < max_steps

            if d < stop_dist
                @printf("  stopped at t = %.4f: gap %.2e < %.0e\n", t, d, stop_dist)
                break
            end
            if abs(xcm[1] - xcm[2]) > max_separation
                @printf("  stopped at t = %.4f: separation %.2e > %.2e\n", t, abs(xcm[1] - xcm[2]), max_separation)
                break
            end

            v, ω = velocity(xcm, θcm)
            next = midpoint_step(velocity, xcm, θcm, v, ω, a, b; safety=safety, dt_min=dt_min,
                                 dt_max=dt_max, dt_cap=t_max - t)
            if next === nothing
                @printf("  stopped at t = %.4f: no step of at least dt_min keeps the gap\n", t)
                break
            end
            xcm, θcm, dt = next

            t += dt
            step += 1
            d = min_gap(xcm, θcm, a, b)
            write_row(io, t, xcm, θcm, d)
            push!(ts, t)
            push!(ds, d)

            show_plots && display(plot(plot_pair(xcm, θcm, t, z_outline, z_arc), plot_gap(ts, ds);
                                       layout=(1, 2), size=(1000, 450)))
            @printf("  t = %10.4f   dt = %.2e   gap = %.2e\n", t, dt, d)

        end

        close(io)
        println("Wrote ", relpath(outfile, @__DIR__), " (", step + 1, " steps)")

    end
end

show_plots && readline()
