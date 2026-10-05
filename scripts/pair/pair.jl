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

aspect_ratios = [10.0]
Δx_range      = [1.65]
orient        = :antiparallel   # or :parallel

# -- Discretization ------------------------------------------------------------

Ns = 128   # arc quadrature points
M  = 32    # Stokes expansion order

t_max          = 500
dt_min, dt_max = 1e-3, 1.0
cfl_scale      = 0.5    # dt ≤ cfl_scale * gap / max speed
stop_dist      = 1e-2   # stop when the boundary gap falls below this
max_separation = Inf    # stop when the center-to-center distance exceeds this

outdir = joinpath(@__DIR__, "output", string(orient))
mkpath(outdir)

# -- Helpers -------------------------------------------------------------------

"Place body-frame points `z` at center `xc` with orientation `θ`."
place(z, xc, θ) = exp(im * θ) .* z .+ xc

"Minimum distance between the boundary points `z_body` of the two placed particles."
function gap(z_body, xcm, θcm)
    z1 = place(z_body, xcm[1], θcm[1])
    z2 = place(z_body, xcm[2], θcm[2])
    return sqrt(minimum(abs2.(z1 .- permutedims(z2))))
end

"Complex translational and angular velocities from the output of `marangoni_solve`."
velocities(U) = ([u[1] + im * u[2] for u in U], [u[3] for u in U])

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

# -- Main loop: aspect ratio and separation ------------------------------------

for aspect_ratio in aspect_ratios

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    γ(s) = a*cos(s) + im*b*sin(s)
    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

    z_body    = vcat(Γ_D.γ, Γ_N.γ)
    z_outline = γ.(range(0, 2π; length=201))
    z_arc     = γ.(range(s₀, 2π - s₀; length=101))

    for Δx in Δx_range

        @printf("a/b = %g, Δx = %g\n", aspect_ratio, Δx)

        xcm = [-Δx/2 + 0im, Δx/2 + 0im]
        θcm = orient == :parallel ? [π/2, π/2] : [-π/2, π/2]

        outfile = joinpath(outdir, @sprintf("a%.4f_b%.4f_d%.1f.dat", a, b, Δx))
        io = open(outfile, "w")

        t, step = 0.0, 0
        d = gap(z_body, xcm, θcm)
        write_row(io, t, xcm, θcm, d)

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

            U, = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
            v, ω = velocities(U)
            dt = min(clamp(cfl_scale * d / maximum(abs, v), dt_min, dt_max), t_max - t)

            U_mid, = marangoni_solve(xcm .+ dt/2 .* v, θcm .+ dt/2 .* ω, Ω, Γ_D, Γ_N, g)
            v_mid, ω_mid = velocities(U_mid)
            xcm .+= dt .* v_mid
            θcm .+= dt .* ω_mid

            t += dt
            step += 1
            d = gap(z_body, xcm, θcm)
            write_row(io, t, xcm, θcm, d)

            display(plot_pair(xcm, θcm, t, z_outline, z_arc))
            @printf("  t = %10.4f   dt = %.2e   gap = %.2e\n", t, dt, d)

        end

        close(io)
        println("Wrote ", relpath(outfile, @__DIR__), " (", step + 1, " steps)")

    end
end

readline()
