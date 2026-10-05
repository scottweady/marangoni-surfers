# Trajectories of two interacting elliptical surfers, integrated with midpoint RK2.
#
# Sweeps aspect ratio over `aspect_ratios` at fixed area and initial separation over
# `Δx_range`, with the pair parallel or antiparallel (`orient`).
#
# Output: output/a<a>_b<b>_d<Δx>.dat, one row per time step
#   t  x1  y1  θ1  x2  y2  θ2  gap
# where gap is the minimum boundary-to-boundary distance.

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀   = π/2
g(s) = -cos(s)

aspect_ratios = [10.0]
Δx_range = [1.65]

# orient = :parallel
orient = :antiparallel

# -- Discretization ------------------------------------------------------------

Ns = 64      # arc quadrature points
M  = 16      # spectral degree for ellipse

# Time integration parameters
t_max = 500  # max simulation time
dt_max = 1.0   # max time step
dt_min = 1e-3    # min time step
nplt = 10        # plot/write every nplt steps
cfl_scale = 0.5  # safety factor in U*dt < d using boundary-to-boundary gap d
stop_dist = 1e-2  # stop when nearest boundary-to-boundary distance is below this
max_separation = Inf  # stop when center-to-center distance exceeds this

# -- Main loop: aspect ratio and separation ------------------------------------

for aspect_ratio in aspect_ratios

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    @printf("(a, b) = (%.4f, %.4f)\n", a, b)

    # Define ellipse and arc boundaries (body frame, θ=0)
    γ(s) = a*cos(s) + im*b*sin(s)
    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

    # Boundary points in body frame for nearest-distance checks
    z_boundary_body = vcat(Γ_D.γ, Γ_N.γ)

    # Densely sampled body-frame curves for plotting: full outline and Dirichlet arc
    z_outline_body = γ.(range(0, 2π; length=201))
    z_arc_body     = γ.(range(s₀, 2π - s₀; length=101))

    for Δx in Δx_range

        @printf("  Δx = %.1f\n", Δx)

        # -- Initial conditions ------------------------------------------------

        xcm  = [-Δx/2 + 0im, Δx/2 + 0im]   # center of mass positions

        if orient == :parallel
            θcm  = [π/2, π/2]
        elseif orient == :antiparallel
            θcm  = [-π/2, π/2]
        end

        U_all = Vector{Vector{Float64}}(undef, 2)

        # Storage for trajectory
        times = Float64[]
        trajectory = Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64}[]

        # Record initial condition (t=0)
        z_boundary_1_0 = exp(im * θcm[1]) .* z_boundary_body .+ xcm[1]
        z_boundary_2_0 = exp(im * θcm[2]) .* z_boundary_body .+ xcm[2]
        dist_0 = sqrt(minimum(abs2.(z_boundary_1_0 .- permutedims(z_boundary_2_0))))
        push!(times, 0.0)
        push!(trajectory, (
            real(xcm[1]), imag(xcm[1]), θcm[1],
            real(xcm[2]), imag(xcm[2]), θcm[2],
            dist_0
        ))

        # Setup file for writing
        outdir = joinpath(@__DIR__, "output")
        mkpath(outdir)
        outfile = joinpath(outdir, @sprintf("a%.4f_b%.4f_d%.1f.dat", a, b, Δx))
        io = open(outfile, "w")

        # -- Time integration --------------------------------------------------

        t = 0.0
        step = 0
        last_written = 0
        max_steps = Int(ceil(t_max / dt_max)) + 1000

        while t < t_max && step < max_steps

            # Solve for velocities at current configuration
            U, σ, Γ, ϕinf = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)

            # Extract velocities: U[i] = [vx, vy, ω]
            for i in 1:2
                U_all[i] = U[i]
            end

            @printf(" xcm1 = (%.15e, %.15e)  θ1 = %.15e\n", real(xcm[1]), imag(xcm[1]), θcm[1])
            @printf(" U1 = (%.15e, %.15e, %.15e) \n", U[1][1], U[1][2], U[1][3])

            # Current boundary-to-boundary gap d (nearest points on both particles)
            z_boundary_1 = exp(im * θcm[1]) .* z_boundary_body .+ xcm[1]
            z_boundary_2 = exp(im * θcm[2]) .* z_boundary_body .+ xcm[2]
            dist = sqrt(minimum(abs2.(z_boundary_1 .- permutedims(z_boundary_2))))

            if dist < stop_dist
                @printf("    Boundary-gap threshold reached (gap=%.6e < %.1e) at t = %.6f\n", dist, stop_dist, t)
                break
            end

            if abs(xcm[1] - xcm[2]) > max_separation
                @printf("    Max separation threshold reached (sep=%.6e > %.1f) at t = %.6f\n", abs(xcm[1] - xcm[2]), max_separation, t)
                break
            end

            # Adaptive time step: enforce U*dt < d with d = boundary-to-boundary gap
            v_max = maximum(sqrt(sum(U_all[i][1:2].^2)) for i in 1:2)
            dt = min(dt_max, cfl_scale * dist / v_max)

            remaining = t_max - t
            dt = min(dt, remaining)

            if dt < dt_min && remaining > dt_min
                dt = dt_min
            end

            # Midpoint RK2 step
            # Predict
            vx1, vy1 = U_all[1][1], U_all[1][2]
            vx2, vy2 = U_all[2][1], U_all[2][2]
            xcm_mid = [
                xcm[1] + 0.5 * dt * (vx1 + im * vy1),
                xcm[2] + 0.5 * dt * (vx2 + im * vy2),
            ]
            θcm_mid  = θcm .+ 0.5 * dt .* [U_all[1][3], U_all[2][3]]

            # Eval at midpoint
            U_mid, _, _, _ = marangoni_solve(xcm_mid, θcm_mid, Ω, Γ_D, Γ_N, g)

            # Update
            for i in 1:2
                xcm[i] += dt * (U_mid[i][1] + im * U_mid[i][2])
                θcm[i] += dt * U_mid[i][3]
            end

            t += dt
            step += 1

            z_boundary_1 = exp(im * θcm[1]) .* z_boundary_body .+ xcm[1]
            z_boundary_2 = exp(im * θcm[2]) .* z_boundary_body .+ xcm[2]
            dist = sqrt(minimum(abs2.(z_boundary_1 .- permutedims(z_boundary_2))))

            # Save trajectory (every step to array)
            push!(times, t)
            push!(trajectory, (
                real(xcm[1]), imag(xcm[1]), θcm[1],
                real(xcm[2]), imag(xcm[2]), θcm[2],
                dist
            ))

            plt = plot(aspect_ratio=:equal, legend=false, xlabel="x", ylabel="y",
                       title=@sprintf("t = %.2f", t))
            xlims!(-20, 20)
            ylims!(-20, 20)

            for i in 1:2
                z_outline = exp(im * θcm[i]) .* z_outline_body .+ xcm[i]
                z_arc     = exp(im * θcm[i]) .* z_arc_body     .+ xcm[i]
                plot!(plt, Shape(real.(z_outline), imag.(z_outline)),
                      fillcolor=:lightgray, linecolor=:black, linewidth=1.5)
                plot!(plt, real.(z_arc), imag.(z_arc), linecolor=:crimson, linewidth=3)
            end
            display(plt)

            @printf("    step=%d  t=%.6f  dt=%.3e  gap=%.6e\n", step, t, dt, dist)

            # Write to file every 10 steps
            if mod(step, nplt) == 0
                # Write lines since the last write
                for i in last_written+1:length(times)
                    t_val = times[i]
                    x1, y1, θ1, x2, y2, θ2, dist_val = trajectory[i]
                    @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                            t_val, x1, y1, θ1, x2, y2, θ2, dist_val)
                end
                last_written = length(times)
                flush(io)
            end

            if dist < stop_dist
                @printf("    Boundary-gap threshold reached (gap=%.6e < %.1e) at t = %.6f\n", dist, stop_dist, t)
                break
            end

            if abs(xcm[1] - xcm[2]) > max_separation
                @printf("    Max separation threshold reached (sep=%.6e > %.1f) at t = %.6f\n", abs(xcm[1] - xcm[2]), max_separation, t)
                break
            end

        end

        # Write any remaining steps to file
        for i in last_written+1:length(times)
            t_val = times[i]
            x1, y1, θ1, x2, y2, θ2, dist_val = trajectory[i]
            @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                    t_val, x1, y1, θ1, x2, y2, θ2, dist_val)
        end

        close(io)

        println("    Wrote $outfile ($(length(times)) time steps)")
    end
end

println("\nDone.")
readline()
