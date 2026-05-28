"""
Generate trajectory data for Figure 3: two interacting ellipses.

Varies aspect ratio a ∈ [1, 4] with fixed area πab = π (so b = 1/a).
Varies initial separation d ∈ [5, 10].
Outputs trajectories to .dat files: fig3_a{a:.1f}_b{b:.1f}_d{d:.1f}.dat

Each file contains columns:
  t  x1_cm  y1_cm  θ1  x2_cm  y2_cm  θ2  separation
"""

using ProjectedSphericalHarmonics
using Printf

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

# ── Parameters ───────────────────────────────────────────────────────────────

α_arc = π/2   # arc angle (Dirichlet domain)
g(θ) = -cos(θ)  # surfactant distribution

aspect_ratios = [1.0]
Δx_range = [10.0]

# ── Discretization ─────────────────────────────────────────────────────────
Ns = 256      # arc quadrature points
M  = 32       # spectral degree for ellipse

# Time integration parameters
t_max = 10000.0  # max simulation time
dt_max = 100.0   # max time step
dt_min = 1e-3    # min time step
cfl_scale = 0.5  # safety factor in U*dt < d using boundary-to-boundary gap d
stop_dist = 1e-2  # stop when nearest boundary-to-boundary distance is below this

# ── Main loop: aspect ratio and separation ───────────────────────────────────

for a in aspect_ratios

    b = 1.0 / a  # fixed area: πab = π
    
    @printf("(a, b) = (%.4f, %.4f)\n", a, b)
    
    # Define ellipse and arc boundaries (body frame, θ=0)
    γ(s)   =  a*cos(s) + im*b*sin(s)
    γs(s)  = -a*sin(s) + im*b*cos(s)
    γss(s) = -a*cos(s) - im*b*sin(s)
    
    Ω = ellipse(a, b, M)
    Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
    Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

    # Boundary points in body frame for nearest-distance checks
    z_boundary_body = vcat(Γ_D.γ, Γ_N.γ)
    
    for Δx in Δx_range
      
        @printf("  Δx = %.1f\n", Δx)
        
        # ── Initial conditions ───────────────────────────────────────────────
        
        xcm  = [-Δx/2 + 0im, Δx/2 + 0im]   # center of mass positions
        θcm  = [π/2, π/2]                 # both oriented vertically
        U_all = Vector{Vector{Float64}}(undef, 2)
        
        # Storage for trajectory
        times = Float64[]
        trajectory = Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64}[]
        
        # Setup file for writing
        out_dir = joinpath(@__DIR__, "fig3")
        mkpath(out_dir)
        filename = @sprintf("fig3_a%.1f_b%.4f_d%.1f.dat", a, b, Δx)
        filepath = joinpath(out_dir, filename)
        io = open(filepath, "w")
        
        # ── Time integration ─────────────────────────────────────────────────
        
        t = 0.0
        step = 0
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

            @printf("    step=%d  t=%.6f  dt=%.3e  gap=%.6e\n", step, t, dt, dist)
            
            # Write to file every 10 steps
            if mod(step, 10) == 0
                # Write last 10 lines to file
                idx_start = max(1, length(times) - 9)
                for i in idx_start:length(times)
                    t_val = times[i]
                    x1, y1, θ1, x2, y2, θ2, dist_val = trajectory[i]
                    @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                            t_val, x1, y1, θ1, x2, y2, θ2, dist_val)
                end
                flush(io)
            end

            if dist < stop_dist
                @printf("    Boundary-gap threshold reached (gap=%.6e < %.1e) at t = %.6f\n", dist, stop_dist, t)
                break
            end

        end
        
        # Write any remaining steps to file
        idx_start = div(step, 10) * 10 + 1
        for i in idx_start:length(times)
            t_val = times[i]
            x1, y1, θ1, x2, y2, θ2, dist_val = trajectory[i]
            @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                    t_val, x1, y1, θ1, x2, y2, θ2, dist_val)
        end
        
        close(io)
        
        println("    Wrote $filename ($(length(times)) time steps)")
    end
end

println("\nDone.")
