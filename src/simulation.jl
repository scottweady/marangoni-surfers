
"""
    suspension_simulate(xcm0, θcm0, Ω, Γ_D, Γ_N, g;
                        t_end=10.0, dt_min=1e-4, dt_max=10.0, safety=0.5,
                        plot_every=10.0, plot_mode=:concentration,
                        output_file="", plot_dir="", tol=1e-6)

Time-integrate a suspension of Marangoni surfers with adaptive-step RK2.

At each step, [`marangoni_solve`](@ref) gives the rigid-body velocity of every particle;
positions and orientations are advanced with a two-stage Runge-Kutta scheme. The next
time step is chosen from the minimum inter-particle gap and the maximum particle speed,
scaled by `safety`, and clamped to `[dt_min, dt_max]`. Frames are written or displayed
every `plot_every` time units, and per-frame particle configurations are appended to
`frame_<label>.dat` files under `plot_dir`.

# Arguments
- `xcm0`: initial complex particle centers.
- `θcm0`: initial particle orientations, in radians.
- `Ω`: Stokes domain boundary (shared by all particles).
- `Γ_D::Boundary`: reference Dirichlet arc boundary.
- `Γ_N::Boundary`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.
- `t_end`: total simulation time (default `10.0`).
- `dt_min`, `dt_max`: bounds on the adaptive time step (default `1e-4`, `10.0`).
- `safety`: safety factor applied to the gap/speed-based time-step estimate (default `0.5`).
- `plot_every`: time interval between saved/displayed frames (default `10.0`).
- `plot_mode`: `:concentration` plots the surfactant field; `:distance` plots the
  two-particle separation instead (default `:concentration`).
- `output_file`: if non-empty, append `(t, separation)` pairs here when `plot_mode ==
  :distance` (default `""`, no file output).
- `plot_dir`: if non-empty, save frames as PNGs and configuration `.dat` files here
  instead of displaying them (default `""`).
- `tol`: FMM/GMRES accuracy tolerance passed to [`marangoni_solve`](@ref) (default `1e-6`).

# Returns
- `(traj_xcm, traj_θcm, traj_t)`: the recorded trajectory of particle centers,
  orientations, and times.
"""
function suspension_simulate(xcm0, θcm0, Ω, Γ_D::Boundary, Γ_N::Boundary, g;
                  t_end=10.0,
                  dt_min=1e-4, dt_max=10.0,
                  safety=0.5,
                  plot_every=10.0,
                  plot_mode=:concentration,
                  output_file="",
                  plot_dir="",
                  tol=1e-6)

    xcm, θcm = copy(xcm0), copy(θcm0)
    Nparticles = length(xcm)

    t = 0.0
    next_plot = plot_every

    traj_xcm = [copy(xcm)]
    traj_θcm = [copy(θcm)]
    traj_t = [t]

    dt = dt_min
    nt = 0

    function save_config(label, t, xcm, θcm, ucm)
        dir  = isempty(plot_dir) ? "." : plot_dir
        path = joinpath(dir, "frame_$(label).dat")
        open(path, "w") do io
            for i in eachindex(xcm)
                @printf(io, "%.6f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                        t, real(xcm[i]), imag(xcm[i]), θcm[i], ucm[i][1], ucm[i][2], ucm[i][3])
            end
        end
    end

    # Plot frame 0
    if !isempty(plot_dir) || plot_mode != :distance
        U0, σ0, Γ0, ϕinf0 = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol)
        p = plot_concentration(xcm, σ0, ϕinf0, Γ0; title="t = 0.00")
        if !isempty(plot_dir)
            mkpath(plot_dir)
            savefig(p, joinpath(plot_dir, "frame_000000.png"))
        else
            display(p)
        end
        save_config("000000", 0.0, xcm, θcm, U0)
    end

    while t < t_end

        dt = min(dt, t_end - t, next_plot - t)

        # Runge-Kutta 2 for particle positions and orientations
        U1, σ, Γ, ϕinf = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol)
        Δc1 = dt .* [U1[i][1] + im*U1[i][2] for i in 1:Nparticles]
        Δψ1 = dt .* [U1[i][3]                for i in 1:Nparticles]

        U2, _, _, _ = marangoni_solve(xcm .+ Δc1, θcm .+ Δψ1, Ω, Γ_D, Γ_N, g; tol=tol)
        Δc2 = dt .* [U2[i][1] + im*U2[i][2] for i in 1:Nparticles]
        Δψ2 = dt .* [U2[i][3]                for i in 1:Nparticles]

        # Adaptive time step based on minimum gap and maximum speed
        function edge_gap(i, j)
            d = xcm[j] - xcm[i]
            d_hat = d / abs(d)
            p_i = vcat(Γ[2i-1].γ, Γ[2i].γ)
            p_j = vcat(Γ[2j-1].γ, Γ[2j].γ)
            x_i = p_i[argmax(real.(conj(d_hat) .* (p_i .- xcm[i])))]
            x_j = p_j[argmax(real.(conj(d_hat) .* (xcm[j] .- p_j)))]
            return abs(x_i - x_j)
        end

        d_min     = Nparticles > 1 ? minimum(edge_gap(i, j) for i in 1:Nparticles for j in i+1:Nparticles) : Inf
        max_speed = maximum(hypot(U1[i][1], U1[i][2]) for i in 1:Nparticles)
        dt_next   = clamp(safety * d_min / max(max_speed, 1e-14), dt_min, dt_max)

        xcm_plot = copy(xcm)
        xcm .+= 0.5 .* (Δc1 .+ Δc2)
        θcm .+= 0.5 .* (Δψ1  .+ Δψ2)
        t += dt
        dt = dt_next
        nt += 1
        push!(traj_xcm, copy(xcm))
        push!(traj_θcm, copy(θcm))
        push!(traj_t,   t)

        @printf("  t = %.4f  dt = %.2e\n", t, dt_next)

        if t >= next_plot
            if plot_mode == :distance
                dists = [abs(traj_xcm[k][2] - traj_xcm[k][1]) for k in eachindex(traj_t)]
                if !isempty(output_file)
                    open(output_file, "a") do io
                        @printf(io, "%.6f %.6f\n", traj_t[end], dists[end])
                    end
                end
                p = plot(traj_t, dists .- 2.0, xlabel="t", ylabel="|r₂ − r₁|",
                         legend=false, title="particle separation")
            else
                p = plot_concentration(xcm_plot, σ, ϕinf, Γ;
                                       title="t = $(round(t, digits=2))")
            end
            frame = lpad(length(traj_t), 6, '0')
            if !isempty(plot_dir)
                mkpath(plot_dir)
                savefig(p, joinpath(plot_dir, "frame_$(frame).png"))
            else
                display(p)
            end
            save_config(frame, t, xcm, θcm, U1)
            next_plot += plot_every
        end
    end

    return traj_xcm, traj_θcm, traj_t

end
