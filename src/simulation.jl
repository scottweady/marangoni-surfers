
"""
    suspension_simulate(xcm0, θcm0, Ω, Γ_D, Γ_N, g;
                        t_end=10.0, dt_min=1e-4, dt_max=10.0, safety=0.5,
                        plot_every=10.0, plot_mode=:concentration,
                        output_file="", plot_dir="", tol=1e-3,
                        colors=fill(:white, length(xcm0)), cmap=:viridis, pad=4.0)

Time-integrate a suspension of Marangoni surfers with adaptive-step RK2.

At each step, [`marangoni_solve`](@ref) gives the rigid-body velocity of every particle;
positions and orientations are advanced with a two-stage Runge-Kutta scheme. The next
time step is chosen so that no pair of particles closes its gap within the step: for each
pair, the relative velocity is projected onto the line joining the particles, and (when
approaching) the step is limited by gap / closing-speed, scaled by `safety` and clamped to
`[dt_min, dt_max]`. Frames are displayed every
`plot_every` time units, and per-frame particle configurations are appended to
`frame_<label>.dat` files under `plot_dir`, if given.

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
- `plot_dir`: if non-empty, save configuration `.dat` files here (default `""`, no file
  output). Frames are displayed regardless.
- `tol`: FMM/GMRES accuracy tolerance passed to [`marangoni_solve`](@ref) (default `1e-3`).
- `colors`: per-particle fill colors (default `:white` for all particles).
- `cmap`: heatmap colormap (default `:viridis`).
- `pad`: margin added around the current particle centers when computing each frame's
  plot extent, so all particles stay in view as they move (default `4.0`).

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
                  tol=1e-3,
                  colors=fill(:white, length(xcm0)),
                  cmap=:viridis,
                  pad=4.0)

    xcm, θcm = copy(xcm0), copy(θcm0)
    Nparticles = length(xcm)

    function plot_lims(xcm_c)
        xl = (minimum(real.(xcm_c)) - pad, maximum(real.(xcm_c)) + pad)
        yl = (minimum(imag.(xcm_c)) - pad, maximum(imag.(xcm_c)) + pad)
        return xl, yl
    end

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
        xlims0, ylims0 = plot_lims(xcm)
        p0 = plot_concentration(xcm, σ0, ϕinf0, Γ0; title="t = 0.00", colors=colors, cmap=cmap, xlims=xlims0, ylims=ylims0)
        display(p0)
        if !isempty(plot_dir)
            mkpath(plot_dir)
            savefig(p0, joinpath(plot_dir, "frame_000000.png"))
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

        # Adaptive time step: shrink dt so no pair of particles closes its gap this step,
        # using the component of relative velocity along the line joining each pair.
        function edge_gap(i, j, d_hat)
            p_i = vcat(Γ[2i-1].γ, Γ[2i].γ)
            p_j = vcat(Γ[2j-1].γ, Γ[2j].γ)
            x_i = p_i[argmax(real.(conj(d_hat) .* (p_i .- xcm[i])))]
            x_j = p_j[argmax(real.(conj(d_hat) .* (xcm[j] .- p_j)))]
            return abs(x_i - x_j)
        end

        dt_next = dt_max
        for i in 1:Nparticles, j in i+1:Nparticles
            d_hat = (xcm[j] - xcm[i]) / abs(xcm[j] - xcm[i])
            gap   = edge_gap(i, j, d_hat)
            v_i   = U1[i][1] + im*U1[i][2]
            v_j   = U1[j][1] + im*U1[j][2]
            closing_speed = real(conj(d_hat) * (v_i - v_j))  # > 0 means approaching
            if closing_speed > 0
                dt_next = min(dt_next, safety * gap / closing_speed)
            end
        end
        dt_next = clamp(dt_next, dt_min, dt_max)

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

        frame = lpad(length(traj_t), 6, '0')
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
                xlims, ylims = plot_lims(xcm_plot)
                p = plot_concentration(xcm_plot, σ, ϕinf, Γ;
                                       title="t = $(round(t, digits=2))", colors=colors, cmap=cmap, xlims=xlims, ylims=ylims)
            end
            display(p)
            if !isempty(plot_dir)
                mkpath(plot_dir)
                savefig(p, joinpath(plot_dir, "frame_$(frame).png"))
            end
            save_config(frame, t, xcm, θcm, U1)
            next_plot += plot_every
        end
    end

    return traj_xcm, traj_θcm, traj_t

end
