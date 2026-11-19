
"""
    rigid_velocities(U)

Complex translational and angular velocities from the output of [`marangoni_solve`](@ref).

# Arguments
- `U`: per-particle rigid-body velocities `U[i] = (ux, uy, ω)`.

# Returns
- `(v, ω)`: complex translational velocities and angular velocities.
"""
rigid_velocities(U) = ([u[1] + im * u[2] for u in U], [u[3] for u in U])

"""
    boundary_gap(a, b, xc1, θ1, xc2, θ2; n=64)

Minimum distance between the boundaries of two elliptical particles with semi-axes `a`, `b`.

The closest pair of `n` uniformly spaced boundary points is refined by damped Newton on the
squared distance in the two ellipse parameters. Overlapping particles have zero gap, up to
rounding.

# Arguments
- `a`, `b`: semi-axes of both particles, in the body frame.
- `xc1`, `θ1`: complex center and orientation of the first particle.
- `xc2`, `θ2`: complex center and orientation of the second particle.
- `n`: number of boundary points used for the initial guess (default `64`).

# Returns
- The boundary-to-boundary distance.
"""
function boundary_gap(a, b, xc1, θ1, xc2, θ2; n=64)

    e1, e2 = cis(θ1), cis(θ2)
    ell(s)  =  a * cos(s) + im * b * sin(s)
    ells(s) = -a * sin(s) + im * b * cos(s)
    inside(z, xc, e) = (w = conj(e) * (z - xc); (real(w) / a)^2 + (imag(w) / b)^2 < 1)

    s = 2π .* (0:n-1) ./ n
    p = xc1 .+ e1 .* ell.(s)
    q = xc2 .+ e2 .* ell.(s)
    (any(z -> inside(z, xc2, e2), p) || any(z -> inside(z, xc1, e1), q)) && return 0.0

    D2 = abs2.(p .- transpose(q))
    k  = argmin(D2)
    s1, s2 = s[k[1]], s[k[2]]
    d2 = D2[k]

    # Damped Newton (Levenberg-Marquardt) on f(s1, s2) = |P(s1) - Q(s2)|²/2: the damping μ
    # is raised until the step decreases f, and relaxed after each accepted step
    rdot(z, w) = real(conj(z) * w)
    μ = 0.0
    for _ in 1:100
        d  = xc1 + e1 * ell(s1) - xc2 - e2 * ell(s2)
        P1, Q1 = e1 * ells(s1), e2 * ells(s2)
        P2, Q2 = -e1 * ell(s1), -e2 * ell(s2)
        g1, g2 = rdot(d, P1), -rdot(d, Q1)
        H11 = abs2(P1) + rdot(d, P2)
        H22 = abs2(Q1) - rdot(d, Q2)
        H12 = -rdot(P1, Q1)
        μ₀  = 1e-12 * (abs(H11) + abs(H22))
        Δ1, Δ2, accepted = 0.0, 0.0, false
        for _ in 1:60
            A11, A22 = H11 + μ, H22 + μ
            det = A11 * A22 - H12^2
            if A11 > 0 && det > 0
                Δ1 = -(A22 * g1 - H12 * g2) / det
                Δ2 = -(A11 * g2 - H12 * g1) / det
                d2_new = abs2(xc1 + e1 * ell(s1 + Δ1) - xc2 - e2 * ell(s2 + Δ2))
                if d2_new <= d2
                    s1, s2, d2 = s1 + Δ1, s2 + Δ2, d2_new
                    μ /= 4
                    accepted = true
                    break
                end
            end
            μ = max(4μ, μ₀)
        end
        (!accepted || abs(Δ1) + abs(Δ2) < 1e-14) && break
    end

    return sqrt(d2)

end

"""
    min_gap(xcm, θcm, a, b)

Smallest boundary-to-boundary distance over all pairs of particles.

# Arguments
- `xcm`: complex particle centers.
- `θcm`: particle orientations, in radians.
- `a`, `b`: semi-axes of the particles.

# Returns
- The minimum gap, or `Inf` for a single particle.
"""
function min_gap(xcm, θcm, a, b)
    R = max(a, b)
    gap = Inf
    for i in eachindex(xcm), j in i+1:length(xcm)
        abs(xcm[i] - xcm[j]) - 2R >= gap && continue
        gap = min(gap, boundary_gap(a, b, xcm[i], θcm[i], xcm[j], θcm[j]))
    end
    return gap
end

"""
    step_size(xcm, θcm, v, ω, a, b; safety=0.5, dt_max=1.0)

Time step over which, at frozen velocities, no pair of particles closes more than a
fraction `safety` of its gap,

    dt = safety * min_{i<j} gap_ij / (|v_i - v_j| + R (|ω_i| + |ω_j|)),

with `R = max(a, b)`, so that the denominator bounds the relative speed of any two boundary
points. The step is capped at `dt_max`.

# Arguments
- `xcm`: complex particle centers.
- `θcm`: particle orientations, in radians.
- `v`, `ω`: complex translational velocities and angular velocities.
- `a`, `b`: semi-axes of the particles.
- `safety`: fraction of a gap that may close in one step (default `0.5`).
- `dt_max`: largest allowed step (default `1.0`).

# Returns
- The time step.
"""
function step_size(xcm, θcm, v, ω, a, b; safety=0.5, dt_max=1.0)
    R = max(a, b)
    τ = dt_max / safety
    for i in eachindex(xcm), j in i+1:length(xcm)
        speed = abs(v[i] - v[j]) + R * (abs(ω[i]) + abs(ω[j]))
        speed == 0 && continue
        (abs(xcm[i] - xcm[j]) - 2R) / speed >= τ && continue
        τ = min(τ, boundary_gap(a, b, xcm[i], θcm[i], xcm[j], θcm[j]) / speed)
    end
    return safety * τ
end

"""
    _gaps_kept(xcm0, θcm0, xcm1, θcm1, a, b, safety)

Whether every pair's gap after a step is at least `1 - safety` times its gap before it.
Pairs are skipped when center-distance bounds on the gaps already guarantee this.
"""
function _gaps_kept(xcm0, θcm0, xcm1, θcm1, a, b, safety)
    R, r = max(a, b), min(a, b)
    for i in eachindex(xcm0), j in i+1:length(xcm0)
        abs(xcm1[i] - xcm1[j]) - 2R >= (1 - safety) * (abs(xcm0[i] - xcm0[j]) - 2r) && continue
        gap0 = boundary_gap(a, b, xcm0[i], θcm0[i], xcm0[j], θcm0[j])
        gap1 = boundary_gap(a, b, xcm1[i], θcm1[i], xcm1[j], θcm1[j])
        gap1 < (1 - safety) * gap0 && return false
    end
    return true
end

"""
    midpoint_step(velocity, xcm, θcm, v1, ω1, a, b;
                  safety=0.5, dt_min=1e-4, dt_max=1.0, dt_cap=Inf)

One explicit midpoint step with gap control.

The trial step is [`step_size`](@ref), capped at `dt_cap`. A step that lets any pair close
more than a fraction `safety` of its gap is rejected and retried with half the step, down
to `min(dt_min, dt_cap)`.

# Arguments
- `velocity`: function `(xcm, θcm) -> (v, ω)` giving the rigid-body velocities.
- `xcm`, `θcm`: current centers and orientations.
- `v1`, `ω1`: velocities at the current configuration.
- `a`, `b`: semi-axes of the particles.
- `safety`: fraction of a gap that may close in one step (default `0.5`).
- `dt_min`: smallest step tried before giving up (default `1e-4`).
- `dt_max`: largest allowed step (default `1.0`).
- `dt_cap`: additional cap on the step, e.g. the time to the next output (default `Inf`).

# Returns
- `(xcm, θcm, dt)` after the step, or `nothing` if no step of at least `min(dt_min, dt_cap)`
  keeps every gap.
"""
function midpoint_step(velocity, xcm, θcm, v1, ω1, a, b;
                       safety=0.5, dt_min=1e-4, dt_max=1.0, dt_cap=Inf)

    dt    = min(step_size(xcm, θcm, v1, ω1, a, b; safety=safety, dt_max=dt_max), dt_cap)
    floor = min(dt_min, dt_cap)

    while dt >= floor
        v2, ω2 = velocity(xcm .+ dt/2 .* v1, θcm .+ dt/2 .* ω1)
        xcm_new = xcm .+ dt .* v2
        θcm_new = θcm .+ dt .* ω2
        _gaps_kept(xcm, θcm, xcm_new, θcm_new, a, b, safety) && return xcm_new, θcm_new, dt
        dt /= 2
    end

    return nothing

end

"""
    suspension_simulate(xcm0, θcm0, Ω, Γ_D, Γ_N, g;
                        t_end=10.0, dt_min=1e-4, dt_max=1.0, safety=0.5, stop_gap=1e-2,
                        plot_every=10.0, plot_mode=:concentration,
                        output_file="", plot_dir="", tol=1e-6,
                        colors=fill(:white, length(xcm0)), cmap=:viridis, pad=4.0,
                        verbose=true, render=true, self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

Time-integrate a suspension of Marangoni surfers with the explicit midpoint method.

At each step, [`marangoni_solve`](@ref) gives the rigid-body velocity of every particle,
and positions and orientations are advanced by [`midpoint_step`](@ref), which limits the
step so that no pair closes more than a fraction `safety` of its gap. The run stops when the
smallest gap falls below `stop_gap`, or when no step of at least `dt_min` keeps the gaps.
Frames are displayed every
`plot_every` time units, and per-frame particle configurations are appended to
`frame_<label>.dat` files under `plot_dir`, if given.

# Arguments
- `xcm0`: initial complex particle centers.
- `θcm0`: initial particle orientations, in radians.
- `Ω`: Stokes domain boundary (shared by all particles).
- `Γ_D::Boundary`: reference Dirichlet arc boundary.
- `Γ_N::Boundary`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.
- `t_end`: total simulation time (default `1.0`).
- `dt_min`, `dt_max`: bounds on the adaptive time step (default `1e-4`, `1.0`).
- `safety`: fraction of a gap that may close in one step (default `0.5`).
- `stop_gap`: stop when the smallest boundary-to-boundary gap falls below this (default `1e-2`).
- `plot_every`: time interval between saved/displayed frames (default `1.0`).
- `plot_mode`: `:concentration` plots the surfactant field; `:distance` plots the
  two-particle separation instead (default `:concentration`).
- `output_file`: if non-empty, append `(t, separation)` pairs here when `plot_mode ==
  :distance` (default `""`, no file output).
- `plot_dir`: if non-empty, save configuration `.dat` files here (default `""`, no file
  output). Frames are displayed regardless.
- `tol`: FMM/GMRES accuracy tolerance passed to [`marangoni_solve`](@ref) (default `1e-6`).
- `colors`: per-particle fill colors (default `:white` for all particles).
- `cmap`: heatmap colormap (default `:viridis`).
- `pad`: margin added around the current particle centers when computing each frame's
  plot extent, so all particles stay in view as they move (default `4.0`).
- `verbose`: print the time and step size after each step (default `true`).
- `render`: render, display and save frames (default `true`); if `false`, only the
  configuration `.dat` files are written. Rendering requires Plots.jl (`using Plots`).
- `self_op`: self-flow operator from [`marangoni_self_operator`](@ref), built once and reused
  at every step.

# Returns
- `(traj_xcm, traj_θcm, traj_t)`: the recorded trajectory of particle centers,
  orientations, and times.
"""
function suspension_simulate(xcm0, θcm0, Ω, Γ_D::Boundary, Γ_N::Boundary, g;
                  t_end=10.0,
                  dt_min=1e-4, dt_max=1.0,
                  safety=0.5,
                  stop_gap=1e-2,
                  plot_every=10.0,
                  plot_mode=:concentration,
                  output_file="",
                  plot_dir="",
                  tol=1e-6,
                  colors=fill(:white, length(xcm0)),
                  cmap=:viridis,
                  pad=4.0,
                  verbose=true,
                  render=true,
                  self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

    render && !_plots_loaded() &&
        throw(ArgumentError("render=true needs Plots.jl: run `using Plots` first, or pass render=false"))

    xcm, θcm = copy(xcm0), copy(θcm0)

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

    a, b = real(Ω.f(1.0 + 0im)), imag(Ω.f(1.0im))
    velocity(x, θ) = rigid_velocities(first(marangoni_solve(x, θ, Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)))
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
        U0, σ0, Γ0, ϕinf0 = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)
        !isempty(plot_dir) && mkpath(plot_dir)
        if render
            xlims0, ylims0 = plot_lims(xcm)
            p0 = plot_concentration(xcm, σ0, ϕinf0, Γ0; title="t = 0.00", colors=colors, cmap=cmap, xlims=xlims0, ylims=ylims0)
            _show_frame(p0, isempty(plot_dir) ? "" : joinpath(plot_dir, "frame_000000.png"))
        end
        save_config("000000", 0.0, xcm, θcm, U0)
    end

    while t < t_end

        U1, σ, Γ, ϕinf = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)
        v1, ω1 = rigid_velocities(U1)

        step = midpoint_step(velocity, xcm, θcm, v1, ω1, a, b; safety=safety, dt_min=dt_min,
                             dt_max=dt_max, dt_cap=min(t_end - t, next_plot - t))
        if step === nothing
            @printf("stopped at t = %.4f: no step of at least dt_min keeps the gaps\n", t)
            break
        end

        xcm_plot = copy(xcm)
        xcm, θcm, dt = step
        t += dt
        nt += 1
        push!(traj_xcm, copy(xcm))
        push!(traj_θcm, copy(θcm))
        push!(traj_t,   t)

        verbose && @printf("t = %10.4f   dt = %.2e\n", t, dt)

        frame = lpad(length(traj_t), 6, '0')
        if t >= next_plot
            if plot_mode == :distance
                dists = [abs(traj_xcm[k][2] - traj_xcm[k][1]) for k in eachindex(traj_t)]
                if !isempty(output_file)
                    open(output_file, "a") do io
                        @printf(io, "%.6f %.6f\n", traj_t[end], dists[end])
                    end
                end
            end
            !isempty(plot_dir) && mkpath(plot_dir)
            if render
                if plot_mode == :distance
                    p = _plot_separation(traj_t, dists .- 2.0)
                else
                    xlims, ylims = plot_lims(xcm_plot)
                    p = plot_concentration(xcm_plot, σ, ϕinf, Γ;
                                           title="t = $(round(t, digits=2))", colors=colors, cmap=cmap, xlims=xlims, ylims=ylims)
                end
                _show_frame(p, isempty(plot_dir) ? "" : joinpath(plot_dir, "frame_$(frame).png"))
            end
            save_config(frame, t, xcm, θcm, U1)
            next_plot += plot_every
        end

        gap = min_gap(xcm, θcm, a, b)
        if gap < stop_gap
            @printf("stopped at t = %.4f: gap %.2e < %.0e\n", t, gap, stop_gap)
            break
        end
    end

    return traj_xcm, traj_θcm, traj_t

end
