
function _rk_increment(U, dt, Np)
    Δc = dt .* [U[i][1] + im*U[i][2] for i in 1:Np]
    Δψ = dt .* [U[i][3]               for i in 1:Np]
    return Δc, Δψ
end

function simulate(centers0, orientations0;
                  t_end=10.0,
                  dt_min=1e-4, dt_max=1e3,
                  tol=0.05, safety=0.9,
                  plot_every=1.0,
                  plot_mode=:concentration,
                  g=g, α=α, Nₛ=Nₛ, D=D, eps=1e-6)

    centers      = copy(centers0)
    orientations = copy(orientations0)
    Np           = length(centers)
    t            = 0.0
    next_plot    = plot_every

    traj_pos  = [copy(centers)]
    traj_ang  = [copy(orientations)]
    traj_t    = [t]

    step = 0
    while t < t_end
        U1, σs, boundaries, ϕinf = marangoni_step(centers, orientations; g=g, α=α, Nₛ=Nₛ, D=D, eps=eps)

        max_speed = maximum(hypot(U1[i][1], U1[i][2]) for i in 1:Np)
        dt = clamp(safety * tol / max(max_speed, 1e-14), dt_min, dt_max)
        dt = min(dt, t_end - t)

        Δc1, Δψ1 = _rk_increment(U1, dt, Np)

        c_pred  = centers      .+ Δc1
        ψ_pred  = orientations .+ Δψ1
        U2, _, _, _ = marangoni_step(c_pred, ψ_pred; g=g, α=α, Nₛ=Nₛ, D=D, eps=eps)
        Δc2, Δψ2 = _rk_increment(U2, dt, Np)

        centers_plot  = copy(centers)
        centers      .+= 0.5 .* (Δc1  .+ Δc2)
        orientations .+= 0.5 .* (Δψ1  .+ Δψ2)
        t    += dt
        step += 1
        push!(traj_pos, copy(centers))
        push!(traj_ang, copy(orientations))
        push!(traj_t,   t)

        @printf("  t = %.4f  dt = %.2e  speed = %.2e", t, dt, max_speed)
        for i in 1:Np
            @printf("   U%d=(%.4f,%.4f)", i, U1[i][1], U1[i][2])
        end
        println()

        if t >= next_plot
            if plot_mode == :distance
                dists = [abs(traj_pos[k][2] - traj_pos[k][1]) for k in eachindex(traj_t)]
                p = plot(traj_t, dists .- 2.0, xlabel="t", ylabel="|r₂ − r₁|",
                         legend=false, title="particle separation")
            else
                p = plot_concentration(centers_plot, σs, boundaries, ϕinf;
                                       title="t = $(round(t, digits=2))")
            end
            display(p)
            next_plot += plot_every
        end
    end

    return traj_pos, traj_ang, traj_t
end
