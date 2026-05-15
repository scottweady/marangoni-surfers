
function simulate(centers0, orientations0, Ω, Γ_D::Boundary, Γ_N::Boundary, g;
                  t_end=10.0,
                  dt_min=1e-4, dt_max=1.0,
                  tol=0.5, safety=0.9,
                  plot_every=1.0,
                  plot_mode=:concentration,
                  output_file="",
                  eps=1e-6)

    centers      = copy(centers0)
    orientations = copy(orientations0)
    Np           = length(centers)
    t            = 0.0
    next_plot    = plot_every

    traj_pos  = [copy(centers)]
    traj_ang  = [copy(orientations)]
    traj_t    = [t]

    dt   = dt_min
    step = 0
    while t < t_end
        dt = min(dt, t_end - t)

        U1, σs, boundaries, ϕinf = marangoni_solve(centers, orientations, Ω, Γ_D, Γ_N, g; eps=eps)
        Δc1 = dt .* [U1[i][1] + im*U1[i][2] for i in 1:Np]
        Δψ1 = dt .* [U1[i][3]                for i in 1:Np]

        U2, _, _, _ = marangoni_solve(centers .+ Δc1, orientations .+ Δψ1, Ω, Γ_D, Γ_N, g; eps=eps)
        Δc2 = dt .* [U2[i][1] + im*U2[i][2] for i in 1:Np]
        Δψ2 = dt .* [U2[i][3]                for i in 1:Np]

        function edge_gap(i, j)
            d     = centers[j] - centers[i]
            d_hat = d / abs(d)
            pi    = vcat(boundaries[2i-1].γ, boundaries[2i].γ)
            pj    = vcat(boundaries[2j-1].γ, boundaries[2j].γ)
            xi    = pi[argmax(real.(conj(d_hat) .* (pi .- centers[i])))]
            xj    = pj[argmax(real.(conj(d_hat) .* (centers[j] .- pj)))]
            return abs(xi - xj)
        end
        d_min     = Np > 1 ? minimum(edge_gap(i, j) for i in 1:Np for j in i+1:Np) : Inf
        max_speed = maximum(hypot(U1[i][1], U1[i][2]) for i in 1:Np)
        dt_next   = clamp(safety * tol * d_min / max(max_speed, 1e-14), dt_min, dt_max)

        centers_plot  = copy(centers)
        centers      .+= 0.5 .* (Δc1  .+ Δc2)
        orientations .+= 0.5 .* (Δψ1  .+ Δψ2)
        t    += dt
        dt    = dt_next
        step += 1
        push!(traj_pos, copy(centers))
        push!(traj_ang, copy(orientations))
        push!(traj_t,   t)

        @printf("  t = %.4f  dt = %.2e  d_min = %.3f  speed = %.2e", t, dt_next, d_min, max_speed)
        for i in 1:Np
            @printf("   U%d=(%.4f,%.4f)", i, U1[i][1], U1[i][2])
        end
        println()

        if t >= next_plot
            if plot_mode == :distance
                dists = [abs(traj_pos[k][2] - traj_pos[k][1]) for k in eachindex(traj_t)]
                if !isempty(output_file)
                    open(output_file, "a") do io
                        @printf(io, "%.6f %.6f\n", traj_t[end], dists[end])
                    end
                end
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
