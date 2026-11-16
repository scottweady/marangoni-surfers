
"""
    _boundary_polygon(Γ_D, Γ_N)

Merge a particle's Dirichlet and Neumann arc boundaries into a single closed polygon,
ordered by arc parameter.

# Arguments
- `Γ_D::Boundary`: Dirichlet arc boundary.
- `Γ_N::Boundary`: Neumann arc boundary.

# Returns
- `(x, y)`: real and imaginary coordinates of the merged, arc-parameter-sorted polygon.
"""
function _boundary_polygon(Γ_D::Boundary, Γ_N::Boundary)
    pts = vcat(Γ_N.γ, Γ_D.γ)[sortperm(vcat(Γ_N.s, Γ_D.s))]
    return real.(pts), imag.(pts)
end

"""
    plot_concentration(xcm, σ, ϕinf, Γ; pad=4.0, Nplt=150, title="")

Heatmap of the surfactant concentration field around a collection of particles, with
each particle boundary drawn as a filled white polygon.

# Arguments
- `xcm`: complex particle centers, used to set the plot extent.
- `σ`: layer densities, ordered `[Γ_D₁, Γ_N₁, Γ_D₂, Γ_N₂, ...]` (as returned by
  [`surfactant_solve`](@ref)/[`marangoni_solve`](@ref)).
- `ϕinf`: far-field concentration constant.
- `Γ`: placed boundaries, ordered `[Γ_D₁, Γ_N₁, Γ_D₂, Γ_N₂, ...]`.
- `pad`: margin added around the particle centers when framing the plot (default `4.0`).
- `Nplt`: number of grid points per axis (default `150`).
- `title`: plot title (default `""`).
- `colors`: per-particle fill colors (default `:white` for all particles).
- `cmap`: heatmap colormap (default `:viridis`).
- `xlims`, `ylims`: fixed plot extents; if `nothing`, computed from `xcm` (default `nothing`).

# Returns
- The `Plots.jl` plot object.
"""
function plot_concentration(xcm, σ, ϕinf, Γ; pad=4.0, Nplt=150, title="", colors=fill(:white, length(xcm)), cmap=:viridis, xlims=nothing, ylims=nothing)

    Nparticles = length(xcm)
    xplt = range(something(xlims, (minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad))..., length=Nplt)
    yplt = range(something(ylims, (minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad))..., length=Nplt)
    z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

    p = heatmap(xplt, yplt, ϕ, aspect_ratio=:equal, colorbar=true, colorbar_title="ϕ", title=title, legend=false, clims=(0, 1), c=cmap)
    for k in 1:Nparticles
        xs, ys = _boundary_polygon(Γ[2k-1], Γ[2k])
        plot!(p, Shape(xs, ys), lc=:black, lw=1.5, fc=colors[k], fa=1)
    end
    return p

end

"""
    concentration_grid(xg, yg, xcm, θcm, σ, ϕinf, Γ, a, b; R=3max(a, b), P=20)

Surfactant concentration on the grid `xg × yg` around elliptical particles: the single
layers of each particle summed directly within `R` of its center and through their
multipole expansion about the center, to order `P`, beyond it. Threaded over `xg`.

# Arguments
- `xg`, `yg`: grid coordinates.
- `xcm`, `θcm`: complex particle centers and orientations.
- `σ`, `ϕinf`, `Γ`: layer densities, far-field constant and placed boundaries, as returned by
  [`surfactant_solve`](@ref).
- `a`, `b`: semi-axes along and across the particle orientation.

# Returns
- `ϕ`: `length(yg) × length(xg)` matrix, `NaN` inside the particles.
"""
function concentration_grid(xg, yg, xcm, θcm, σ, ϕinf, Γ, a, b; R=3max(a, b), P=20)
    N = length(xcm)
    q = [vcat(σ[2k-1] .* Γ[2k-1].dΓ, σ[2k] .* Γ[2k].dΓ) for k in 1:N]   # weighted densities
    y = [vcat(Γ[2k-1].γ, Γ[2k].γ) for k in 1:N]                        # nodes
    M = [[sum(q[k] .* (y[k] .- xcm[k]).^p) for p in 0:P] for k in 1:N]  # moments about the center

    ϕ = fill(ϕinf, length(yg), length(xg))
    Threads.@threads for ix in eachindex(xg)
        for iy in eachindex(yg)
            z = xg[ix] + im * yg[iy]
            u = 0.0
            for k in 1:N
                ζ = z - xcm[k]
                if abs(ζ) < R
                    zb = exp(-im * θcm[k]) * ζ
                    if (real(zb) / a)^2 + (imag(zb) / b)^2 < 1
                        u = NaN
                        break
                    end
                    u += sum(q[k][j] * log(abs(z - y[k][j])) for j in eachindex(q[k]))
                else
                    w = M[k][1] * log(ζ)
                    ζp = one(ζ)
                    for p in 1:P
                        ζp *= ζ
                        w -= M[k][p+1] / (p * ζp)
                    end
                    u += real(w)
                end
            end
            ϕ[iy, ix] += -u / (2π)
        end
    end
    return ϕ
end
