"""
    MarangoniSurfersPlotsExt

Plots.jl backend for MarangoniSurfers, loaded automatically with `using Plots`: methods for
[`plot_concentration`](@ref) and for the frame rendering of [`suspension_simulate`](@ref).
"""
module MarangoniSurfersPlotsExt

using MarangoniSurfers
using MarangoniSurfers: _boundary_polygon
using Plots

function MarangoniSurfers.plot_concentration(xcm, σ, ϕinf, Γ; pad=4.0, Nplt=150, title="", colors=fill(:white, length(xcm)), cmap=:viridis, xlims=nothing, ylims=nothing)

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

MarangoniSurfers._plot_separation(t, d) =
    plot(t, d, xlabel="t", ylabel="|r₂ − r₁|", legend=false, title="particle separation")

function MarangoniSurfers._show_frame(p, path)
    display(p)
    isempty(path) || savefig(p, path)
    return nothing
end

end
