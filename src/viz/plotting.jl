
using Plots

function _boundary_polygon(Γ_D::Boundary, Γ_N::Boundary)
    pts = vcat(Γ_N.γ, Γ_D.γ)[sortperm(vcat(Γ_N.s, Γ_D.s))]
    return real.(pts), imag.(pts)
end

function plot_concentration(xcm, σ, ϕinf, Γ; pad=4.0, Nplt=150, title="")

    Nparticles = length(xcm)
    xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
    yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
    z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

    p = heatmap(xplt, yplt, ϕ, aspect_ratio=:equal, colorbar=true, colorbar_title="ϕ", title=title, legend=false, clims=(0, 1))
    for k in 1:Nparticles
        xs, ys = _boundary_polygon(Γ[2k-1], Γ[2k])
        plot!(p, Shape(xs, ys), lc=:black, lw=1.5, fc=:white, fa=1)
    end
    return p

end
