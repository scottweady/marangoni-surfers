
function _boundary_polygon(Γ_D::Boundary, Γ_N::Boundary)
    idx = sortperm(vcat(Γ_N.s, Γ_D.s))
    pts = vcat(Γ_N.γ, Γ_D.γ)[idx]
    return real.(pts), imag.(pts)
end

function plot_concentration(centers, σs, boundaries, ϕinf; pad=4.0, Nplt=150, title="")
    Np   = length(centers)
    # xplt = range(minimum(real.(centers)) - pad, maximum(real.(centers)) + pad, length=Nplt)
    # yplt = range(minimum(imag.(centers)) - pad, maximum(imag.(centers)) + pad, length=Nplt)
    xplt = range(-pad, pad, length=Nplt)
    yplt = range(-pad, pad, length=Nplt)
    ζ    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    ϕ = sum(𝒱(ζ, σs[2k-1], boundaries[2k-1]) .+ 𝒱(ζ, σs[2k], boundaries[2k]) for k in 1:Np) .+ ϕinf

    p = heatmap(xplt, yplt, ϕ, aspect_ratio=:equal, colorbar_title="ϕ", title=title, legend=false)
    for k in 1:Np
        xs, ys = _boundary_polygon(boundaries[2k-1], boundaries[2k])
        plot!(p, Shape(xs, ys), lc=:black, lw=1.5, fc=:white, fa=1)
    end
    return p
end
