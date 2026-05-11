
function plot_concentration(centers, σs, boundaries, ϕinf; pad=2.5, Nplt=150, title="")
    Np   = length(centers)
    xplt = range(minimum(real.(centers)) - pad, maximum(real.(centers)) + pad, length=Nplt)
    yplt = range(minimum(imag.(centers)) - pad, maximum(imag.(centers)) + pad, length=Nplt)
    ζ    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    ϕ = sum(𝒱(ζ, σs[2k-1], boundaries[2k-1]) .+ 𝒱(ζ, σs[2k], boundaries[2k]) for k in 1:Np) .+ ϕinf
    for z0 in centers
        ϕ[abs.(ζ .- z0) .< 1.0] .= NaN
    end

    p = heatmap(xplt, yplt, ϕ, aspect_ratio=:equal, colorbar_title="ϕ", title=title, legend=false)
    θ = range(0, 2π, length=64)
    for i in 1:Np
        plot!(p, real(centers[i]) .+ cos.(θ), imag(centers[i]) .+ sin.(θ), label=false, lc=:black, lw=2)
    end
    return p
end
