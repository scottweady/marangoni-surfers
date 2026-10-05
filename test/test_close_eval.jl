# Close evaluation of the single layer (𝒱_close) and adjoint double layer (K′_close) on
# an open star-shaped arc, at a target just off the arc along the outward normal.
# Compared against direct quadrature on a high-N reference arc.

println("\nOpen arc: close evaluation of 𝒱 and K′")
let
    a      = 0.3
    γ(s)   = (1 + a*cos(5s)) * exp(im*s)
    γs(s)  = im*(1 + a*cos(5s))*exp(im*s) - 5a*sin(5s)*exp(im*s)
    γss(s) = -(1 + a*cos(5s))*exp(im*s) - 10a*im*sin(5s)*exp(im*s) - 25a*cos(5s)*exp(im*s)
    dom    = [-π/2, π/2]
    σ_fun(s) = exp(cos(s)) + sin(5s) + 0.1

    s₀  = 0.0
    n̂   = (imag(γs(s₀)) - im*real(γs(s₀))) / abs(γs(s₀))   # tangent rotated by -90°
    ζ   = γ(s₀) + 0.01 * n̂

    Γref = openarc(γ, γs, γss, 2^12; dom=dom)
    σref = σ_fun.(Γref.s)
    ref_V  = 𝒱(ζ, σref, Γref)
    ref_K′ = -(1/(2π)) * sum(dot.(Ref(n̂), ζ .- Γref.γ) ./ abs2.(ζ .- Γref.γ) .* σref .* Γref.dΓ)

    println("  𝒱_close (single layer)")
    for N in Nspan
        Γ   = openarc(γ, γs, γss, N; dom=dom)
        σ   = σ_fun.(Γ.s)
        val = 𝒱_close(ζ, σ, Γ)
        printerror(N, abs(val - ref_V))
    end

    println("  K′_close (adjoint double layer)")
    for N in Nspan
        Γ   = openarc(γ, γs, γss, N; dom=dom)
        σ   = σ_fun.(Γ.s)
        val = K′_close(ζ, n̂, σ, Γ)
        printerror(N, abs(val - ref_K′))
    end
end
