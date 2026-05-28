
# ── Curve (ellipse) Dirichlet test ───────────────────────────────────────────
# Solve ∫_Γ G(ζ,γ') σ(γ') ds' + ϕ∞ = g on an ellipse using the `curve` constructor.
# Exact solution: u(x,y) = x  (harmonic everywhere).
# After solving, evaluate 𝒱(ζ, σ, Γ) + ϕ∞ at interior points and compare to Re(ζ).
println("\nCurve (star γ=(1+a·cos5s)eⁱˢ, a=0.3): laplace_solve Dirichlet solve, u=x")
let
    a = 0.3
    γ(s)   = (1 + a * cos(5s)) * exp(im*s)
    γs(s)  = im * (1 + a * cos(5s))*exp(im*s) - 5a * sin(5s) * exp(im*s)
    γss(s) = -(1 + a * cos(5s)) * exp(im*s) - 10a * im * sin(5s) * exp(im*s) - 25a * cos(5s) * exp(im*s)

    ζ_pts = [0.0 + 0.0im, 0.2 + 0.1im, -0.3 + 0.2im, 0.1 - 0.2im]
    u_exact = real.(ζ_pts)

    println("  interior potential")
    for N in Nspan
        Γ = curve(γ, γs, γss, N)
        g = real.(Γ.γ)
        σs, ϕ∞ = laplace_solve([Γ], [:D], [g])
        σ = σs[1]
        u_num = [𝒱(ζ, σ, Γ) + ϕ∞ for ζ in ζ_pts]
        err = norm(u_num .- u_exact, Inf)
        printerror(N, err)
    end

    # Interior normal derivative: ∂u/∂n = -σ/2 + K′σ
    # Exact: ∂(Re ζ)/∂n = Re(outward normal) = Re(i γt/|γt|) = -Im(γt)/|γt|
    println("  boundary normal derivative")
    for N in Nspan
        Γ = curve(γ, γs, γss, N)
        g = real.(Γ.γ)
        σs, _ = laplace_solve([Γ], [:D], [g])
        σ = σs[1]
        dudn_num   = σ / 2 + V′(σ, Γ)
        dudn_exact = imag.(Γ.γt) ./ abs.(Γ.γt)
        err = norm(dudn_num .- dudn_exact, Inf)
        printerror(N, err)
    end
end
