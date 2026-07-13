
# ── Open arc (star) Dirichlet self-refinement study ──────────────────────────
# Solve ∫_Γ G(ζ,γ') σ(γ') ds' + ϕ∞ = g on the upper half of the star curve.
# Convergence measured via the potential at an off-curve point against a high-N reference.
println("\nOpen arc (star (1+0.3·cos5s)eⁱˢ, s∈[0,π]): laplace_solve Dirichlet, self-refinement")
let
    a = 0.3
    γ(s)   = (1 + a * cos(5s)) * exp(im * s)
    γs(s)  = im * (1 + a * cos(5s)) * exp(im * s) - 5a * sin(5s) * exp(im * s)
    γss(s) = -(1 + a * cos(5s)) * exp(im * s) - 10a * im * sin(5s) * exp(im * s) - 25a * cos(5s) * exp(im * s)

    ζ = 3.0 + 1.0im

    function solve_and_eval(N)
        Γ = openarc(γ, γs, γss, N; dom=[0.0, π])
        g = real.(Γ.γ)
        σs, ϕ∞ = laplace_solve([Γ], [:D], [g])
        return 𝒱(ζ, σs[1], Γ) + ϕ∞
    end

    Nref = 500
    uref = solve_and_eval(Nref)

    for N in Nspan
        u = solve_and_eval(N)
        err = abs(u - uref)
        printerror(N, err)
    end
end
