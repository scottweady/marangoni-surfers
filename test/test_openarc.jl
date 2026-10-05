# Dirichlet solve on an open arc (part of a star-shaped curve). Convergence is measured by
# the potential at an off-curve point against a high-N reference.

println("\nOpen arc: Dirichlet solve, self-refinement")
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
