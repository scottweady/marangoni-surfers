# Single-layer Dirichlet solve on a straight segment [-L, L]. The weighted single layer
# maps Chebyshev polynomials to multiples of themselves,
#
#     -1/(2π) ∫₋₁¹ log|x - y| Tₙ(y) / √(1 - y²) dy = Tₙ(x) / (2n),   n ≥ 1,
#
# so Dirichlet data L T₂ / 4 has exact density T₂.

println("\nSegment: single-layer Dirichlet solve")
let
    L = 1.0
    σ_exact(s) = 2s^2 - 1
    λ = L/4

    for N in Nspan
        Γ = segment(N; dom=[-L, L])
        g = λ * σ_exact.(Γ.s)
        σ, _ = laplace_solve([Γ], [:D], [g])
        err = norm(σ[1] .- σ_exact.(Γ.s), Inf)
        printerror(N, err)
    end
end
