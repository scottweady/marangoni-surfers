# Single-layer Dirichlet solve on a straight segment [-L, L]. The weighted single layer
# maps Chebyshev polynomials to multiples of themselves,
#
#     -1/(2π) ∫₋₁¹ log|x - y| Tₙ(y) / √(1 - y²) dy = Tₙ(x) / (2n),   n ≥ 1,
#
# so Dirichlet data L T₂ / 4 has exact density T₂, recovered at every resolution.

@testset "Segment: single-layer Dirichlet solve" begin
    section("Segment: single-layer Dirichlet solve, refining boundary resolution N")
    Nspan = 2:6
    L = 1.0
    σ_exact(s) = 2s^2 - 1
    μ = L/4

    errs = map(Nspan) do N
        Γ = segment(N; dom=[-L, L])
        σ, _ = laplace_solve([Γ], [:D], [μ * σ_exact.(Γ.s)])
        norm(σ[1] .- σ_exact.(Γ.s), Inf)
    end
    check_errors(Nspan, errs; every=true)
end
