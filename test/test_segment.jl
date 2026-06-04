# Exact solution for -1/(2π) ∫₋₁¹ log|x-y| α(y) / √(1-y²) dy = λₙTₙ(y), x ∈ [-1,1]
# Smooth density: σ_exact = T_n
L = 1.0
println("Segment: single-layer solve convergence (σ = T₂(s), L = 2)")

σ_exact(s) = 2s^2 - 1
λ = L/4

for N in Nspan
    Γ = segment(N; dom=[-L, L])
    g = λ * σ_exact.(Γ.s)
    σ, _ = laplace_solve([Γ], [:D], [g])
    err = norm(σ[1] .- σ_exact.(Γ.s), Inf)
    printerror(N, err)
end
