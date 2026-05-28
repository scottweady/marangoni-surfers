
# Exact solution for -1/(2π) ∫₋₁¹ log|x-y| α(y) / √(1-y²) dy = λₙTₙ(y), x ∈ [-1,1]
# Smooth density: α_exact = T_n
L = 1.0
println("Segment: single-layer solve convergence (α = T₂(s), L = 2)")

α_exact(s) = 2s^2 - 1
λ = L/4

for N in Nspan
    Γ = segment(N; dom=[-L, L])
    g = λ * α_exact.(Γ.s)
    α, _ = V⁻¹(g, Γ; decay=true)
    err = norm(α .- α_exact.(Γ.s), Inf)
    printerror(N, err)
end
