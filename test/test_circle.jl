
# K′[1] = 1/2 exactly on any closed curve (interior jump condition)
println("\nCircle: adjoint double-layer applied to 1 (exact = 0.5)")
for N in Nspan
    Γ   = circle(N)
    err = norm(Γ.V′ * ones(N) .+ 0.5, Inf)
    printerror(N, err)
end

# Manufactured density σ = cos(s), solve V⁻¹(V σ) and recover σ
println("\nCircle: single-layer solve convergence (σ = cos(s))")
for N in Nspan
    Γ = circle(N)
    σ_exact = 2 * cos.(Γ.s)
    g = cos.(Γ.s)
    σ, _ = V⁻¹(g, Γ)
    err = norm(σ .- σ_exact, Inf)
    printerror(N, err)
end
