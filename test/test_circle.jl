# Closed-curve operators on the unit circle, where both checks have exact answers.

# On the circle K′ coincides with the double layer, so K′ 1 = -1/2
println("\nCircle: adjoint double layer applied to a constant")
for N in Nspan
    Γ   = circle(N)
    err = norm(Γ.K′ * ones(N) .+ 0.5, Inf)
    printerror(N, err)
end

# Dirichlet data g = cos(s) has exact single-layer density σ = 2 cos(s)
println("\nCircle: single-layer Dirichlet solve")
for N in Nspan
    Γ = circle(N)
    σ_exact = 2 * cos.(Γ.s)
    g = cos.(Γ.s)
    σs, _ = laplace_solve([Γ], [:D], [g])
    err = norm(σs[1] .- σ_exact, Inf)
    printerror(N, err)
end
