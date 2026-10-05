# Circular surfer with g(s) = cos(2s): convergence in the arc resolution `Ns` against the
# exact swimming speed (π/16)(1 + cos s₀) sin² s₀.

println("\nCircular surfer: convergence in Ns against exact speed")
let
    s₀ = π/2
    g(s) = cos(2s)
    Uexact = (π/16) * (1 + cos(s₀)) * sin(s₀)^2
    M = 32
    aspect_ratio = 1.0
    a = sqrt(aspect_ratio)
    b = 1.0 / a
    Ω_ref = ellipse(a, b, M)
    for Ns in Nspan
        Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)
        U, _ = marangoni_solve(Ω_ref, Γ_D, Γ_N, g)
        err = abs(U[1] - Uexact)
        printerror(Ns, err)
    end
end
