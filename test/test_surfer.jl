
# Exact swimming speed for a unit disk with g(θ) = cos(2θ), α = π/2.
println("\nDisk surfer: convergence in Ns")
let
    α = π/2
    g(θ) = cos(2θ)
    Uexact = (π/16) * (1 + cos(α)) * sin(α)^2
    M = 32
    a, b = 1.0, 1.0
    γ(s)   =  a * cos.(s) + im * b * sin.(s)
    γs(s)  = -a * sin.(s) + im * b * cos.(s)
    γss(s) = -a * cos.(s) - im * b * sin.(s)
    Ω_ref = ellipse(a, b, M)
    for Ns in Nspan
        Γ_D = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
        Γ_N = arc(γ, γs, γss, Ns; dom=[-α, α])
        U, _ = marangoni_solve(Ω_ref, Γ_D, Γ_N, g)
        err = abs(U[1] - Uexact)
        printerror(Ns, err)
    end
end
