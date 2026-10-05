# Elliptical surfer: convergence in the Stokes expansion order `M`, with the arcs fixed at
# `Ns_ref`, against the semi-analytical velocity
#
#     u_{cm,i} = (m_i - m_0) / (2π a b) ∫_Γ xᵢ σ dΓ,
#
# where m_0 = b/2 K(ε²), m_1 = b/(2ε²) [K(ε²) - E(ε²)], ε = √(1 - b²/a²).

using SpecialFunctions: ellipk, ellipe

println("\nElliptical surfer: convergence in M against semi-analytical speed")
let
    s₀ = π/2
    g(s) = -cos(s)
    aspect_ratio = 4.0
    a = sqrt(aspect_ratio)
    b = 1.0 / a

    ε  = sqrt(1 - b^2/a^2)
    m0 = b/2 * ellipk(ε^2)
    m1 = b/(2ε^2) * (ellipk(ε^2) - ellipe(ε^2))
    mobility = (m1 - m0) / (2π * a * b)

    semi_analytical(σ, Γ) = mobility * (integral(real.(Γ[1].γ) .* σ[1], Γ[1]) +
                                        integral(real.(Γ[2].γ) .* σ[2], Γ[2]))

    Ns_ref = 1024
    M_ref = 48

    Ω_ref, Γ_D_ref, Γ_N_ref = surfer(a, b, s₀, Ns_ref, M_ref)
    _, σ_ref, Γ_ref, _ = marangoni_solve(Ω_ref, Γ_D_ref, Γ_N_ref, g)
    Uref = semi_analytical(σ_ref, Γ_ref)

    for M in Nspan
        Ω = ellipse(a, b, M)
        U, _ = marangoni_solve(Ω, Γ_D_ref, Γ_N_ref, g)
        err = abs(U[1] - Uref)
        printerror(M, err)
    end
end
