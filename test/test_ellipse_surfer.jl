
using SpecialFunctions: ellipk, ellipe
using Plots

# Ellipse surfer: convergence against the semi-analytical velocity formula
#   u_{cm,i} = (m_i - m_0) / (2π a b) * ∫_Γ xᵢ σ dΓ
# where m_0 = b/2 K(ε), m_1 = b/(2ε) [K(ε) - E(ε)], ε = sqrt(1 - b²/a²).
# Only the Laplace BIM (arc quadrature) contributes to the error.
println("\nEllipse surfer (a=2.0, b=1.0): convergence vs semi-analytical solution")
let
    α = π/2
    g(s) = cos(2s)
    a, b = 2.0, 1.0
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)
    M  = 32
    Ω  = ellipse(a, b, M)

    ε  = sqrt(1 - b^2/a^2)
    m0 = b/2 * ellipk(ε^2)
    m1 = b/(2ε) * (ellipk(ε^2) - ellipe(ε^2))
    mobility = (m1 - m0) / (2π * a * b)

    semi_analytical(σ, Γ) = mobility * (integral(real.(Γ[1].γ) .* σ[1], Γ[1]) +
                                        integral(real.(Γ[2].γ) .* σ[2], Γ[2]))

    Ns_ref = 1024
    Γ_D_ref = arc(γ, γs, γss, Ns_ref; dom=[α, 2π-α])
    Γ_N_ref = arc(γ, γs, γss, Ns_ref; dom=[-α, α])
    _, σ_ref, Γ_ref, _ = marangoni_solve(Ω, Γ_D_ref, Γ_N_ref, g)
    Uref = semi_analytical(σ_ref, Γ_ref)
    
    for Ns in Nspan
        Γ_D = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
        Γ_N = arc(γ, γs, γss, Ns; dom=[-α, α])
        _, σ, Γ, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
        Uexact = semi_analytical(σ, Γ)
        err = abs(Uexact - Uref)
        printerror(Ns, err)
    end
end
