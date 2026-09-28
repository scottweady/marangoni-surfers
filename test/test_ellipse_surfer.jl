
using SpecialFunctions: ellipk, ellipe

# Ellipse surfer: convergence against the semi-analytical velocity formula
#   u_{cm,i} = (m_i - m_0) / (2π a b) * ∫_Γ xᵢ σ dΓ
# where m_0 = b/2 K(ε), m_1 = b/(2ε) [K(ε) - E(ε)], ε = sqrt(1 - b²/a²).
# Only the Laplace BIM (arc quadrature) contributes to the error.
println("\nEllipse surfer (a=2.0, b=1.0): convergence vs semi-analytical solution")
let
    α = π/2
    g(s) = -cos(s)
    a = sqrt(4.0)
    b = 1.0/a
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    ε  = sqrt(1 - b^2/a^2)
    m0 = b/2 * ellipk(ε^2)
    m1 = b/(2ε^2) * (ellipk(ε^2) - ellipe(ε^2))
    mobility = (m1 - m0) / (2π * a * b)

    semi_analytical(σ, Γ) = mobility * (integral(real.(Γ[1].γ) .* σ[1], Γ[1]) +
                                        integral(real.(Γ[2].γ) .* σ[2], Γ[2]))

    Ns_ref = 1024
    M_ref = 48

    Γ_D_ref = openarc(γ, γs, γss, Ns_ref; dom=[α, 2π-α])
    Γ_N_ref = openarc(γ, γs, γss, Ns_ref; dom=[-α, α], quad=:open)
    Ω_ref = ellipse(a, b, M_ref)
    _, σ_ref, Γ_ref, _ = marangoni_solve(Ω_ref, Γ_D_ref, Γ_N_ref, g)
    Uref = semi_analytical(σ_ref, Γ_ref)
    
    for Ns in Nspan
        Ω = ellipse(a, b, Ns)
        Γ_D = openarc(γ, γs, γss, Ns_ref; dom=[α, 2π-α])
        Γ_N = openarc(γ, γs, γss, Ns_ref; dom=[-α, α], quad=:open)
        U, σ, Γ, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
        err = abs(U[1] - Uref)
        printerror(Ns, err)
    end
end
