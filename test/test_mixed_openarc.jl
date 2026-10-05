# Mixed Dirichlet/Neumann solve on an ellipse split into two arcs, with smooth Dirichlet
# data on Γ_D and zero Neumann data on Γ_N. Convergence is measured by the potential at
# an off-curve point against a high-N reference.

println("\nMixed arcs: Dirichlet/Neumann solve, self-refinement")
let
    s₀ = π/2
    a, b = 2.0, 1.0
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    ζ = 3.0 + 0.0im

    function solve_and_eval(N)
        Γ_D = openarc(γ, γs, γss, N; dom=[s₀, 2π-s₀])
        Γ_N = openarc(γ, γs, γss, N; dom=[-s₀, s₀])
        σ, ϕ∞ = laplace_solve([Γ_D, Γ_N], [:D, :N], [-cos.(3.21*Γ_D.s).*exp.(-Γ_D.s.^2), zeros(N)])
        return ϕ∞ + 𝒱(ζ, σ[1], Γ_D) + 𝒱(ζ, σ[2], Γ_N)
    end

    Nref = 500
    uref = solve_and_eval(Nref)

    for N in Nspan
        u = solve_and_eval(N)
        err = abs(u - uref)
        printerror(N, err)
    end
end
