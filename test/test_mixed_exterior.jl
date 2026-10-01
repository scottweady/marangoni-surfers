
# ── Mixed exterior problem (ellipse) ─────────────────────────────────────────
# Exact exterior solution u = Re Σ cₖ/(z - zₖ), poles inside the ellipse, so u → 0
# (ϕ∞ = 0). Dirichlet data u on Γ_D, exterior Neumann data ∂u/∂n on Γ_N.
# Catches a wrong-sided jump term: the error is O(1) and does not converge.
println("\nMixed exterior (ellipse a=2,b=1): manufactured solution, potential at far points")
let
    a, b = 2.0, 1.0
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    poles = [(0.7 + 0.2im, 1.0 + 0.5im), (-0.9 - 0.1im, -0.6 + 0.3im)]
    f(z)  = sum(c / (z - z0) for (z0, c) in poles)
    fp(z) = sum(-c / (z - z0)^2 for (z0, c) in poles)

    ζ_pts = [3.0 + 0.0im, 0.0 + 2.0im, -2.5 - 1.5im, 1.0 - 2.5im]

    err = Inf
    for N in (16, 32, 64, 128, 256)
        Γ_D = openarc(γ, γs, γss, N; dom=[π/2, 3π/2])
        Γ_N = openarc(γ, γs, γss, N; dom=[-π/2, π/2], quad=:open)
        n̂   = (imag.(Γ_N.γt) .- im .* real.(Γ_N.γt)) ./ abs.(Γ_N.γt)
        σ, ϕ∞ = laplace_solve([Γ_D, Γ_N], [:D, :N], [real.(f.(Γ_D.γ)), real.(fp.(Γ_N.γ) .* n̂)])
        u_num = [ϕ∞ + real(𝒱(ζ, σ[1], Γ_D) + 𝒱(ζ, σ[2], Γ_N)) for ζ in ζ_pts]
        err = norm(u_num .- real.(f.(ζ_pts)), Inf)
        printerror(N, err)
    end
end
