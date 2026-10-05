# Mixed exterior problem on an ellipse split into two arcs, with the exact solution
# u = Re Σ cₖ/(z - zₖ) from poles inside the ellipse, so u → 0 and ϕ∞ = 0. Dirichlet data
# u on Γ_D and exterior Neumann data ∂u/∂n on Γ_N. A wrong-sided jump term shows up as an
# O(1) error that does not converge.

println("\nMixed arcs: exterior manufactured solution, potential at far points")
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
