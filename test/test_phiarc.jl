# Single layer on φ-parametrized elliptical arcs (`ellipse_arc`), for a physical density
# σ = q/sin φ with q smooth in ψ, as at a Dirichlet-Neumann junction. Compared with adaptive
# quadrature in φ, at the arc's own nodes (self) and at the nodes of the adjoining arc
# (cross, nearly singular near the junctions).

using QuadGK

@testset "φ-arc: single layer" begin
    section("φ-arc: single layer of a junction-type density, refining boundary resolution N")
    Nspan = sweep(32:24:128, [96])

    a, b, s₀ = 2.0, 0.5, π/2
    γ(s) = a * cos(s) + im * b * sin(s)
    speed(s) = abs(-a * sin(s) + im * b * cos(s))
    q(ψ) = exp(ψ) * cos(2ψ) + 0.5 + 0.3 * sin(3ψ)

    # ∫ log|x - γ| σ dΓ over the arc [sa, sb], in φ, split at φx when x is on the arc
    function 𝒱_exact(x, sa, sb; φx=nothing)
        J = (sb - sa) / 2
        integrand(φ) = (s = sa + J * (cos(φ) + 1); r = abs(x - γ(s));
                        (r == 0 ? 0.0 : -log(r) / (2π)) * q(1 - 2φ/π) * speed(s) * J)
        pts = φx === nothing ? (0.0, π) : (0.0, φx, π)
        return quadgk(integrand, pts...; rtol=1e-13, atol=1e-15)[1]
    end

    errs_self, errs_cross = Float64[], Float64[]
    for N in Nspan
        Γ_D = ellipse_arc(a, b, s₀, 2π - s₀, N)
        Γ_N = ellipse_arc(a, b, -s₀, s₀, N)
        σ_D = q.(Γ_D.t) ./ sin.(π .* (1 .- Γ_D.t) ./ 2)
        σ_N = q.(Γ_N.t) ./ sin.(π .* (1 .- Γ_N.t) ./ 2)
        φ_D = π .* (1 .- Γ_D.t) ./ 2

        u_self = Γ_D.V * σ_D
        push!(errs_self, maximum(abs(u_self[i] - 𝒱_exact(Γ_D.γ[i], s₀, 2π - s₀; φx=φ_D[i])) for i in 1:N))

        V_DN = _laplace2d_single_layer_phi_cross(a, b, -s₀, s₀, Γ_N.t, Γ_N.w, Γ_D.s)
        u_cross = V_DN * σ_N
        push!(errs_cross, maximum(abs(u_cross[i] - 𝒱_exact(Γ_D.γ[i], -s₀, s₀)) for i in 1:N))
    end

    subsection("self, at the arc's nodes")
    check_errors(Nspan, errs_self)

    subsection("cross, at the adjoining arc's nodes")
    check_errors(Nspan, errs_cross)
end
