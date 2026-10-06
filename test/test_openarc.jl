# Dirichlet solve on an open elliptical arc with the `:open` quadrature (smooth density), via
# the product-integration single layer. The density is prescribed as a smooth function of
# the arc's reference variable t, the Dirichlet data are its single layer computed by
# adaptive quadrature, and the potential at an off-arc point is compared with the exact one.

using QuadGK

@testset "Open arc: smooth-density Dirichlet solve" begin
    section("Open arc: smooth-density Dirichlet solve, refining boundary resolution N")
    Nspan = sweep(16:8:48, [24])

    a, b = 2.0, 1.0
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)
    sa, sb = -π/2, π/2

    t_of(s) = 2(s - sa) / (sb - sa) - 1
    σ_exact(t) = exp(t) * cos(2t)

    # Single layer of the exact density at x, split at the parameter of x when it lies on
    # the arc
    function 𝒱_exact(x, split=nothing)
        integrand(s) = (r = abs(x - γ(s)); r == 0 ? 0.0 : -log(r) / (2π)) * σ_exact(t_of(s)) * abs(γs(s))
        pts = split === nothing ? (sa, sb) : (sa, split, sb)
        return quadgk(integrand, pts...; rtol=1e-13, atol=1e-15)[1]
    end

    ζ = 3.0 + 1.0im
    u_ζ = 𝒱_exact(ζ)

    errs = map(Nspan) do N
        Γ = openarc(γ, γs, γss, N; dom=[sa, sb], quad=:open)
        g = [𝒱_exact(Γ.γ[i], Γ.s[i]) for i in 1:N]
        σ = Γ.V \ g
        abs(𝒱(ζ, σ, Γ) - u_ζ)
    end
    check_errors(Nspan, errs)
end
