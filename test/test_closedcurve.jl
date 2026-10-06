# Dirichlet solve on a closed star-shaped curve via `closedcurve`, with exact solution
# u(x, y) = x. Checks the interior potential 𝒱σ + ϕ∞ against Re(ζ), and the boundary
# normal derivative σ/2 + K′σ against the exact Re(n̂).

@testset "Closed curve: Dirichlet solve" begin
    section("Closed curve: Dirichlet solve, u = x, refining boundary resolution N")
    Nspan = sweep(32:32:160, [96])

    a = 0.3
    γ(s)   = (1 + a * cos(5s)) * exp(im*s)
    γs(s)  = im * (1 + a * cos(5s))*exp(im*s) - 5a * sin(5s) * exp(im*s)
    γss(s) = -(1 + a * cos(5s)) * exp(im*s) - 10a * im * sin(5s) * exp(im*s) - 25a * cos(5s) * exp(im*s)

    ζ_pts = [0.0 + 0.0im, 0.2 + 0.1im, -0.3 + 0.2im, 0.1 - 0.2im]
    u_exact = real.(ζ_pts)

    function solve(N)
        Γ = closedcurve(γ, γs, γss, N)
        σs, ϕ∞ = laplace_solve([Γ], [:D], [real.(Γ.γ)])
        return Γ, σs[1], ϕ∞
    end

    subsection("interior potential")
    errs = map(Nspan) do N
        Γ, σ, ϕ∞ = solve(N)
        u_num = [𝒱(ζ, σ, Γ) + ϕ∞ for ζ in ζ_pts]
        norm(u_num .- u_exact, Inf)
    end
    check_errors(Nspan, errs)

    # Outward normal n̂ = -i γt/|γt|, so ∂(Re ζ)/∂n = Im(γt)/|γt|
    subsection("boundary normal derivative")
    errs = map(Nspan) do N
        Γ, σ, _ = solve(N)
        dudn_num   = σ / 2 + Γ.K′ * σ
        dudn_exact = imag.(Γ.γt) ./ abs.(Γ.γt)
        norm(dudn_num .- dudn_exact, Inf)
    end
    check_errors(Nspan, errs)
end
