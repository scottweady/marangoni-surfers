
# Single-layer solve on an open arc of a star curve γ=(1+0.3cos5s)eⁱˢ, dom=[-π/2, π/2].
# Solve V σ = g, then evaluate the potential at an off-arc point against a high-N reference.
println("\nArc (star, dom=[-π/2,π/2]): single-layer solve, potential at ζ=3")
let
    a = 0.3
    γ(s)   = (1 + a * cos(5s)) * exp(im*s)
    γs(s)  = im * (1 + a * cos(5s)) * exp(im*s) - 5a * sin(5s) * exp(im*s)
    γss(s) = -(1 + a * cos(5s)) * exp(im*s) - 10a * im * sin(5s) * exp(im*s) - 25a * cos(5s) * exp(im*s)
    g(s)   = 0.2 * cos(2s) + 0.5 * sin(3s)
    ζ      = 3.0 + 0.0im

    function solve_and_eval(N)
        Γ = arc(γ, γs, γss, N; dom=[-π/2, π/2])
        σ, uinf = V⁻¹(g.(Γ.s), Γ)
        return 𝒱(ζ, σ, Γ) + uinf
    end

    val_ref = solve_and_eval(512)
    for N in Nspan
        err = abs(solve_and_eval(N) - val_ref)
        printerror(N, err)
    end
end
