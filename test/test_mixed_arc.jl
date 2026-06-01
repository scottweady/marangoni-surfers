
using Plots
using FFTW

# Mixed Dirichlet/Neumann solve on a split ellipse arc.
# g_D = cos(2s) on left half, g_N = 0 on right half [-α, α].
# Convergence measured via the potential at an off-curve point against a high-N reference.
println("\nMixed arc (ellipse a=2,b=1): g_D=cos(2s), g_N=0, potential at ζ=3")
let
    α = π/2
    a, b = 2.0, 1.0
    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    ζ = 3.0 + 0.0im

    function solve_and_eval(N; close_eval=true)
        Γ_D = arc(γ, γs, γss, N; dom=[α, 2π-α])
        Γ_N = arc(γ, γs, γss, N; dom=[-α, α])
        σ, ϕ∞ = laplace_solve([Γ_D, Γ_N], [:D, :N], [cos.(2 .* Γ_D.s) .+ exp.(Γ_D.s), zeros(N)]; close_eval=close_eval)
        return 𝒱(ζ, σ[1], Γ_D) + 𝒱(ζ, σ[2], Γ_N) + ϕ∞, σ
    end

    uref, _ = solve_and_eval(1024; close_eval=false)

    for N in Nspan
        u, _ = solve_and_eval(N; close_eval=true)
        err = abs(u - uref)
        printerror(N, err)
    end
end
