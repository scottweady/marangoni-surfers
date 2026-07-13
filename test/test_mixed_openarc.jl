
using FFTW
using Statistics

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

    function solve_and_eval(N)
        Γ_D = openarc(γ, γs, γss, N; dom=[α, 2π-α])
        Γ_N = openarc(γ, γs, γss, N; dom=[-α, α])
        σ, ϕ∞ = laplace_solve([Γ_D, Γ_N], [:D, :N], [-cos.(3.21*Γ_D.s).*exp.(-Γ_D.s.^2), zeros(N)])
        val = ϕ∞ .+ 𝒱(ζ, σ[1], Γ_D) + 𝒱(ζ, σ[2], Γ_N)
                   
        return val, σ, Γ_D, Γ_N
    end

    Nref = 500
    uref, σref, Γ_D_ref, Γ_N_ref = solve_and_eval(Nref)

    # p = plot(Γ_D_ref.s, σref[1], marker=:circle, label="σ_D", xlabel="s", ylabel="σ")
    # plot!(Γ_N_ref.s, σref[2], marker=:circle, label="σ_N")
    # display(p)
    # readline()
    # σ̂_D = _cheb_coeffs(σref[1])
    # σ̂_N = _cheb_coeffs(σref[2])

    # # Fit tail slope on log-log scale: |σ̂_D[k]| ~ C * k^slope
    # tail_start = Nref ÷ 4
    # ks = tail_start:Nref
    # log_k = log10.(float.(ks))
    # log_σ = log10.(abs.(σ̂_D[ks]))
    # valid = isfinite.(log_σ)
    # slope_D, _ = [log_k[valid] ones(sum(valid))] \ log_σ[valid]
    # println("σ̂_D tail slope (log-log): ", round(slope_D; digits=3))
    # intercept_D = mean(log_σ[valid] .- slope_D .* log_k[valid])
    # fit_D = 10 .^ (slope_D .* log10.(float.(ks)) .+ intercept_D)

    # p = plot(1:2:Nref, abs.(σ̂_D[1:2:end]), marker=:circle, xscale=:identity, yscale=:log10, label="σ̂_D", xlabel="N", ylabel="|σ̂|")
    # plot!(1:2:Nref, abs.(σ̂_N[1:2:end]), marker=:circle, label="σ̂_N")
    # plot!(ks, fit_D, linestyle=:dash, label="slope $(round(slope_D; digits=2))")
    # ylims!(1e-16, 1e0)
    # yticks!(10 .^ (-16.0:4.0:0.0))
    # display(p)
    # readline()

    for N in Nspan
        u, _, _, _ = solve_and_eval(N)
        err = abs(u - uref)
        printerror(N, err)
    end
    
end
