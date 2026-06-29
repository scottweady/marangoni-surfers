
# Test: interior Dirichlet via double-layer + hypersingular normal derivative.
#
# u = 𝒦_Γ[σ]  (double-layer representation, interior)
# u = Re(z) on Γ  →  (-½I + K) σ = Re(γ)
# ∂u/∂n = T[σ]  (T = -Dₛᵀ V Dₛ, hypersingular)
#
# Exact solution: u = Re(z)  →  ∂u/∂n = ∇x · n̂ = Im(γ') / |γ'|

include(joinpath(@__DIR__, "src/marangoni.jl"))

function run_test(label, γ_fun, γt_fun, γtt_fun; Ns=[16, 32, 64, 128])
    println("\n── $label ──────────────────────────────")
    @printf("  %6s   %-20s  %-20s\n", "N", "‖σ residual‖∞", "‖∂u/∂n error‖∞")
    for N in Ns
        t   = collect(range(0, 2π, length=N+1)[1:N])
        γ   = γ_fun.(t)
        γt  = γt_fun.(t)
        γtt = γtt_fun.(t)
        w   = fill(2π/N, N)

        f      = real.(γ)                        # Dirichlet BC: u = Re(z)
        dudn_exact = imag.(γt) ./ abs.(γt)       # ∂u/∂n = Im(γ') / |γ'|

        K_mat  = _laplace2d_double_layer(t, γ, γt, γtt, w)
        σ      = (-0.5*I(N) + K_mat) \ f

        T      = _laplace2d_hypersingular(t, γ, γt, γtt)
        dudn   = T * σ

        res    = maximum(abs.((-0.5*I(N) + K_mat) * σ .- f))
        err    = maximum(abs.(dudn .- dudn_exact))
        @printf("  %6d   %-20.2e  %-20.2e\n", N, res, err)
    end
end

# Unit circle
run_test("Unit circle",
    t -> exp(im*t),
    t -> im*exp(im*t),
    t -> -exp(im*t))

# Ellipse (a=2, b=1)
run_test("Ellipse a=2, b=1",
    t -> 2cos(t) + im*sin(t),
    t -> -2sin(t) + im*cos(t),
    t -> -2cos(t) - im*sin(t))

# Starfish: r(t) = 1 + 0.3 cos(5t)
# γ = r e^{it},  γ' = (r' + ir) e^{it},  γ'' = (r''-r + 2ir') e^{it}
run_test("Starfish (0.3 cos 5t)",
    t -> (1 + 0.3cos(5t)) * exp(im*t),
    t -> (-1.5sin(5t) + im*(1 + 0.3cos(5t))) * exp(im*t),
    t -> (-(1 + 7.8cos(5t)) - 3im*sin(5t)) * exp(im*t))

# Kite (Kress 1991): γ = (cos t + 0.65 cos 2t - 0.65) + 1.5i sin t
run_test("Kite",
    t -> (cos(t) + 0.65cos(2t) - 0.65) + 1.5im*sin(t),
    t -> (-sin(t) - 1.3sin(2t)) + 1.5im*cos(t),
    t -> (-cos(t) - 2.6cos(2t)) - 1.5im*sin(t))
