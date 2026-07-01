
# Test: interior Dirichlet via double-layer + hypersingular normal derivative.
#
# u = 𝒦_Γ[σ]  (double-layer representation, interior)
# u = Re(z) on Γ  →  (-½I + K) σ = Re(γ)
# ∂u/∂n = T[σ]  (T = -Dₛᵀ V Dₛ, hypersingular)
#
# Exact solution: u = Re(z)  →  ∂u/∂n = ∇x · n̂ = Im(γ') / |γ'|

curves = Dict(
    "Unit circle" => (
        t -> exp(im*t),
        t -> im*exp(im*t),
        t -> -exp(im*t)),
    "Ellipse a=2, b=1" => (
        t -> 2cos(t) + im*sin(t),
        t -> -2sin(t) + im*cos(t),
        t -> -2cos(t) - im*sin(t)),
    "Starfish (0.3 cos 5t)" => (
        t -> (1 + 0.3cos(5t)) * exp(im*t),
        t -> (-1.5sin(5t) + im*(1 + 0.3cos(5t))) * exp(im*t),
        t -> (-(1 + 7.8cos(5t)) - 3im*sin(5t)) * exp(im*t)),
    "Kite" => (
        t -> (cos(t) + 0.65cos(2t) - 0.65) + 1.5im*sin(t),
        t -> (-sin(t) - 1.3sin(2t)) + 1.5im*cos(t),
        t -> (-cos(t) - 2.6cos(2t)) - 1.5im*sin(t)),
)

for (label, (γ_fun, γt_fun, γtt_fun)) in curves
    println("\nHypersingular normal derivative: $label")
    for N in Nspan
        t   = collect(range(0, 2π, length=N+1)[1:N])
        γ   = γ_fun.(t)
        γt  = γt_fun.(t)
        γtt = γtt_fun.(t)
        w   = fill(2π/N, N)

        f          = real.(γ)                     # Dirichlet BC: u = Re(z)
        dudn_exact = imag.(γt) ./ abs.(γt)        # ∂u/∂n = Im(γ') / |γ'|

        K_mat = _laplace2d_double_layer(t, γ, γt, γtt, w)
        σ     = (-0.5*I(N) + K_mat) \ f

        T    = _laplace2d_hypersingular(t, γ, γt, γtt)
        dudn = T * σ

        err = norm(dudn .- dudn_exact, Inf)
        printerror(N, err)
    end
end
