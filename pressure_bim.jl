
include(joinpath(@__DIR__, "src/marangoni.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

N   = 128
δ   = 0.1
a11 = 1.0
a12 = 0.0

Γ = circle(N)

# ── Step 1: solve H_Γ[σ] = 1 on Γ ───────────────────────────────────────────
# ℋ has kernel K₀(k|x-y|)/(2π) with k = 1/δ.

k = 1 / δ
H = _helmholtz2d_single_layer_kress(Γ.t, Γ.γ, Γ.γt; k=k)
σ = H \ ones(N)

# Check: max residual of H_Γ[σ] - 1 on Γ
@show maximum(abs.(H * σ .- 1))

# ── Step 2: solve p = 0 on Γ for ν ───────────────────────────────────────────
# p = V_Γ[ν] - a12*δ²*H_Γ[σ] - a11*|x|²/4 = 0 on Γ
# Since H_Γ[σ] = 1 on Γ, rhs = a12*δ² + a11*|γ|²/4.
# V_Γ[ν] = rhs  →  Γ.V * ν = rhs.

rhs = a12 * δ^2 .+ a11 .* abs2.(Γ.γ) ./ 4

# Augmented system: enforce ∫_Γ ν ds = 0 to pin the null space.
A_aug = [Γ.V       ones(N);
         Γ.dΓ'     0.0    ]
b_aug = [rhs; 0.0]
ν_aug = A_aug \ b_aug
ν     = ν_aug[1:N]
p0 = ν_aug[N+1]

# ── Evaluate p in the interior ────────────────────────────────────────────────

x = range(-0.9, 0.9, length=200)
ζ = complex.(x', x)           # 200×200 grid, rows=y cols=x

p = 𝒱(ζ, ν, Γ) .- a12 .* δ^2 .* ℋ(ζ, σ, Γ; k=k) .- a11 .* abs2.(ζ) ./ 4 .+ p0

# Mask exterior points
inside = abs.(ζ) .< 0.9
p_plot = ifelse.(inside, p, NaN)

# ── Plot ──────────────────────────────────────────────────────────────────────

using Plots

fig = heatmap(x, x, p_plot;
    aspect_ratio = :equal,
    colorbar      = true,
    colorbar_title = "p",
    xlabel = "x", ylabel = "y"
)
plot!(fig, real.(Γ.γ), imag.(Γ.γ); lc=:black, lw=1.5, label=false)
display(fig)
readline()