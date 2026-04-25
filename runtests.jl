
using LinearAlgebra

include("boundary.jl")
include("operators.jl")

function printerror(N, err)

  # If below 1e-12, green; if above 1e-6, red; otherwise, yellow
  color = if err < 1e-12
      :green
  elseif err > 1e-6
      :red
  else
      :yellow
  end

  printstyled("N = $N, error = $err\n", color=color)
  
end

Ns = [8, 16, 32, 64, 128]

# ── Segment test ──────────────────────────────────────────────────────────────
# Exact solution for -1/(2π) ∫₋₁¹ log|x-y| / √(1-y²) dy = 1, x ∈ [-1,1]
# Smooth density: α_exact = 2/log(2) (constant)
println("Segment: single-layer, constant RHS")
α_exact = 2 / log(2)
for N in Ns
    Γ   = segment(N)
    α   = Γ.V \ ones(N)
    err = norm(α .- α_exact) / abs(α_exact)
    printerror(N, err)
end

# ── Circle test 1 ─────────────────────────────────────────────────────────────
# K′[1] = 1/2 exactly on any closed curve (interior jump condition)
println("\nCircle: adjoint double-layer applied to 1 (exact = 0.5)")
for N in Ns
    Γ   = circle(N)
    err = norm(Γ.V′ * ones(N) .+ 0.5) / (N * 0.5)
    printerror(N, err)
end

# ── Circle test 2 ─────────────────────────────────────────────────────────────
# Manufactured density σ = cos(s), solve V⁻¹(V σ) and recover σ
println("\nCircle: single-layer solve convergence (σ = cos(s))")
for N in Ns
    Γ        = circle(N)
    σ_exact  = cos.(Γ.s)
    g        = Γ.V * σ_exact
    σ        = V⁻¹(g, Γ)
    err      = norm(σ .- σ_exact) / norm(σ_exact)
    printerror(N, err)
end