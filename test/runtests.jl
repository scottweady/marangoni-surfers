
using LinearAlgebra
using Plots
using ProjectedSphericalHarmonics

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

function printerror(N, err)

  # If below 1e-12, green; if above 1e-6, red; otherwise, yellow
  color = if err < 1e-12
      :green
  elseif err > 1e-4
      :red
  else
      :yellow
  end

  printstyled("N = $N, error = $err\n", color=color)
  
end

Nspan = [4, 8, 16, 32, 64]

# ── Segment test ──────────────────────────────────────────────────────────────
# Exact solution for -1/(2π) ∫₋₁¹ log|x-y| α(y) / √(1-y²) dy = λₙTₙ(y), x ∈ [-1,1]
# Smooth density: α_exact = T_n
L = 1.0
println("Segment: single-layer solve convergence (α = T₂(s), L = 2)")

# α_exact(s) = 1.0
# λ = log(2)/2

α_exact(s) = 2s^2 - 1
λ = L/4

for N in Nspan
    Γ = segment(N; dom=[-L, L])
    g = λ * α_exact.(Γ.s)
    α, _ = V⁻¹(g, Γ; decay=true)
    # display(α)
    err = norm(α .- α_exact.(Γ.s), Inf)
    printerror(N, err)
end

# ── Circle test 1 ─────────────────────────────────────────────────────────────
# K′[1] = 1/2 exactly on any closed curve (interior jump condition)
println("\nCircle: adjoint double-layer applied to 1 (exact = 0.5)")
for N in Nspan
    Γ   = circle(N)
    err = norm(Γ.V′ * ones(N) .+ 0.5, Inf)
    printerror(N, err)
end

# ── Circle test 2 ─────────────────────────────────────────────────────────────
# Manufactured density σ = cos(s), solve V⁻¹(V σ) and recover σ
println("\nCircle: single-layer solve convergence (σ = cos(s))")
for N in Nspan
    Γ = circle(N)
    σ_exact = 2 * cos.(Γ.s)
    g = cos.(Γ.s)
    σ, _ = V⁻¹(g, Γ)
    err = norm(σ .- σ_exact, Inf)
    printerror(N, err)
end

# ── Arc test ───────────────────────────────────────────────────────────────────
# Single-layer solve on the upper semicircle arc dom=[-π/2, π/2].
# RHS g(θ) = cos(2θ); convergence measured via the potential at an off-arc point.
println("\nArc (unit circle, dom=[-π/2,π/2]): single-layer solve, potential at z=2")
let
    γ(s)   = exp(1im * s)
    γs(s)  = 1im * exp(1im * s)
    γss(s) = -exp(1im * s)
    g(θ) = 0.2 * cos(2θ) + 0.5 * sin(3θ)

    Nref  = 512
    Γref  = arc(γ, γs, γss, Nref; dom=[-π/2, π/2])

    solve(Γref) = V(g.(angle.(Γref.γ)), Γref)#; decay=true)
    obj(Γ) = sum(Γ.s .* solve(Γ) .* Γ.w)
    aref = obj(Γref)

    for N in Nspan
        Γ   = arc(N; dom=[-π/2, π/2])
        a = obj(Γ)
        err = abs(aref - a)
        printerror(N, err)
    end
end

include("test_curve.jl")

# ── Circle surfer test ────────────────────────────────────────────────────────
# Exact swimming speed for a unit disk with g(θ) = cos(2θ), α = π/2.
println("\nCircle surfer: marangoni_mobility convergence in Ns")
let
    α = π/2
    g(θ) = cos(2θ)
    Uexact = (π/16) * (1 + cos(α)) * sin(α)^2
    Nq = 32
    for Ns in Nspan
        U, _ = marangoni_mobility(Ns, Nq, 1.0, 1.0, α, g)
        err = abs(U[1] - Uexact)
        printerror(Ns, err)
    end
end