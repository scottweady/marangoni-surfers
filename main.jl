

using ProjectedSphericalHarmonics
using Plots

include("boundary.jl")
include("operators.jl")

M = 16
N = 16

D = disk(M)
ζ = D.ζ

Γ = circle(N)

ϕ = 1.0 .+ 0.5 .* cos.(Γ.s)
σ = V⁻¹(ϕ, Γ)
fσ = ∇ℬ(ζ, σ, Γ)

U, _ = stokes_mobility_solve(fσ, (0.0, 0.0), 0.0, D)
err = abs(U[1] + π/16)
println("U = $U, error = $err")

x, y = real.(ζ), imag.(ζ)
plt = scatter(vec(x), vec(y), zcolor=vec(fσ[1]), marker=:circle, ms=10, aspect_ratio=1)
display(plt)
readline()