
using Plots
using ProjectedSphericalHarmonics

include("boundary.jl")
include("operators.jl")

function marangoni_velocity(α, g_D, M, N)

  D = disk(M)
  ζ = D.ζ

  Γ_D = arc(N; dom=[α, 2π-α])
  Γ_N = arc(N; dom=[-α, α])

  g_D = g_D.(Γ_D.θ)
  g_N = zeros(length(Γ_N))
  α_D, α_N = solve_mixed(Γ_D, Γ_N, g_D, g_N)

  fx_D, fy_D = ∇ℬ_arc(ζ, α_D, Γ_D)
  fx_N, fy_N = ∇ℬ_arc(ζ, α_N, Γ_N)
  fx = fx_D .+ fx_N
  fy = fy_D .+ fy_N

  U, _ = stokes_mobility_solve((fx, fy), (0.0, 0.0), 0.0, D)
  return U, α_D, α_N, Γ_D, Γ_N

end

α = π/4
g(θ) = 1.0 + cos(2*θ)
U, α_D, α_N, Γ_D, Γ_N = marangoni_velocity(α, g, 32, 32)
c = -(1/4) * (1 + cos(α)) * sin(α)^2

Uexact = -(π/4) * c

println("U = $U")
println("Uexact = $Uexact")
println("Error = $(abs(U[1] - Uexact))")

# plt = plot(Γ_D.s, α_D)
# plot!(Γ_N.s .+ 2.0, α_N)
# display(plt)
# readline()
