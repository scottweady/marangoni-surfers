
using LinearAlgebra
using Plots

include("boundary.jl")
include("operators.jl")

# Manufactured solution: u = r cos θ = Re(z), harmonic
# On unit circle: u = cos θ = Re(x),  ∂u/∂n = cos θ = Re(x)

function solve_mixed(N)
    Γ_D = arc(N; dom=[0.0, π])
    Γ_N = arc(N; dom=[π,  2π])
    g_D = real.(Γ_D.x)
    g_N = real.(Γ_N.x)
    α_D, α_N = solve_mixed(Γ_D, Γ_N, g_D, g_N)
    return Γ_D, Γ_N, α_D, α_N
end

function eval_potential(z, α_D, α_N, Γ_D, Γ_N)
    N_D, N_N = length(Γ_D), length(Γ_N)
    u = sum(j -> -(π/N_D)/(2π) * log(abs(z - Γ_D.x[j])) * α_D[j] * abs(Γ_D.xs[j]), 1:N_D)
    u += sum(j -> -(π/N_N)/(2π) * log(abs(z - Γ_N.x[j])) * α_N[j] * abs(Γ_N.xs[j]), 1:N_N)
    return u
end

z₀ = 0.5 + 0.0im

Ns   = [8, 16, 32, 64, 128]
errs = map(Ns) do N
    Γ_D, Γ_N, α_D, α_N = solve_mixed(N)
    abs(eval_potential(z₀, α_D, α_N, Γ_D, Γ_N) - real(z₀))
end

for (N, err) in zip(Ns, errs)
    println("N = $N   error = $err")
end

Γ_D, Γ_N, α_D, α_N = solve_mixed(128)
w_D = sqrt.(1 .- Γ_D.s.^2)
w_N = sqrt.(1 .- Γ_N.s.^2)
plt = plot(Γ_D.s, α_D ./ w_D)
plot!(Γ_N.s.+2, α_N ./ w_N)
display(plt)
readline()