
include("../src/bim/boundary.jl")
include("../src/bim/laplace.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

using FFTW
N = 8192
Nclose = 512

function V_close(x, σ, Γ; t0=nothing)
    N = length(σ)
    val = 0.0
    for j = 1:N
        σ_j = zeros(N)
        σ_j[j] = σ[j]
        val += 𝒱_close(x, σ_j, Γ; t0=t0)
    end
    return val
end

# --- Test on a segment ---
println("=== Segment ===")
Γ = segment(N)
Γclose = segment(Nclose)
σ_fun(s) = exp(cos(s)) + sin(5s) .+ 0.1
σ = σ_fun.(Γ.s)
σclose = σ_fun.(Γclose.s)
z = 0.3 + 1e-3im
val = 𝒱(z, σ, Γ)
val_close = 𝒱_close(z, σclose, Γclose)
println("  direct:     ", val)
println("  close-eval: ", val_close)
println("  error:      ", abs(val - val_close))

# --- Test on an arc ---
println("\n=== Arc ===")
a = 0.3
γ(s)   = (1 + a * cos(5s)) * exp(im*s)
γs(s)  = im * (1 + a * cos(5s))*exp(im*s) - 5a * sin(5s) * exp(im*s)
γss(s) = -(1 + a * cos(5s)) * exp(im*s) - 10a * im * sin(5s) * exp(im*s) - 25a * cos(5s) * exp(im*s)

dom = [-π/2, π/2]
Γ = arc(γ, γs, γss, N; dom=dom)
σ = σ_fun.(Γ.s)

Γclose = arc(γ, γs, γss, Nclose; dom=dom)
σclose = σ_fun.(Γclose.s)

# target close to the arc near s=0 endpoint
z = γ(0.0) + 1e-3

val = 𝒱(z, σ, Γ)
val_close = 𝒱_close(z, σclose, Γclose)

println("  direct:     ", val)
println("  close-eval: ", val_close)
println("  error:      ", abs(val - val_close))