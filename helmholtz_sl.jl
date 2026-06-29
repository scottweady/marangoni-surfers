
include(joinpath(@__DIR__, "src/marangoni.jl"))
include(joinpath(@__DIR__, "local/normal_derivatives.jl"))
include(joinpath(@__DIR__, "local/evolve.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

k   = 3.0
ξ₁, γ₁ = 1.0, 0.0
ξ₂, γ₂ = 1.0, 0.0
δ   = 0.1
β₁₂ =  1.0
β₂₁ = -1.0
R₂  = 3.0
N   = 256

params = (; ξ₁, γ₁, ξ₂, γ₂, δ, β₁₂, β₂₁, k)

# ── Initial conditions ────────────────────────────────────────────────────────

a = 0.1 # perturbation amplitude

Γ₁₀ = curve(
    s -> (1 + a*cos(5s)) * exp(im*s),
    s -> exp(im*s) * (im*(1 + a*cos(5s)) - 5a*sin(5s)),
    s -> exp(im*s) * (-(1 + 26a*cos(5s)) - 10im*a*sin(5s)),
    N)
Γ₂₀ = curve(s -> R₂*exp(im*s), s -> im*R₂*exp(im*s), s -> -R₂*exp(im*s), N)

θ̂₁, L₁, X₀₁ = boundary_to_tangent_angle(Γ₁₀)
θ̂₂, L₂, X₀₂ = boundary_to_tangent_angle(Γ₂₀)

# ── Euler time-stepping ───────────────────────────────────────────────────────

dt    = 0.5
T     = 1000.0
nstep = round(Int, T/dt)

pl = plot(
    real.(Γ₁₀.γ), imag.(Γ₁₀.γ); lc=:gray, ls=:dash, label="t=0",
    aspect_ratio=:equal, legend=:topright
)
plot!(pl, real.(Γ₂₀.γ), imag.(Γ₂₀.γ); lc=:gray, ls=:dash, label=false)
display(pl)

for step in 1:nstep
    dθ̂₁, dL₁, dX₀₁, dθ̂₂, dL₂, dX₀₂ = evolve_rhs(θ̂₁, L₁, X₀₁, θ̂₂, L₂, X₀₂; params...)
    global θ̂₁  = θ̂₁  .+ dt .* dθ̂₁
    global L₁   = L₁   +  dt  * dL₁
    # global X₀₁  = X₀₁  +  dt  * dX₀₂
    global θ̂₂  = θ̂₂  .+ dt .* dθ̂₂
    global L₂   = L₂   +  dt  * dL₂
    # global X₀₂  = X₀₂  +  dt  * dX₀₂

    if step % 10 == 0
        t  = step * dt
        Γ₁ = tangent_angle_to_boundary(θ̂₁, L₁, X₀₁)
        Γ₂ = tangent_angle_to_boundary(θ̂₂, L₂, X₀₂)
        α  = t / T
        plot!(pl, real.(Γ₁.γ), imag.(Γ₁.γ); lc=RGB(0, 0, α), label=false)
        plot!(pl, real.(Γ₂.γ), imag.(Γ₂.γ); lc=RGB(α, 0, 0), label=false)
        display(pl)
        @printf("t = %.4f\n", t)
    end
end
readline()
