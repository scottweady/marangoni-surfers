
include(joinpath(@__DIR__, "src/marangoni.jl"))
using Printf, Plots

# ── Parameters ────────────────────────────────────────────────────────────────

N     = 64
N_ext = 64
R_ext = 5.0
δ     = 0.1
ξ1    = 1.0;  ξ2 = 2.0
a11   =  1.0;  a12 = -3.0
a21   = 2.0;  a22 = 0.0
dt    = 0.001
T     = 5.0

# ── Rebuild a Boundary from new complex node positions ────────────────────────
# Uses spectral (FFT) differentiation with respect to the uniform parameter t.

function _rebuild_boundary(γ::Vector{ComplexF64}, Γ_ref::Boundary)
    N    = length(γ)
    t    = Γ_ref.t
    w    = fill(2π / N, N)
    m    = Float64.(fftfreq(N) .* N)   # frequencies: 0,1,...,N/2,-N/2+1,...,-1
    γ̂    = fft(γ)
    γt   = ifft(im .* m .* γ̂)
    γtt  = ifft(-m.^2 .* γ̂)
    dΓ   = abs.(γt) .* w
    V    = _laplace2d_single_layer_kress(t, γ, γt)
    K′   = _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    return Boundary(:kress, Γ_ref.θ, t, t, γ, γt, γtt, w, dΓ, V, K′)
end

_normals(Γ::Boundary) = (imag.(Γ.γt) .- im.*real.(Γ.γt)) ./ abs.(Γ.γt)

# ── Initial conditions ────────────────────────────────────────────────────────

ε = 0.1; mode = 5
Γ = curve(
    s -> (1 + ε*cos(mode*s)) * exp(im*s),
    s -> exp(im*s) * (im*(1 + ε*cos(mode*s)) - ε*mode*sin(mode*s)),
    s -> exp(im*s) * (-(1 + (mode^2+1)*ε*cos(mode*s)) - 2im*ε*mode*sin(mode*s)),
    N
)
Γ_ext = curve(s -> R_ext*exp(im*s), s -> im*R_ext*exp(im*s), s -> -R_ext*exp(im*s), N_ext)

params = (; ξ1, ξ2, a11, a12, a21, a22, δ)

# ── Time loop (RK2 / midpoint) ────────────────────────────────────────────────

nsteps = round(Int, T / dt)

pl = plot(real.(Γ.γ), imag.(Γ.γ);
    aspect_ratio=:equal, lc=:gray, ls=:dash, label="t=0", legend=:topright)
plot!(pl, real.(Γ_ext.γ), imag.(Γ_ext.γ); lc=:black, lw=1, label="Γ_ext")
display(pl)

for step in 1:nsteps

    # Stage 1
    U1, _ = coupled_pressure_solve(Γ, Γ_ext; params...)
    n̂     = _normals(Γ)
    γ_mid = Γ.γ .+ (dt/2) .* U1 .* n̂
    Γ_mid = _rebuild_boundary(γ_mid, Γ)

    # Stage 2
    U1_mid, _ = coupled_pressure_solve(Γ_mid, Γ_ext; params...)
    n̂_mid     = _normals(Γ_mid)
    γ_new     = Γ.γ .+ dt .* U1_mid .* n̂_mid
    global Γ  = _rebuild_boundary(γ_new, Γ)

    if step % 10 == 0
        t = step * dt
        α = t / T
        plot!(pl, real.(Γ.γ), imag.(Γ.γ); lc=RGB(0, 0, α), label=false)
        display(pl)
        @printf("t = %.3f  |Γ| = %.6f\n", t, sum(Γ.dΓ))
    end

end

readline()
