
include(joinpath(@__DIR__, "src/marangoni.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

N     = 128
N_ext = 128
R_ext = 5.0       # radius of outer circle Γ_ext
δ     = 0.5
ξ1    = 1.0;  ξ2 = 2.0
a11   = 1.0;  a12 = 0.5
a21   = 0.5;  a22 = 1.0

Γ     = circle(N)
Γ_ext = curve(s -> R_ext*exp(im*s), s -> im*R_ext*exp(im*s), s -> -R_ext*exp(im*s), N_ext)

k = 1 / δ

# ── Step 0: solve H_Γ[σ] = 1 on Γ ───────────────────────────────────────────

H_self = _helmholtz2d_single_layer_kress(Γ.t, Γ.γ, Γ.γt; k=k)
σ      = H_self \ ones(N)
@show maximum(abs.(H_self * σ .- 1))

# ── Precompute outward normals ────────────────────────────────────────────────

n̂1    = (imag.(Γ.γt)     .- im.*real.(Γ.γt))     ./ abs.(Γ.γt)
n̂_ext = (imag.(Γ_ext.γt) .- im.*real.(Γ_ext.γt)) ./ abs.(Γ_ext.γt)

# ── Cross-boundary operator matrices ─────────────────────────────────────────

# Single-layer: sources on Γ_ext, targets on Γ
V_ext_on_Γ   = _laplace2d_single_layer(Γ_ext, Γ)

# Adjoint DL (normal at target): sources on Γ_ext, normals at Γ
K′_ext_on_Γ  = _laplace2d_double_layer_adjoint(Γ_ext, Γ)

# Adjoint DL: sources on Γ, normals at Γ_ext
K′_Γ_on_ext  = _laplace2d_double_layer_adjoint(Γ, Γ_ext)

# Helmholtz adjoint DL on Γ (∂/∂n of ℋ, limit from same boundary)
K′_H_self = _helmholtz2d_double_layer_adjoint(Γ.t, Γ.γ, Γ.γt, Γ.γtt, Γ.w; k=k)

# Helmholtz adjoint DL: sources on Γ, normals at Γ_ext (smooth cross-matrix)
K′_H_Γ_on_ext = zeros(N_ext, N)
for i in 1:N_ext, j in 1:N
    r   = Γ_ext.γ[i] - Γ.γ[j]
    ρ   = abs(r)
    K′_H_Γ_on_ext[i, j] = -(k/(2π)) * besselk(1, k*ρ) * real(conj(n̂_ext[i]) * r) / ρ * Γ.dΓ[j]
end

# ── Known forcing terms on each boundary ─────────────────────────────────────

# f1|_Γ = -ξ1*a12*δ²*H[σ] - ξ1*a11/4*|γ|²  (H[σ]=1 on Γ)
g1_Γ  = @. -ξ1*a12*δ^2 - ξ1*a11/4 * abs2(Γ.γ)

# f2|_Γ = -ξ2*a21*δ²*H[σ] - ξ2*a22/4*|γ|²
g2_Γ  = @. -ξ2*a21*δ^2 - ξ2*a22/4 * abs2(Γ.γ)

# ∂f1/∂n1|_Γ from inside Ω1 (- limit): uses ∂H/∂n|_- = (-½+K'_H)*σ
dn1_H_inside = (-0.5 .* σ) .+ K′_H_self * σ
dg1_dn1 = @. -ξ1*a12*δ^2 * dn1_H_inside - ξ1*a11/2 * real(conj(n̂1) * Γ.γ)

# ∂f2/∂n1|_Γ from outside Ω2 (+ limit): uses ∂H/∂n|_+ = (+½+K'_H)*σ
dn1_H_outside = (0.5 .* σ) .+ K′_H_self * σ
dg2_dn1 = @. -ξ2*a21*δ^2 * dn1_H_outside - ξ2*a22/2 * real(conj(n̂1) * Γ.γ)

# ∂f2/∂n2|_Γext (smooth cross): ∂H/∂n_ext = K'_H_Γ_on_ext * σ
dn2_H_ext  = K′_H_Γ_on_ext * σ
dg2_dn_ext = @. -ξ2*a21*δ^2 * dn2_H_ext - ξ2*a22/2 * real(conj(n̂_ext) * Γ_ext.γ)

# ── Block system: unknowns [mu1 (N), mu2 (N), tau (N_ext)] ───────────────────
#
# Ansatz:  p1 = V_Γ[mu1] + f1                     in Ω1
#          p2 = V_Γ[mu2] + V_Γext[tau] + f2        in Ω2
#
# (1) p1 = p2 on Γ:
#       V_Γ*mu1 - V_Γ*mu2 - V_ext_on_Γ*tau = g2_Γ - g1_Γ
#
# (2) (1/ξ1)∂p1/∂n1 = (1/ξ2)∂p2/∂n1 on Γ:
#       (1/ξ1)(-½+K′_Γ)*mu1 - (1/ξ2)(½+K′_Γ)*mu2 - (1/ξ2)K′_ext_on_Γ*tau
#         = (1/ξ2)dg2_dn1 - (1/ξ1)dg1_dn1
#
# (3) ∂p2/∂n2 = 0 on Γ_ext:
#       K′_Γ_on_ext*mu2 + (-½+K′_Γext)*tau = -dg2_dn_ext

K′_Γ  = Γ.K′
K′_ext = Γ_ext.K′

A = [
    Γ.V                                  -Γ.V                                    -V_ext_on_Γ;
    (1/ξ1).*(-0.5.*I(N) .+ K′_Γ)        -(1/ξ2).*(0.5.*I(N) .+ K′_Γ)          -(1/ξ2).*K′_ext_on_Γ;
    zeros(N_ext, N)                        K′_Γ_on_ext                             -0.5.*I(N_ext) .+ K′_ext
]
b = [
    g2_Γ .- g1_Γ;
    (1/ξ2).*dg2_dn1 .- (1/ξ1).*dg1_dn1;
    -dg2_dn_ext
]

sol  = A \ b
mu1  = sol[1:N]
mu2  = sol[N+1:2N]
tau  = sol[2N+1:end]

# ── Evaluate and plot p1 inside Γ and p2 between Γ and Γ_ext ─────────────────

using Plots

Np = 300
xp = range(-R_ext, R_ext, length=Np)
ζ  = complex.(xp', xp)

r_inner = abs.(ζ)
in1 = r_inner .< 1.0
in2 = (r_inner .>= 1.0) .& (r_inner .< R_ext)

H_grid  = ℋ(ζ, σ, Γ; k=k)
f1_grid = @. -ξ1*a12*δ^2 * H_grid - ξ1*a11/4 * abs2(ζ)
f2_grid = @. -ξ2*a21*δ^2 * H_grid - ξ2*a22/4 * abs2(ζ)

p1 = 𝒱(ζ, mu1, Γ) .+ f1_grid
p2 = 𝒱(ζ, mu2, Γ) .+ 𝒱(ζ, tau, Γ_ext) .+ f2_grid

p_plot = fill(NaN, Np, Np)
p_plot[in1] .= p1[in1]
p_plot[in2] .= p2[in2]

fig = heatmap(xp, xp, p_plot;
    aspect_ratio=:equal, colorbar=true, colorbar_title="p",
    xlabel="x", ylabel="y",
)
plot!(fig, real.(Γ.γ),     imag.(Γ.γ);     lc=:black, lw=1.5, label=false)
plot!(fig, real.(Γ_ext.γ), imag.(Γ_ext.γ); lc=:black, lw=1.5, label=false)
display(fig)

# ── Verify BCs ────────────────────────────────────────────────────────────────

p1_Γ = Γ.V * mu1 .+ g1_Γ
p2_Γ = Γ.V * mu2 .+ V_ext_on_Γ * tau .+ g2_Γ
@show maximum(abs.(p1_Γ .- p2_Γ))   # should be ~0

dp1_dn_Γ = (-0.5 .* mu1 .+ K′_Γ * mu1) .+ dg1_dn1
dp2_dn_Γ = ( 0.5 .* mu2 .+ K′_Γ * mu2) .+ K′_ext_on_Γ * tau .+ dg2_dn1
@show maximum(abs.(dp1_dn_Γ./ξ1 .- dp2_dn_Γ./ξ2))   # should be ~0

dp2_dn_ext = K′_Γ_on_ext * mu2 .+ (-0.5 .* tau .+ K′_ext * tau) .+ dg2_dn_ext
@show maximum(abs.(dp2_dn_ext))   # should be ~0

readline()
