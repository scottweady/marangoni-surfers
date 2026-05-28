
using ProjectedSphericalHarmonics
using Printf

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

# ── Setup ─────────────────────────────────────────────────────────────────────

Ns    = 256
M     = 32
α_arc = π/2
a, b  = 1.0, 1.0
g(s)  = -cos(s)

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Ω   = ellipse(a, b, M)
Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Single-particle multipoles (body frame, θ = 0) ────────────────────────────

_, _, α_D_sp, α_N_sp, Γ_D_sp, Γ_N_sp, _ = marangoni_mobility(Ω, Γ_D, Γ_N, g)

M0_sp = integral(α_D_sp, Γ_D_sp) + integral(α_N_sp, Γ_N_sp)
m1_sp = [integral(real.(Γ_D_sp.γ) .* α_D_sp, Γ_D_sp) + integral(real.(Γ_N_sp.γ) .* α_N_sp, Γ_N_sp),
         integral(imag.(Γ_D_sp.γ) .* α_D_sp, Γ_D_sp) + integral(imag.(Γ_N_sp.γ) .* α_N_sp, Γ_N_sp)]

@printf("Single-particle: M0 = %+.6f  m1 = (%+.6f, %+.6f)\n", M0_sp, m1_sp...)

# Rotate body-frame dipole to lab frame for particle at orientation θ
function rotate_dipole(m1, θ)
    c, s = cos(θ), sin(θ)
    return [c*m1[1] - s*m1[2], s*m1[1] + c*m1[2]]
end

# ── Two-particle solve ────────────────────────────────────────────────────────

d   = 20.0
xcm = [-d/2 + 0im, d/2 + 0im]
θcm = [0.0, 0.0]   # orientation of each particle; vary to test θ = 0 vs π/2

# ── Force density on each particle's Stokes domain ───────────────────────────

function particle_force_density(k, xcm, θcm, U_all, σ, Γ, Ω)
    pts = exp(im*θcm[k]) .* Ω.z .+ xcm[k]
    xp, yp = real.(pts), imag.(pts)
    Uk, ωk = U_all[k][1:2], U_all[k][3]

    bx = zeros(size(pts))
    by = zeros(size(pts))
    for j in eachindex(σ)
        bxj, byj = ∇ℬ(pts, σ[j], Γ[j])
        bx .+= bxj
        by .+= byj
    end

    slip_x = Uk[1] .- yp .* ωk .- bx
    slip_y = Uk[2] .+ xp .* ωk .- by
    f = 𝒮_st⁻¹((slip_x, slip_y), Ω)

    # Placed domain: shift quadrature nodes for off-boundary 𝒮_st evaluation
    Ω_placed = (z = pts, dz = Ω.dz)
    return f, Ω_placed
end


# ── Max error over circle of radius r ────────────────────────────────────────
# Columns: r  max_err_th0  max_err_thpi2

rvals  = 10.0 .^ range(2, 4, length=40)
θcm_configs = [[0.0, 0.0], [π/2, π/2]]
Nθ_circle   = 64   # angular samples on each circle

max_errors = zeros(length(rvals), 2)

for (icfg, θcm_cfg) in enumerate(θcm_configs)

    # Re-solve for this orientation
    U_cfg, σ_cfg, Γ_cfg, _ = marangoni_solve(xcm, θcm_cfg, Ω, Γ_D, Γ_N, g)
    f_cfg = Vector{Any}(undef, 2)
    Ω_cfg = Vector{Any}(undef, 2)
    for k in 1:2
        f_cfg[k], Ω_cfg[k] =
            particle_force_density(k, xcm, θcm_cfg, U_cfg, σ_cfg, Γ_cfg, Ω)
    end

    for (ir, r) in enumerate(rvals)
        θ_circle = range(0, 2π, length=Nθ_circle+1)[1:end-1]
        err_max  = 0.0

        for θ in θ_circle
            ζ = ComplexF64(r * exp(im * θ))

            ux_num = uy_num = 0.0
            for k in 1:2
                ux_st, uy_st = 𝒮_st(ζ, f_cfg[k], Ω_cfg[k])
                ux_D, uy_D   = ∇ℬ(ζ, σ_cfg[2k-1], Γ_cfg[2k-1])
                ux_N, uy_N   = ∇ℬ(ζ, σ_cfg[2k],   Γ_cfg[2k])
                ux_num += ux_st + ux_D + ux_N
                uy_num += uy_st + uy_D + uy_N
            end

            ux_mp = uy_mp = 0.0
            for k in 1:2
                m1_k   = rotate_dipole(m1_sp, θcm_cfg[k])
                δ      = ζ - xcm[k]
                rk     = abs(δ)
                rk_hat = [real(δ)/rk, imag(δ)/rk]
                proj   = rk_hat[1]*m1_k[1] + rk_hat[2]*m1_k[2]
                u_k    = (M0_sp/(4π)) .* rk_hat .- (1/(4π*rk)) .* (m1_k .- proj .* rk_hat)
                ux_mp += u_k[1]; uy_mp += u_k[2]
            end

            err_max = max(err_max, hypot(ux_num - ux_mp, uy_num - uy_mp))
        end

        max_errors[ir, icfg] = err_max
    end
end

all_errors = max_errors

open(joinpath(@__DIR__, "multipole_error_2p.dat"), "w") do io
    for ir in eachindex(rvals)
        @printf(io, "%.6f  %.10e  %.10e\n",
                rvals[ir], all_errors[ir, 1], all_errors[ir, 2])
    end
end

println("Wrote multipole_error_2p.dat")
