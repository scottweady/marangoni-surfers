
include("../../src/marangoni.jl")

# ── Geometry ──────────────────────────────────────────────────────────────────

Ns    = 1024
M     = 32
α_arc = π/2
a, b  = 2.0, 1.0
g(s)  = -cos(s)

γ_c(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Ω   = ellipse(a, b, M)
Γ_D = arc(γ_c, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = arc(γ_c, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Single-particle constants (body frame, θ = 0, p = e1) ────────────────────

U_sp, σ_sp, Γ_sp, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
σ_D_sp, σ_N_sp = σ_sp
Γ_D_sp, Γ_N_sp = Γ_sp

U0 = U_sp[1]   # self-propulsion speed

# β: ∫ x σ dΓ = -β p, with p = e1 → β = -m1[1]
β = -(integral(real.(Γ_D_sp.γ) .* σ_D_sp, Γ_D_sp) +
      integral(real.(Γ_N_sp.γ) .* σ_N_sp, Γ_N_sp))

γ_B = ((a/b)^2 - 1.0) / ((a/b)^2 + 1.0)

@printf("U0 = %.6f   β = %.6f   γ = %.6f\n", U0, β, γ_B)

# ── Far-field formulas ────────────────────────────────────────────────────────

function ucm_ff(i, xcm, θcm)
    p_i = [cos(θcm[i]), sin(θcm[i])]
    u   = U0 .* p_i
    for j in eachindex(xcm)
        j == i && continue
        p_j   = [cos(θcm[j]), sin(θcm[j])]
        r_vec = [real(xcm[i] - xcm[j]), imag(xcm[i] - xcm[j])]
        rn    = norm(r_vec);  r̂ = r_vec / rn
        u   .+= (I(2) - r̂ * r̂') * (β * p_j) / (4π * rn)
    end
    return u
end

function ωcm_ff(i, xcm, θcm)
    p_i  = [cos(θcm[i]), sin(θcm[i])]
    pp_i = [-p_i[2], p_i[1]]
    ω    = 0.0
    for j in eachindex(xcm)
        j == i && continue
        p_j   = [cos(θcm[j]), sin(θcm[j])]
        r_vec = [real(xcm[i] - xcm[j]), imag(xcm[i] - xcm[j])]
        rn    = norm(r_vec);  r̂ = r_vec / rn
        T     = (I(2) - 3*(r̂*r̂')) * sum(r̂ .* p_j) + r̂*p_j' + p_j*r̂'
        ω    += -γ_B * sum(pp_i .* (β * T * p_i)) / (4π * rn^2)
    end
    return ω
end

# ── Sweep separations ─────────────────────────────────────────────────────────

separations   = [8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
θcm_configs   = [[-π/2, π/2]]
config_labels = ["θ=0", "θ=π/2"]

open(joinpath(@__DIR__, "farfield_ucm.dat"), "w") do io
    @printf(io, "# d   err_u_th0   err_ω_th0   err_u_thpi2   err_ω_thpi2\n")

    for (icfg, θcm) in enumerate(θcm_configs)
        @printf("\n%s\n", config_labels[icfg])
        @printf("%-10s  %-12s  %-12s\n", "d", "|err u|", "|err ω|")

        for d in separations
            xcm  = [-d/2 + 0im, d/2 + 0im]
            U_bim, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)

            err_u = err_ω = 0.0
            for i in 1:2
                u_bim = U_bim[i][1:2]
                ω_bim = U_bim[i][3]
                err_u = max(err_u, norm(u_bim - ucm_ff(i, xcm, θcm)))
                err_ω = max(err_ω, abs(ω_bim - ωcm_ff(i, xcm, θcm)))
            end

            @printf("%-10.2f  %-12.4e  %-12.4e\n", d, err_u, err_ω)
            @printf(io, "%.6f  %.6e  %.6e  %d\n", d, err_u, err_ω, icfg)
        end
    end
end

println("\nWrote farfield_ucm.dat")
