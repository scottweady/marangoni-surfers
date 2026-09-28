

include("../../src/marangoni.jl")

# ── Setup ─────────────────────────────────────────────────────────────────────

Ns    = 1024
M     = 32
α_arc = π/2
a  = 0.25
b = 1.0/a
g(s)  = -cos(s)

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Ω   = ellipse(a, b, M)
Γ_D = openarc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

U, σ, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)

ucm, ωcm = U[1:2], U[3]
σ_D, σ_N = σ
# ── Force density via the Stokes inverse ──────────────────────────────────────

x_Ω, y_Ω = real.(Ω.z), imag.(Ω.z)
bx_Ω = ∇ℬ(Ω.z, σ_D, Γ_D)[1] .+ ∇ℬ(Ω.z, σ_N, Γ_N)[1]
by_Ω = ∇ℬ(Ω.z, σ_D, Γ_D)[2] .+ ∇ℬ(Ω.z, σ_N, Γ_N)[2]
f = 𝒮_st⁻¹((ucm[1] .- y_Ω*ωcm .- bx_Ω, ucm[2] .+ x_Ω*ωcm .- by_Ω), Ω)

@printf("Force monopole = (%+.2e, %+.2e)\n",
        real(sum(f[1] .* Ω.dz)), real(sum(f[2] .* Ω.dz)))

function eval_stokes_sl(ζ, f, Ω)
    ux = uy = 0.0
    for i in eachindex(Ω.z)
        δz  = ζ - Ω.z[i]
        r   = abs(δz)
        rhx = real(δz)/r;  rhy = imag(δz)/r
        w   = real(Ω.dz[i])
        ux += ((1 + rhx^2)*f[1][i] + rhx*rhy*f[2][i]) * w / (4π*r)
        uy += (rhx*rhy*f[1][i] + (1 + rhy^2)*f[2][i]) * w / (4π*r)
    end
    return ux, uy
end

# ── Multipole moments ─────────────────────────────────────────────────────────

M0 = integral(σ_D, Γ_D) + integral(σ_N, Γ_N)
m1 = [integral(real.(Γ_D.γ) .* σ_D, Γ_D) + integral(real.(Γ_N.γ) .* σ_N, Γ_N),
      integral(imag.(Γ_D.γ) .* σ_D, Γ_D) + integral(imag.(Γ_N.γ) .* σ_N, Γ_N)]

@printf("Swimming speed = (%+.6f, %+.6f)\n", ucm...)      
@printf("Monopole M0 = %+.6f\n", M0)
@printf("Dipole   m1 = (%+.6f, %+.6f)\n", m1...)

# ── Compare full solution vs multipole at increasing r ────────────────────────

rvals = [8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
θvals = [0.0, π/4, π/2, π]

out_dir = joinpath(@__DIR__, "output")
mkpath(out_dir)

open(joinpath(out_dir, "multipole_1p.dat"), "w") do io
    @printf(io, "# r   err_th0   err_thpi4   err_thpi2   err_thpi\n")

    @printf("%-10s  %-12s  %-12s  %-12s  %-12s\n", "r", "err(θ=0)", "err(θ=π/4)", "err(θ=π/2)", "err(θ=π)")

    for r in rvals

        errs = Float64[]
        for θ in θvals
            r_hat = [cos(θ), sin(θ)]
            ζ = r * exp(im * θ)

            # Full numerical velocity: Stokes SL + ∇ℬ
            ux_st, uy_st = eval_stokes_sl(ζ, f, Ω)
            ux_D, uy_D   = ∇ℬ(ζ, σ_D, Γ_D)
            ux_N, uy_N   = ∇ℬ(ζ, σ_N, Γ_N)
            u_num = [ux_st + ux_D + ux_N, uy_st + uy_D + uy_N]

            # Multipole formula
            proj = r_hat[1]*m1[1] + r_hat[2]*m1[2]
            u_mp = (M0/(4π)) .* r_hat .- (1/(4π*r)) .* (m1 .- proj .* r_hat)

            push!(errs, hypot(u_num[1] - u_mp[1], u_num[2] - u_mp[2]))
        end

        @printf("%-10.2f  %-12.4e  %-12.4e  %-12.4e  %-12.4e\n", r, errs...)
        @printf(io, "%.6e  %.6e  %.6e  %.6e  %.6e\n", r, errs...)
    end
end

println("\nWrote output/multipole_1p.dat")
