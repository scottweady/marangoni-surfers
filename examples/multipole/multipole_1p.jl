

include("../../src/marangoni.jl")

# ── Setup ─────────────────────────────────────────────────────────────────────

Ns    = 1024
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

@printf("Monopole M0 = %+.6f\n", M0)
@printf("Dipole   m1 = (%+.6f, %+.6f)\n", m1...)

# ── Compare full solution vs multipole at increasing r ────────────────────────

rvals  = 10.0 .^ range(1, 3, length=40)
θvals  = [0.0, π/4, π/2, π]
colors = [:blue, :red, :green, :orange]

p = plot(xlabel="r", ylabel="|u_num - u_mp|", xscale=:log10, yscale=:log10,
         title="Multipole expansion error (full solution)")

e0 = Ref(NaN)

for (θ, col) in zip(θvals, colors)
    errors = Float64[]
    r_hat  = [cos(θ), sin(θ)]

    for r in rvals
        ζ = r * exp(im * θ)

        # Full numerical velocity: Stokes SL + ∇ℬ
        ux_st, uy_st = eval_stokes_sl(ζ, f, Ω)
        ux_D, uy_D   = ∇ℬ(ζ, σ_D, Γ_D)
        ux_N, uy_N   = ∇ℬ(ζ, σ_N, Γ_N)
        u_num = [ux_st + ux_D + ux_N, uy_st + uy_D + uy_N]

        # Multipole formula
        proj  = r_hat[1]*m1[1] + r_hat[2]*m1[2]
        u_mp  = (M0/(4π)) .* r_hat .- (1/(4π*r)) .* (m1 .- proj .* r_hat)

        push!(errors, hypot(u_num[1] - u_mp[1], u_num[2] - u_mp[2]))
    end

    isnan(e0[]) && (e0[] = errors[1])
    plot!(p, rvals, errors; color=col, lw=2, label="θ = $(round(θ/π, digits=2))π")
end

plot!(p, rvals, e0[] .* (rvals[1] ./ rvals).^2; ls=:dash, color=:black, label="O(r⁻²)")

display(p)
readline()
