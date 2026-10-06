# Accuracy of the far-field rigid-body model for a pair of surfers: compares the
# translational and angular velocities from the full solve with the far-field formulas,
# parametrized by the self-propulsion speed U, surfactant dipole strength η, and
# Bretherton constant β, over `separations` for each orientation in `θcm_configs`.
#
# Output: output/farfield_dynamics.dat
#   d  err_u  err_ω  config

include("../../src/marangoni.jl")

# -- Geometry ------------------------------------------------------------------

Ns   = 128
M    = 32
s₀   = π/2
g(s) = -cos(s)

aspect_ratio = 4.0
a = sqrt(aspect_ratio)
b = 1.0 / a

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)
self_op = marangoni_self_operator(Ω, Γ_D, Γ_N)

# -- Single-particle constants -------------------------------------------------

U_sp, σ_sp, Γ_sp, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)
σ_D_sp, σ_N_sp = σ_sp
Γ_D_sp, Γ_N_sp = Γ_sp

U = U_sp[1]   # self-propulsion speed

# η: ∫ x σ dΓ = -η p, with p = e1 → η = -m1[1]
η = -(integral(real.(Γ_D_sp.γ) .* σ_D_sp, Γ_D_sp) +
      integral(real.(Γ_N_sp.γ) .* σ_N_sp, Γ_N_sp))

β = ((a/b)^2 - 1.0) / ((a/b)^2 + 1.0)

@printf("%-18sU = %.6f   η = %.6f   β = %.6f\n", "Single particle:", U, η, β)

# -- Far-field formulas --------------------------------------------------------

function ucm_ff(i, xcm, θcm)
    p_i = [cos(θcm[i]), sin(θcm[i])]
    u   = U .* p_i
    for j in eachindex(xcm)
        j == i && continue
        p_j   = [cos(θcm[j]), sin(θcm[j])]
        r_vec = [real(xcm[i] - xcm[j]), imag(xcm[i] - xcm[j])]
        rn    = norm(r_vec);  r̂ = r_vec / rn
        u   .+= (I(2) - r̂ * r̂') * (η * p_j) / (4π * rn)
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
        ω    += -β * sum(pp_i .* (η * T * p_i)) / (4π * rn^2)
    end
    return ω
end

# -- Sweep separations ---------------------------------------------------------

separations = [8.0, 16.0, 32.0, 64.0, 128.0, 256.0]
θcm_configs = [[-π/2, π/2]]

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
outfile = joinpath(outdir, "farfield_dynamics.dat")

open(outfile, "w") do io
    @printf(io, "# d   err_u   err_ω   config\n")

    for (icfg, θcm) in enumerate(θcm_configs)
        @printf("\nθcm = (%.2fπ, %.2fπ)\n", (θcm ./ π)...)
        @printf("%8s%14s%14s\n", "d", "|Δu|", "|Δω|")

        for d in separations
            xcm  = [-d/2 + 0im, d/2 + 0im]
            U_bim, _, _, _ = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; self_op=self_op)

            err_u = err_ω = 0.0
            for i in 1:2
                u_bim = U_bim[i][1:2]
                ω_bim = U_bim[i][3]
                err_u = max(err_u, norm(u_bim - ucm_ff(i, xcm, θcm)))
                err_ω = max(err_ω, abs(ω_bim - ωcm_ff(i, xcm, θcm)))
            end

            @printf("%8.1f%14.2e%14.2e\n", d, err_u, err_ω)
            @printf(io, "%.6f  %.6e  %.6e  %d\n", d, err_u, err_ω, icfg)
        end
    end
end

println("\nWrote ", relpath(outfile, @__DIR__))
