
# Compare ∇ℬ_fmm against direct ∇ℬ for increasing numbers of particles.
# Particles placed on a ring of radius R; each has N boundary nodes.
# All-to-all evaluation: for each particle i, compute ∇ℬ at its nodes
# from all other particles' sources.

println("\n∇ℬ_fmm vs ∇ℬ (direct): scaling with number of particles (N=64 nodes each)")

N = 256
Γ0 = circle(N)

function translate(Γ::Boundary, z0::ComplexF64)
    Boundary(Γ.open, Γ.t, Γ.s, Γ.γ .+ z0, Γ.γt, Γ.γtt, Γ.w, Γ.V, Γ.V′)
end

function run_scaling_test(Npart)
    z0 = [4.0*k + 0im for k in 0:Npart-1]
    Γs = [translate(Γ0, z) for z in z0]
    σs = [cos.(Γs[k].s) for k in 1:Npart]

    # Interior evaluation points: 0.5-radius nodes centered on each particle.
    # These are distinct from source nodes, mirroring marangoni_solve where
    # targets (Stokes boundary) never coincide with sources (arc nodes).
    ζs = [0.5 .* Γ0.γ .+ z0[i] for i in 1:Npart]

    # Direct: all-to-all cross-particle ∇ℬ
    t_dir = @elapsed begin
        for i in 1:Npart
            ζ = ζs[i]
            bx = zeros(N);  by = zeros(N)
            for j in 1:Npart
                j == i && continue
                bxj, byj = ∇ℬ(ζ, σs[j], Γs[j])
                bx .+= bxj;  by .+= byj
            end
        end
    end

    # FMM: single batched call (all sources → all targets), then subtract self
    offs  = [0; cumsum(fill(N, Npart))]
    ζ_all = vcat(ζs...)
    t_fmm = @elapsed begin
        bx_all, by_all = ∇ℬ_fmm(ζ_all, σs, Γs)
        for i in 1:Npart
            ri = offs[i]+1:offs[i+1]
            bxs, bys = ∇ℬ(ζ_all[ri], σs[i], Γs[i])
            bx_all[ri] .-= bxs;  by_all[ri] .-= bys
        end
    end

    # Error check: FMM cross vs direct cross for particle 1
    ζ = ζs[1]
    cross = [j for j in 1:Npart if j != 1]
    bx_dir = zeros(N);  by_dir = zeros(N)
    for j in cross
        bxj, byj = ∇ℬ(ζ, σs[j], Γs[j])
        bx_dir .+= bxj;  by_dir .+= byj
    end
    bx_fmm = bx_all[1:N]
    err = norm(bx_fmm .- bx_dir, Inf) / max(norm(bx_dir, Inf), 1e-15)

    @printf("  Npart = %4d  |  direct = %7.2f ms  |  FMM = %7.2f ms  |  rel err = %.2e\n",
            Npart, t_dir*1e3, t_fmm*1e3, err)
end

using Printf
for Npart in [2, 4, 8, 16, 32, 64, 128]
    run_scaling_test(Npart)
end
