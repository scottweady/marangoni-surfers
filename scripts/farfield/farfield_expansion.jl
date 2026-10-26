# Far-field model of a large hedgehog patch (concentric rings of surfers pointing outward,
# as in scripts/suspension/suspension.jl) with the dipole-only model or, if `closure` is set,
# the monopole closure, with the outward orientations optionally perturbed by Gaussian noise
# of standard deviation `δθ`. Stops at
# `t_end` or on contact, unless `repulsion` is set, in which case particles are points with
# a soft repulsion of range `r_rep` and run through contact.
#
# Output: <outdir>/frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

using Random

include("../../src/script_params.jl")
include("../../src/farfield.jl")

aspect_ratio = script_param(:aspect_ratio, 1.0)
N            = script_param(:N, 10000)      # rings are added until the patch holds at least N
d            = script_param(:d, 8.0)        # ring spacing and in-ring spacing
t_end        = script_param(:t_end, 200.0)
save_every   = script_param(:save_every, 5.0)
dt_max       = script_param(:dt_max, 0.5)
δθ           = script_param(:δθ, 0.0)
seed         = script_param(:seed, 1)
repulsion    = script_param(:repulsion, 0.0)  # point particles with soft repulsion when positive
U_factor     = script_param(:U_factor, 1.0)   # scales the self-propulsion speed only
r_rep        = script_param(:r_rep, 2.0)
closure      = script_param(:closure, false)
outdir       = script_param(:outdir, joinpath(@__DIR__, "output",
                                             "expansion_$(aspect_ratio)_d$(d)" *
                                             (closure ? "_closure" : "") *
                                             (δθ > 0 ? "_δθ$(δθ)" : "") *
                                             (repulsion > 0 ? "_rep$(repulsion)" : "") *
                                             (U_factor != 1 ? "_U$(U_factor)" : "")))

xcm, θcm = ComplexF64[], Float64[]
k = 0
while length(xcm) < N
    global k += 1
    r   = k * d
    n_k = round(Int, 2π * r / d)
    for j in 1:n_k
        φ = 2π * (j - 1) / n_k
        push!(xcm, r * cis(φ))
        push!(θcm, φ)
    end
end

θcm .+= δθ .* randn(Xoshiro(seed), length(θcm))

U, η, β = farfield_constants(aspect_ratio)
@printf("a/b = %g   N = %d   rings = %d   U = %.6f   η = %.6f   β = %.6f\n",
        aspect_ratio, length(xcm), k, U_factor * U, η, β)

isdir(outdir) && rm(outdir; recursive=true)
farfield_simulate(xcm, θcm, aspect_ratio; t_end=t_end, dt_max=dt_max, save_every=save_every,
                  repulsion=repulsion, r_rep=r_rep, U_factor=U_factor, closure=closure,
                  outdir=outdir)
