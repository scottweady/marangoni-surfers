# Far-field model of a uniformly polarized patch of surfers (square lattice in a disk, all
# pointing along e_x) with the orientations perturbed by a single plane wave
# δθ = ε cos(k·x), k = k (cos φ, sin φ), and by Gaussian noise of standard deviation `δθ`, to
# compare with the linear stability of the aligned mean-field state. Particles are points
# with a soft repulsion; with `closure` set, the monopole closure is used.
#
# Output: <outdir>/frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

using Random

include("../../src/script_params.jl")
include("../../src/farfield.jl")

aspect_ratio = script_param(:aspect_ratio, 4.0)
N            = script_param(:N, 10000)       # approximate number of particles
d            = script_param(:d, 24.0)        # lattice spacing
nλ           = script_param(:nλ, 4.0)        # perturbation wavelength in lattice spacings
φ            = script_param(:φ, 0.82)        # wavevector angle (fastest growing for a/b = 4)
ε            = script_param(:ε, 0.02)        # plane-wave amplitude
δθ           = script_param(:δθ, 0.0)        # noise amplitude
seed         = script_param(:seed, 1)
closure      = script_param(:closure, false)
t_end        = script_param(:t_end, 3000.0)
save_every   = script_param(:save_every, 25.0)
dt_max       = script_param(:dt_max, 5.0)
outdir       = script_param(:outdir, joinpath(@__DIR__, "output",
                                             "polarized_$(aspect_ratio)_N$(N)_d$(d)" *
                                             (ε > 0 ? "_nλ$(nλ)" : "") *
                                             (δθ > 0 ? "_δθ$(δθ)" : "") *
                                             (closure ? "_closure" : "")))

R = d * sqrt(N / π)
n = ceil(Int, R / d)
xcm = ComplexF64[d * (i + im * j) for i in -n:n, j in -n:n if d * hypot(i, j) <= R]
k = 2π / (nλ * d) * cis(φ)
θcm = ε .* cos.(real.(conj(k) .* xcm)) .+ δθ .* randn(Xoshiro(seed), length(xcm))

U, η, β = farfield_constants(aspect_ratio)
η̄ = η / d^2
@printf("a/b = %g   N = %d   R = %.1f   U = %.6f   η̄ = %.3e   k = %.4f   k′ = %.3f\n",
        aspect_ratio, length(xcm), R, U, η̄, abs(k), abs(k) * U / η̄)

isdir(outdir) && rm(outdir; recursive=true)
farfield_simulate(xcm, θcm, aspect_ratio; t_end=t_end, dt_max=dt_max, save_every=save_every,
                  repulsion=1.0, r_rep=2.0, closure=closure, outdir=outdir)
