# Far-field model, with the monopole closure, of a tandem chain of surfers: `N` particles
# spaced `d` apart on the x axis, all pointing along e_x, with the orientations perturbed by
# a wave δθᵢ = ε cos(2π m i/N) of `m` wavelengths along the chain and by Gaussian noise of
# standard deviation `δθ`.
#
# Output: <outdir>/frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

using Random

include("../../src/script_params.jl")
include("../../src/farfield.jl")

aspect_ratio = script_param(:aspect_ratio, 4.0)
N            = script_param(:N, 8)
d            = script_param(:d, 10.0)       # initial spacing
m            = script_param(:m, 0)          # wavelengths along the chain
ε            = script_param(:ε, 0.0)        # wave amplitude
δθ           = script_param(:δθ, 0.01)      # orientation noise amplitude
seed         = script_param(:seed, 1)
t_end        = script_param(:t_end, 4000.0)
save_every   = script_param(:save_every, 10.0)
dt_max       = script_param(:dt_max, 1.0)
outdir       = script_param(:outdir, joinpath(@__DIR__, "output",
                                             "tandem_$(aspect_ratio)_N$(N)_d$(d)" *
                                             (ε > 0 ? "_m$(m)_ε$(ε)" : "") *
                                             (δθ > 0 ? "_δθ$(δθ)" : "")))

xcm = ComplexF64.(d .* (0:N-1))
θcm = ε .* cos.(2π * m .* (0:N-1) ./ N) .+ δθ .* randn(Xoshiro(seed), N)

U, η, β = farfield_constants(aspect_ratio)
@printf("a/b = %g   N = %d   d = %.1f   U = %.6f   η = %.6f   β = %.4f\n",
        aspect_ratio, N, d, U, η, β)

isdir(outdir) && rm(outdir; recursive=true)
farfield_simulate(xcm, θcm, aspect_ratio; t_end=t_end, dt_max=dt_max, save_every=save_every,
                  closure=true, outdir=outdir, verbose=false)
