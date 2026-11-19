# Direct numerical simulation of a tandem chain of surfers: `N` particles spaced `d` apart on
# the x axis, all pointing along e_x, with the orientations perturbed by a wave
# δθᵢ = ε cos(2π m i/N) of `m` wavelengths along the chain and by Gaussian noise of standard
# deviation `δθ`.
#
# Output: <outdir>/frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

using Random

using MarangoniSurfers

# -- Parameters ----------------------------------------------------------------

Ns   = 64
M    = 16
s₀   = π/2
g(s) = -cos(s)

aspect_ratio = script_param(:aspect_ratio, 4.0)
N            = script_param(:N, 20)
d            = script_param(:d, 10.0)       # initial spacing
m            = script_param(:m, 10)         # wavelengths along the chain
ε            = script_param(:ε, 0.01)       # wave amplitude
δθ           = script_param(:δθ, 0.0)       # orientation noise amplitude
seed         = script_param(:seed, 1)
t_end        = script_param(:t_end, 1000.0)
save_every   = script_param(:save_every, 10.0)
render       = script_param(:plot, true)    # render, display and save each frame
render && @eval using Plots   # loaded only when plotting
outdir       = script_param(:outdir, joinpath(@__DIR__, "output",
                                             "tandem_$(aspect_ratio)_N$(N)_d$(d)" *
                                             (ε > 0 ? "_m$(m)_ε$(ε)" : "") *
                                             (δθ > 0 ? "_δθ$(δθ)" : "")))
a = sqrt(aspect_ratio)
b = 1.0 / a

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Run -----------------------------------------------------------------------

xcm = ComplexF64.(d .* (0:N-1))
θcm = ε .* cos.(2π * m .* (0:N-1) ./ N) .+ δθ .* randn(Xoshiro(seed), N)

isdir(outdir) && rm(outdir; recursive=true)
suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g; t_end=t_end, dt_max=1.0, dt_min=1e-3,
                    plot_every=save_every, plot_dir=outdir, render=render)
