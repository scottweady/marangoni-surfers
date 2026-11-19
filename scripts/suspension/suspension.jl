# Time-dependent simulation of a suspension of surfers from a line, hedgehog (concentric
# rings pointing outward), or perturbed ring configuration, selected by `config`.
#
# Output: output/frames_<a/b>/frame_<n>.png and frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

using MarangoniSurfers

# -- Parameters ----------------------------------------------------------------

Ns   = 64
M    = 16
s₀   = π/2
g(s) = -cos(s)

aspect_ratio = script_param(:aspect_ratio, 4.0)
t_end        = script_param(:t_end, 200.0)
render       = script_param(:plot, true)   # render, display and save each frame
a = sqrt(aspect_ratio)
b = 1.0 / a

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Run -----------------------------------------------------------------------

config = script_param(:config, :hedgehog)   # :line, :hedgehog, or :ring

d = 8.0

if config == :line
    Nparticles = 13
    xcm = [d*(i-1-(Nparticles-1)/2) .+ 0im for i in 1:Nparticles]
    θcm = (π/2) * ones(Nparticles)

elseif config == :hedgehog
    n_rings = 5
    xs = ComplexF64[]
    θs  = Float64[]
    for k in 1:n_rings
        r   = k * d
        n_k = round(Int, 2π * r / d)
        for j in 1:n_k
            φ = 2π * (j-1) / n_k
            push!(xs, r * exp(im * φ))
            push!(θs, φ)                 # pointing radially outward
        end
    end
    xcm = xs
    θcm = θs

elseif config == :ring
    Nparticles = 13
    R          = 8.0
    pert       = 0.1     # perturbation amplitude
    k_mode     = 5       # azimuthal wavenumber of the perturbation

    φ  = [2π*(j-1)/Nparticles for j in 1:Nparticles]
    r  = R .* (1 .+ pert .* cos.(k_mode .* φ))
    rp = -R .* pert .* k_mode .* sin.(k_mode .* φ)   # dr/dφ

    xcm = r .* exp.(im .* φ)

    # outward normal to the curve x(φ) = r(φ)*[cos(φ), sin(φ)]
    Nx  = rp .* sin.(φ) .+ r .* cos.(φ)
    Ny  = -rp .* cos.(φ) .+ r .* sin.(φ)
    θcm = atan.(Ny, Nx)

end

outdir = script_param(:outdir, joinpath(@__DIR__, "output", "frames_$(round(a/b, digits=4))"))
isdir(outdir) && rm(outdir; recursive=true)
traj_xcm, traj_θcm, traj_t = suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g;
                                                 t_end=t_end, dt_max=1.0, dt_min=1e-3,plot_every=5.0,
                                                 plot_dir=outdir, render=render)
