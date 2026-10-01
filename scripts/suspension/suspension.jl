
include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Ns = 64
M = 16
α = π/2

aspect_ratio = 0.25
a = sqrt(aspect_ratio)
b = 1.0 / a
g(θ)  = -cos(θ)

Ω   = ellipse(a, b, M)
γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = openarc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α, α])

# ── Run ───────────────────────────────────────────────────────────────────────

config = :hedgehog   # :line, :hedgehog, :courtroom, or :ring

d = 5.0

pastel_red    = RGB(1.0, 0.71, 0.76)
pastel_blue   = RGB(0.68, 0.85, 1.0)
pastel_purple = RGB(0.80, 0.70, 0.95)

colors = Any[]
cmap = :viridis

if config == :line
    Nparticles = 13
    xcm = [d*(i-1-(Nparticles-1)/2) .+ 0im for i in 1:Nparticles]
    θcm = (π/2) * ones(Nparticles)
    colors = fill(:white, Nparticles)

elseif config == :hedgehog
    n_rings = 3
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
    Nparticles = length(xcm)
    colors = fill(:white, Nparticles)

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
    colors = fill(:white, Nparticles)

end

frames_dir = joinpath(@__DIR__, "frames_$(round(a/b, digits=4))")
isdir(frames_dir) && rm(frames_dir; recursive=true)
println("Running simulation...")
traj_xcm, traj_θcm, traj_t = suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g; t_end=200.0, plot_mode=:concentration, plot_every=5.0, output_file="hedgehog.dat", plot_dir=frames_dir, dt_max=5.0, colors=colors, cmap=cmap)
