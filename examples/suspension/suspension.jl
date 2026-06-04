
include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Ns = 64
M = 16
α = π/2
a = 1.5
b = 1.0 / a
g(θ)  = -cos(θ)

Ω   = ellipse(a, b, M)
γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α, α])

# ── Run ───────────────────────────────────────────────────────────────────────

Nparticles = 31
d = 5.0
xcm = [d*(i-1-(Nparticles-1)/2) .+ 0im for i in 1:Nparticles]   # center of mass position
θcm = (π/2) * ones(Nparticles)   # tangent orientation (counterclockwise)

isdir("frames") && rm("frames"; recursive=true)
println("Running simulation...")
traj_xcm, traj_θcm, traj_t = suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g; t_end=1.0e8, plot_mode=:concentration, plot_every=10.0, output_file="distances.dat", plot_dir="frames", dt_max=10.0)

readline()
