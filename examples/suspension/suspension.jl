
include("../../src/marangoni.jl")

# ── Parameters ────────────────────────────────────────────────────────────────

Ns = 128
M = 16
α = π/2
a = sqrt(1.0)
b = 1.0 / a
g(θ)  = -cos(θ)

Ω   = ellipse(a, b, M)
γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)
Γ_D = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α, α])

# ── Run ───────────────────────────────────────────────────────────────────────

config = :hedgehog   # :line or :hedgehog

d = 5.0

if config == :line
    Nparticles = 13
    xcm = [d*(i-1-(Nparticles-1)/2) .+ 0im for i in 1:Nparticles]
    θcm = (π/2) * ones(Nparticles)

elseif config == :hedgehog
    n_rings = 2
    xs = ComplexF64[]
    θs  = Float64[]
    for k in 1:n_rings
        r   = k * d
        n_k = round(Int, 2π * r / d)
        for j in 1:n_k
            φ = 2π * (j-1) / n_k
            push!(xs, r * exp(im * φ))
            push!(θs, 2π*rand())                 # pointing radially outward
        end
    end
    xcm = xs
    θcm = θs
    Nparticles = length(xcm)
end

isdir("data") && rm("data"; recursive=true)
println("Running simulation...")
traj_xcm, traj_θcm, traj_t = suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g; t_end=500.0, plot_mode=:concentration, plot_every=10.0, output_file="distances.dat", plot_dir="data", dt_max=10.0)

readline()
