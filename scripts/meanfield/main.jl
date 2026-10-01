
include(joinpath(@__DIR__, "meanfield.jl"))

using Plots, ColorSchemes
using SpecialFunctions: ellipk, ellipe

# Dimensionless parameters
a, b = 2.0, 0.5
λ  = (a^2 - b^2) / (a^2 + b^2)

U = 1.0
β = 10.0

dT = 0.005
dR = 0.005
L  = 2π

# Initial condition, :patch or :vonmises
config = :vonmises
tplt = 0.05
ζ  = 5.0
θ0 = 0.0

# Grid parameters
N  = 64
Nθ = 32
G  = SpectralGrid(N, N, Nθ; Lx=L, Ly=L)
tf = 200.0

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

const frame_count = Ref(0)

# Visualization
function render(snap, G; savedir="")

    format = (aspect_ratio = :equal, xlims = (0, G.Lx), ylims = (0, G.Ly), xlabel = "x", ylabel = "y")
    cformat = (c = cgrad(:algae,rev=true), clims = (0.0, 2.0))
    pρ = heatmap(G.x, G.y, snap.ρ'; cformat..., format...)
    p = plot(pρ)
    display(p)

    if !isempty(savedir)
        mkpath(savedir)
        savefig(p, joinpath(savedir, "frame_$(lpad(frame_count[], 6, '0')).png"))
    end

    frame_count[] += 1

    return p

end

# Main
if config == :patch
    R0 = 1.0;
    ρ0 = L^2/(π * R0^2)
    ψ0 = ic_patch(G; ρ0=ρ0, R=R0,width=10.0/N)
elseif config == :vonmises
    ψ0 = ic_polar(G; ρ0=1.0, θ0=θ0, ζ=ζ, modes = random_modes(12; nmin=1, nmax=4, ερ=1e-2, εm=0.0, seed=1))
end

savedir = joinpath(outdir, "frames")
simulate(ψ0, G; tf=tf, U=U, λ=λ, β=β, dT=dT, dR=dR, tplt=tplt, callback=(s, g) -> render(s, g; savedir=savedir))
