# Mean-field simulation (meanfield.jl) from a patch, perturbed von Mises, or band
# initial condition, selected by `config`.
#
# Output: output/frames/frame_<n>.png, density with surface velocity alongside the
# polarization fluctuation.

include(joinpath(@__DIR__, "meanfield.jl"))

using Plots

# Dimensionless parameters
aspect_ratio = 4.0
a = sqrt(aspect_ratio)
b = 1.0 / a
β = (a^2 - b^2) / (a^2 + b^2)

U = 1.0
η = 10.0

dT = 0.005
dR = 0.005
L  = 2π

# Initial condition, :patch, :vonmises, or :band
config = :vonmises
tplt = 0.1
ζ  = 10.0
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
    s  = 4
    is = 1:s:G.Nx; js = 1:s:G.Ny
    X  = [G.x[i] for i in is, j in js]; Y = [G.y[j] for i in is, j in js]
    sc = 2 * s * G.dx

    # Left: density with surface velocity u, longest arrow ≈ 2 arrow spacings
    pρ = heatmap(G.x, G.y, snap.ρ'; format..., title="ρ, u")
    U1 = snap.u1[is, js]; U2 = snap.u2[is, js]
    umax = max(maximum(hypot.(U1, U2)), eps())
    ku = vec(hypot.(U1, U2) .> 0.05 * umax)
    quiver!(pρ, vec(X)[ku], vec(Y)[ku], quiver=(sc .* vec(U1)[ku] ./ umax, sc .* vec(U2)[ku] ./ umax), color=:white, lw=0.75)

    # Right: polarization fluctuation p - ⟨p⟩ and its magnitude, where p = m/ρ and ⟨p⟩ = ∫m/∫ρ;
    # longest arrow ≈ 2 arrow spacings
    δp = hypot.(snap.m1 ./ snap.ρ .- sum(snap.m1) / sum(snap.ρ), snap.m2 ./ snap.ρ .- sum(snap.m2) / sum(snap.ρ))
    δp[snap.ρ .< 0.05 * maximum(snap.ρ)] .= NaN   # p is ill-defined where ρ ≈ 0
    pp = heatmap(G.x, G.y, δp'; format..., title="|p − ⟨p⟩|, p − ⟨p⟩")
    R  = snap.ρ[is, js]
    P1 = snap.m1[is, js] ./ R .- sum(snap.m1) / sum(snap.ρ)
    P2 = snap.m2[is, js] ./ R .- sum(snap.m2) / sum(snap.ρ)
    pmax = max(maximum(hypot.(P1, P2)), eps())
    kp = vec((R .> 0.05 * maximum(snap.ρ)) .& (hypot.(P1, P2) .> 0.05 * pmax))
    quiver!(pp, vec(X)[kp], vec(Y)[kp], quiver=(sc .* vec(P1)[kp] ./ pmax, sc .* vec(P2)[kp] ./ pmax), color=:white, lw=0.75)

    p = plot(pρ, pp; layout=(1, 2), size=(1000, 420))
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
    ψ0 = ic_polar(G; ρ0=1.0, θ0=θ0, ζ=ζ, modes = random_modes(12; nmin=1, nmax=4, nθmax=4, ε=1e-2, seed=1))
elseif config == :band
    W  = 1.0
    ρ0 = L/W   # mean density 1
    # centerline displacement: random phases on wavenumbers 1-4 in y
    rng = MersenneTwister(1)
    ε  = 0.05
    φ  = 2π .* rand(rng, 4)
    δ(y) = ε * sum(cos(2π * n * y / L + φ[n]) for n in 1:4) / 4
    ψ0 = ic_band(G; ρ0=ρ0, W=W, width=10.0/N, θ0=θ0, ζ=ζ, δ=δ)
end

savedir = joinpath(outdir, "frames")
simulate(ψ0, G; tf=tf, U=U, β=β, η=η, dT=dT, dR=dR, tplt=tplt, callback=(s, g) -> render(s, g; savedir=savedir))
