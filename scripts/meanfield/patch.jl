# Mean-field model (meanfield.jl) of an expanding patch, the continuum counterpart of the
# particle simulations in scripts/suspension/suspension.jl. Lengths and times are in the
# units of the particle simulations: U, η and β are the single-particle constants for the
# given aspect ratio, and the annulus Rin < r < R0 carries the number density and extent of
# the particle patch, whose rings sit at r = 5, 10, 15. Orientations are von Mises distributed about the outward radial direction. With
# `perturb=:orientation` the local mean orientation is rotated by a sum of random plane
# waves in (x, y) of amplitude ε radians, leaving the density unperturbed; with
# `perturb=:azimuthal` it is rotated by the single azimuthal mode ε cos(m φ), with φ the
# polar angle about the patch center, again leaving the density unperturbed; with
# `perturb=:modes` ψ is multiplied by 1 + Σ ε cos(k·x + nθ θ + φ), which for a sharply
# peaked distribution perturbs mostly the density.
#
# Output: <outdir>/frame_<n>.dat (see `simulate`) and <outdir>/patch.dat, one row per frame
#   t  R/R0  ⟨r⟩/⟨r⟩0  ⟨u_r⟩/U  ⟨p·r̂⟩  min ψ  |α̂_m|  |ρ̂_m|
# where R is the radius containing 99% of the mass, ⟨r⟩ the mass-weighted mean radius,
# ⟨u_r⟩ the mass-weighted radial velocity (self-propulsion plus surface flow), ⟨p·r̂⟩ the
# mass-weighted radial order, and |α̂_m| and |ρ̂_m| the amplitudes of azimuthal mode m of the
# polarization tilt α = arg(m e^{-iφ}) and of the density, both mass-weighted. With
# `plot=true`, the density, the local radial order m·r̂/ρ and the tilt α are displayed at
# every saved frame. A tilt changes m·r̂/ρ only at second order, so perturbations to the
# orientation show up in α.

include(joinpath(@__DIR__, "meanfield.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "script_params.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "farfield.jl"))

using Plots

aspect_ratio = script_param(:aspect_ratio, 1.0)
t_end        = script_param(:t_end, 200.0)
Nparticles   = script_param(:Nparticles, 10)
R0           = script_param(:R0, 1.0)      # rings at 5, 10, 15 with spacing 5
Rin          = script_param(:Rin, 0.1)      # annulus avoids the orientation singularity at r = 0
L            = script_param(:L, 5.0)
N            = script_param(:N, 256)
Nθ           = script_param(:Nθ, 32)
ζ            = script_param(:ζ, 5.0)
dT           = script_param(:dT, 1e-6)
dR           = script_param(:dR, 1e-6)
ε            = script_param(:ε, 1e-1)
nmodes       = script_param(:nmodes, 24)
seed         = script_param(:seed, 1)
perturb      = script_param(:perturb, :azimuthal)   # or :azimuthal, :modes
m            = script_param(:m, 5)                      # azimuthal mode for :azimuthal
tplt         = script_param(:tplt, 0.01)
outdir       = script_param(:outdir, joinpath(@__DIR__, "output", "patch_$(aspect_ratio)"))
show_plots   = script_param(:plot, false)

U, η, β = farfield_constants(aspect_ratio)
U = 1.0
η = 100.0
G  = SpectralGrid(N, N, Nθ; Lx=L, Ly=L)
ρ0 = Nparticles / (π * (R0^2 - Rin^2))

# perturbations on the scale of the patch: wavelengths L/nmax to L/nmin
nmin, nmax = round(Int, L / (2R0)), round(Int, L / 5)
if perturb == :orientation
    waves = random_modes(nmodes; nmin=nmin, nmax=nmax, nθmax=0, ε=ε, seed=seed)
    δθ(x, y) = sum(w -> w.ε * cos(2π * w.nx * x / L + 2π * w.ny * y / L + w.phase), waves)
    ψ0 = ic_patch(G; ρ0=ρ0, R=R0, Rin=Rin, width=2G.dx, ζ=ζ, δθ=δθ)
elseif perturb == :azimuthal
    ψ0 = ic_patch(G; ρ0=ρ0, R=R0, Rin=Rin, width=2G.dx, ζ=ζ,
                  δθ=(x, y) -> ε * cos(m * atan(y - L/2, x - L/2)))
else
    modes = random_modes(nmodes; nmin=nmin, nmax=nmax, nθmax=2, ε=ε, seed=seed)
    ψ0 = ic_patch(G; ρ0=ρ0, R=R0, Rin=Rin, width=2G.dx, ζ=ζ, modes=modes)
end

@printf("a/b = %g   U = %.6f   η = %.6f   β = %.6f   ρ0 = %.4f   dx = %.3f\n", aspect_ratio, U, η, β, ρ0, G.dx)

isdir(outdir) && rm(outdir; recursive=true)
mkpath(outdir)

xc, yc = G.Lx / 2, G.Ly / 2
rr = [hypot(G.x[i] - xc, G.y[j] - yc) for i in 1:G.Nx, j in 1:G.Ny]
er = [((G.x[i] - xc) + im * (G.y[j] - yc)) / max(rr[i,j], eps()) for i in 1:G.Nx, j in 1:G.Ny]

"Radius of the disk about the center containing a fraction `q` of the mass."
function mass_radius(ρ; q=0.99)
    k = sortperm(vec(rr))
    c = cumsum(vec(ρ)[k])
    return vec(rr)[k][findfirst(c .>= q * c[end])]
end

R_init, r_init = Ref(NaN), Ref(NaN)
io = open(joinpath(outdir, "patch.dat"), "w")
println(io, "# t  R/R0  <r>/<r>0  <u_r>/U  <p.r>  minψ  |α_m|  |ρ_m|")

"Tilt of the polarization from the outward radial direction, masked where `ρ` is small."
function tilt(snap)
    α = angle.(conj.(er) .* (snap.m1 .+ im .* snap.m2))
    α[snap.ρ .< 0.05 * maximum(snap.ρ)] .= NaN
    return α
end

function diagnostics(snap, _)
    R = mass_radius(snap.ρ)
    r̄ = sum(snap.ρ .* rr) / sum(snap.ρ)
    isnan(R_init[]) && (R_init[] = R; r_init[] = r̄)
    mr  = real.(conj.(er) .* (snap.m1 .+ im .* snap.m2))     # m · r̂
    ur  = real.(conj.(er) .* (snap.u1 .+ im .* snap.u2))     # u · r̂
    M   = sum(snap.ρ)
    order = sum(mr) / M
    urad  = (U * sum(mr) + sum(snap.ρ .* ur)) / M / U
    α     = tilt(snap)
    keep  = .!isnan.(α)
    em    = conj.(er[keep]) .^ m                                 # e^{-imφ}
    α̂m    = 2 * abs(sum(snap.ρ[keep] .* α[keep] .* em)) / sum(snap.ρ[keep])
    ρ̂m    = 2 * abs(sum(snap.ρ .* conj.(er) .^ m)) / M
    @printf(io, "%.4f  %.6f  %.6f  %.6f  %.6f  %+.3e  %.6e  %.6e\n", snap.t, R / R_init[], r̄ / r_init[], urad, order, snap.minψ, α̂m, ρ̂m)
    flush(io)
    show_plots && display(plot_patch(snap, 1.3R))
end

"Density, local radial order `m·r̂/ρ` and tilt `α` (masked where `ρ` is small) within radius `w`."
function plot_patch(snap, w)
    x, y = G.x .- xc, G.y .- yc
    order = real.(conj.(er) .* (snap.m1 .+ im .* snap.m2)) ./ snap.ρ
    order[snap.ρ .< 0.05 * maximum(snap.ρ)] .= NaN
    α  = tilt(snap)
    αmax = max(maximum(abs, filter(!isnan, α)), eps())
    format = (aspect_ratio=:equal, xlims=(-w, w), ylims=(-w, w), framestyle=:box)
    pρ = heatmap(x, y, snap.ρ'; format..., c=:dense, title=@sprintf("ρ, t = %.2f", snap.t))
    po = heatmap(x, y, order'; format..., c=:viridis, clims=(0, 1), title="m·r̂/ρ")
    pα = heatmap(x, y, α'; format..., c=:balance, clims=(-αmax, αmax), title="tilt α")
    return plot(pρ, po, pα; layout=(1, 3), size=(1500, 450))
end

simulate(ψ0, G; tf=t_end, U=U, β=β, η=η, dT=dT, dR=dR, tplt=tplt,
         output_dir=outdir, callback=diagnostics)
close(io)

show_plots && readline()
