
include(joinpath(@__DIR__, "meanfield.jl"))

using Plots, ColorSchemes
using SpecialFunctions: ellipk, ellipe

"""
    propulsion_ratio(a, b)

The geometric ratio `χ = U0/β = (m₀ - m₁)/(2π a b)` for an ellipse with semi-axis `a`
along the propulsion axis and semi-axis `b`.

The manuscript formulas assume `a ≥ b`; for an oblate particle (`a < b`) the propulsion
axis is the minor axis, so the axes are relabelled `1 ↔ 2` and the minor-axis coefficient
`m₂ = a/(2ε²)[E - (1-ε²)K]` with `ε² = 1 - (a/b)²` replaces `m₁`. Both branches share the
disk limit `m₁ - m₀ = -π min(a,b)/8`, taken directly to avoid the `0/0` form.
"""
function propulsion_ratio(a, b)
    Δm = if abs(a - b) < 1e-8 * (a + b)
        -π * min(a, b) / 8
    elseif a > b
        ε² = 1 - (b/a)^2
        b/(2ε²) * (ellipk(ε²) - ellipe(ε²)) - b/2 * ellipk(ε²)
    else
        ε² = 1 - (a/b)^2
        a/(2ε²) * (ellipe(ε²) - (1 - ε²) * ellipk(ε²)) - a/2 * ellipk(ε²)
    end
    return -Δm / (2π * a * b)
end

aspect_ratio = 8.0
a = sqrt(aspect_ratio)
b = 1 / sqrt(aspect_ratio)
U0 = 1.0
χ  = propulsion_ratio(a, b)
β  = U0/χ
λ  = (aspect_ratio^2 - 1) / (aspect_ratio^2 + 1)
display(β)
dT = 1e-3
dR = 0.01
L  = 2π

config = :isotropic # :isotropic or :patch
tplt = 0.1

# config = :patch # :isotropic or :patch
# tplt = 1e-3

N  = 128
Nθ = 32
t_end = 200.0

G  = MeanFieldGrid(N, N, Nθ; Lx=L, Ly=L)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

const frame_count = Ref(0)

function render(snap, G; nquiver=20, savedir="")

    mmag  = hypot.(snap.m1, snap.m2)
    ρsafe = max.(snap.ρ, 1e-12 * maximum(abs, snap.ρ))
    order = mmag ./ ρsafe
    mang  = atan.(snap.m2, snap.m1)

    common = (aspect_ratio = :equal,
              xlims = (0, G.Lx), ylims = (0, G.Ly),
              xlabel = "x", ylabel = "y")

    # scale from zero, not over the range: ρ is uniform to roundoff when ερ = 0
    pρ = heatmap(G.x, G.y, snap.ρ';
                 c = cgrad(:algae,rev=true), title = "ρ", clims=(0.5,1.5),common...)

    pm = heatmap(G.x, G.y, order';
                 c = :viridis, clims = (0, max(maximum(order), eps())),
                 title = "|m|/ρ", common...)

    pa = heatmap(G.x, G.y, mang';
                 c = :phase, clims = (-π, π), title = "arg m",
                 colorbar_ticks = ([-π, -π/2, 0, π/2, π],
                                   ["-π", "-π/2", "0", "π/2", "π"]),
                 common...)

    st = max(1, G.Nx ÷ nquiver)
    xs = Float64[]; ys = Float64[]; us = Float64[]; vs = Float64[]
    for j in 1:st:G.Ny, i in 1:st:G.Nx
        push!(xs, G.x[i]); push!(ys, G.y[j])
        push!(us, snap.m1[i,j]); push!(vs, snap.m2[i,j])
    end

    # p = plot(pρ, pm, pa; layout = (1, 3), size = (1500, 450), top_margin = 8Plots.mm, 
            # plot_title = "t = $(round(snap.t, digits=2))")
    p = plot(pρ)

    display(p)
    if !isempty(savedir)
        mkpath(savedir)
        savefig(p, joinpath(savedir, "frame_$(lpad(frame_count[], 6, '0')).png"))
    end
    frame_count[] += 1

    return p

end

if config == :isotropic
    ψ0 = ic_isotropic_perturbed(G; ρ0=1.0, modes = random_modes(12; nmin=1, nmax=4, εm=1e-2, seed=1))
elseif config == :patch
    R0 = 1.0;
    ρ0 = L^2/(π * R0^2)
    ψ0 = ic_patch(G; ρ0=ρ0, R=R0,width=10.0/N)
end

@printf("grid: %d x %d x %d on [0,%.4f]²,  dx = %.5f\n", N, N, Nθ, L, G.dx)

ψ, snaps = meanfield_simulate(ψ0, G; t_end=t_end, U0=U0, λ=λ, β=β, dT=dT, dR=dR, tplt=tplt,callback=(s, g) -> render(s, g; savedir=joinpath(outdir, "frames")))
