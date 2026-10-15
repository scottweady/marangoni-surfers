# Mean-field model (meanfield.jl) of an aligned suspension on a doubly periodic box, for
# oblate, disk-shaped and prolate particles started from the same initial condition.
# Orientations are von Mises distributed about e_x with uniform density, perturbed by a
# single plane wave with mode numbers (nx, ny): with `perturb=:orientation` the local mean
# orientation is tilted by ε cos(k·x), leaving the density unperturbed; with
# `perturb=:density` the density is multiplied by 1 + ε cos(k·x), leaving the orientation
# distribution unperturbed.
#
# Output: <outdir>/aligned_<a/b>.dat, one row per frame
#   t  std(ρ)/⟨ρ⟩  rms|u|/U  ⟨m_x⟩/⟨ρ⟩  |ρ̂_k|/⟨ρ⟩  |m̂_y,k|/⟨ρ⟩  min ψ
# where ⟨·⟩ is the box average and ·̂_k the amplitude of the perturbed mode;
# <outdir>/final_<a/b>.dat, the fields at the final time, one row per grid point
#   x  y  ρ  m1  m2  u1  u2
# <outdir>/aligned.png, the diagnostics against time for every aspect ratio; and
# <outdir>/final.png, the final density and polarization for every aspect ratio; and
# <outdir>/snapshots_<a/b>.png, the density and polarization every `tsnap`, whose fields
# are written to <outdir>/snap_<a/b>_<n>.dat in the format of final_<a/b>.dat.

include(joinpath(@__DIR__, "meanfield.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "script_params.jl"))

using Plots

aspect_ratios = script_param(:aspect_ratios, [0.25, 1.0, 4.0])
t_end         = script_param(:t_end, 20.0)
U             = script_param(:U, 1.0)
η             = script_param(:η, 10.0)
ρ0            = script_param(:ρ0, 1.0)
L             = script_param(:L, 2π)
N             = script_param(:N, 64)
Nθ            = script_param(:Nθ, 32)
ζ             = script_param(:ζ, 10.0)
dT            = script_param(:dT, 0.005)
dR            = script_param(:dR, 0.005)
nx            = script_param(:nx, 0)
ny            = script_param(:ny, 4)
ε             = script_param(:ε, 1e-2)
perturb       = script_param(:perturb, :orientation)   # or :density
tplt          = script_param(:tplt, 0.05)
tsnap         = script_param(:tsnap, t_end / 5)
outdir        = script_param(:outdir, joinpath(@__DIR__, "output", "aligned"))
show_plots    = script_param(:plot, false)

G = SpectralGrid(N, N, Nθ; Lx=L, Ly=L)
wave(x, y) = cos(2π * nx * x / L + 2π * ny * y / L)
ψ0 = perturb == :orientation ?
     ic_vonmises(G, (x, y) -> ρ0, (x, y) -> ε * wave(x, y); ζ=ζ) :
     ic_vonmises(G, (x, y) -> ρ0 * (1 + ε * wave(x, y)), (x, y) -> 0.0; ζ=ζ)

isdir(outdir) && rm(outdir; recursive=true)
mkpath(outdir)

ek = [exp(-im * (2π * nx * G.x[i] / L + 2π * ny * G.y[j] / L)) for i in 1:G.Nx, j in 1:G.Ny]
mode_amplitude(f) = 2 * abs(sum(f .* ek)) / length(f)

"Write the fields of `s`, one row per grid point."
function write_fields(io, s)
    @printf(io, "# t = %.8f   Nx = %d   Ny = %d\n", s.t, G.Nx, G.Ny)
    println(io, "# x y rho m1 m2 u1 u2")
    for j in 1:G.Ny, i in 1:G.Nx
        @printf(io, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", G.x[i], G.y[j], s.ρ[i,j], s.m1[i,j], s.m2[i,j], s.u1[i,j], s.u2[i,j])
    end
end

results = Dict{Float64, Matrix{Float64}}()
finals  = Dict{Float64, NamedTuple}()
snaps   = Dict{Float64, Vector{NamedTuple}}()
for aspect_ratio in aspect_ratios
    a = sqrt(aspect_ratio)
    β = (a^2 - 1/a^2) / (a^2 + 1/a^2)
    @printf("a/b = %g   β = %+.4f   k = (%d, %d)   %s perturbation\n", aspect_ratio, β, nx, ny, perturb)

    rows = Vector{Float64}[]
    snaps[aspect_ratio] = NamedTuple[]
    function diagnostics(snap, _)
        ρ̄ = sum(snap.ρ) / length(snap.ρ)
        push!(rows, [snap.t,
                     sqrt(sum(abs2, snap.ρ .- ρ̄) / length(snap.ρ)) / ρ̄,
                     sqrt(sum(abs2, snap.u1) + sum(abs2, snap.u2)) / sqrt(length(snap.u1)) / U,
                     sum(snap.m1) / sum(snap.ρ),
                     mode_amplitude(snap.ρ) / ρ̄,
                     mode_amplitude(snap.m2) / ρ̄,
                     snap.minψ])
        finals[aspect_ratio] = snap
        snap.t >= length(snaps[aspect_ratio]) * tsnap - 1e-9 && push!(snaps[aspect_ratio], snap)
    end
    simulate(ψ0, G; tf=t_end, U=U, β=β, η=η, dT=dT, dR=dR, tplt=tplt,
             callback=diagnostics, verbose=false)

    D = reduce(vcat, permutedims.(rows))
    results[aspect_ratio] = D
    open(joinpath(outdir, "aligned_$(aspect_ratio).dat"), "w") do io
        println(io, "# t  std(rho)/<rho>  rms|u|/U  <m_x>/<rho>  |rho_k|/<rho>  |m_y,k|/<rho>  min psi")
        foreach(r -> println(io, join([@sprintf("%.6e", v) for v in r], "  ")), eachrow(D))
    end
    s = finals[aspect_ratio]
    open(joinpath(outdir, "final_$(aspect_ratio).dat"), "w") do io
        write_fields(io, s)
    end
    for (n, sn) in enumerate(snaps[aspect_ratio])
        open(io -> write_fields(io, sn), joinpath(outdir, "snap_$(aspect_ratio)_$(lpad(n-1, 3, '0')).dat"), "w")
    end
    @printf("  t = %g   std(ρ)/⟨ρ⟩ = %.3e   rms|u|/U = %.3e   ⟨m_x⟩/⟨ρ⟩ = %.4f   min ψ over run = %+.3e\n", D[end,1], D[end,2], D[end,3], D[end,4], minimum(D[:,7]))
end

titles = ["std(ρ)/⟨ρ⟩", "rms|u|/U", "⟨m_x⟩/⟨ρ⟩", "|m̂_y,k|/⟨ρ⟩"]
cols   = [2, 3, 4, 6]
logy   = [true, true, false, true]
panels = map(1:4) do n
    p = plot(; xlabel="t", title=titles[n], yscale=logy[n] ? :log10 : :identity, legend=n == 1 ? :bottomright : false)
    for aspect_ratio in aspect_ratios
        D = results[aspect_ratio]
        y = logy[n] ? max.(D[:, cols[n]], 1e-16) : D[:, cols[n]]
        plot!(p, D[:, 1], y; label="a/b = $(aspect_ratio)", lw=1.5)
    end
    p
end
fig = plot(panels...; layout=(2, 2), size=(1000, 750))
savefig(fig, joinpath(outdir, "aligned.png"))

"Final density with the local polarization m/ρ, every `s`-th grid point."
function plot_final(snap, title; s=4)
    is, js = 1:s:G.Nx, 1:s:G.Ny
    X = [G.x[i] for i in is, j in js];  Y = [G.y[j] for i in is, j in js]
    P1 = snap.m1[is, js] ./ snap.ρ[is, js];  P2 = snap.m2[is, js] ./ snap.ρ[is, js]
    sc = 0.8 * s * G.dx
    p = heatmap(G.x, G.y, snap.ρ' ./ ρ0; aspect_ratio=:equal, xlims=(0, L), ylims=(0, L),
                c=:dense, framestyle=:box, title=title, xlabel="x", ylabel="y")
    quiver!(p, vec(X), vec(Y); quiver=(sc .* vec(P1), sc .* vec(P2)), color=:black, lw=0.5)
    return p
end
pf = plot([plot_final(finals[r], @sprintf("ρ/ρ₀, a/b = %g, t = %.0f", r, finals[r].t)) for r in aspect_ratios]...;
          layout=(1, length(aspect_ratios)), size=(450 * length(aspect_ratios), 420))
savefig(pf, joinpath(outdir, "final.png"))
for r in aspect_ratios
    S = snaps[r]
    nc = ceil(Int, sqrt(length(S)));  nr = ceil(Int, length(S) / nc)
    ps = plot([plot_final(s, @sprintf("a/b = %g, t = %.0f", r, s.t)) for s in S]...;
              layout=(nr, nc), size=(450 * nc, 420 * nr))
    savefig(ps, joinpath(outdir, "snapshots_$(r).png"))
end
show_plots && (display(fig); display(pf); readline())
