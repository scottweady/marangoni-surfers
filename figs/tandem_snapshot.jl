# Surfactant concentration around the tandem chains (scripts/suspension/tandem.jl) at snapshot
# times, for plotting in MATLAB (fig_tandem.m). At the saved frame nearest each time the full
# surfactant problem is solved for the particle configuration and the concentration is
# evaluated on a grid in the frame of the chain centroid. All snapshots share the x range;
# the y range of each aspect ratio covers both of its snapshots, and at least `yh_min`.
#
# Output: <outdir>/tandem_snapshot_<aspect_ratio>_t<t>.dat
#   Line 1: # Nx=... Ny=... t=...
#   Line 2: x - x_c values (Nx floats)
#   Line 3: y - y_c values (Ny floats)
#   Lines 4..(3+Ny): ϕ matrix, Ny rows × Nx columns (NaN inside particles)

using DelimitedFiles, Statistics

using MarangoniSurfers
using Printf

aspect_ratios = script_param(:aspect_ratios, [0.25, 1.0, 4.0])
t_snaps       = script_param(:t_snaps, [0.0, 1000.0])
h_grid        = script_param(:h_grid, 0.5)       # grid spacing
yh_min        = script_param(:yh_min, 40.0)      # least half-height of the grid
tandemdir     = script_param(:tandemdir, joinpath(@__DIR__, "..", "data", "tandem"))   # holds the runs
outdir        = script_param(:outdir, tandemdir)

s₀   = π/2
g(s) = -cos(s)
Ns   = 64

datadir(ab) = joinpath(tandemdir, "tandem_$(ab)_N50_d10.0_m5_ε0.01")
function load_snapshot(ab, t_snap)
    files  = sort(filter(f -> startswith(f, "frame_") && endswith(f, ".dat"), readdir(datadir(ab))))
    frames = [readdlm(joinpath(datadir(ab), f)) for f in files]   # columns: t x y θ vx vy ω
    return frames[argmin([abs(F[1,1] - t_snap) for F in frames])]
end

snaps = Dict((ab, t) => load_snapshot(ab, t) for ab in aspect_ratios, t in t_snaps)
xc(F) = mean(F[:,2]) + im * mean(F[:,3])

pad = 4 * maximum(sqrt.(aspect_ratios))
xl  = extrema(vcat([F[:,2] .- real(xc(F)) for F in values(snaps)]...)) .+ (-pad, pad)
xg  = range(xl..., step=h_grid)

mkpath(outdir)
for ab in aspect_ratios
    a = sqrt(ab);  b = 1 / a
    Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)
    yh = max(yh_min, maximum(maximum(abs.(snaps[(ab, t)][:,3] .- imag(xc(snaps[(ab, t)])))) for t in t_snaps) + pad)
    yg = range(-yh, yh, step=h_grid)
    for t_snap in t_snaps
        F   = snaps[(ab, t_snap)]
        xcm = F[:,2] .+ im .* F[:,3]
        θcm = F[:,4]
        σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)
        ϕ = concentration_grid(xg .+ real(xc(F)), yg .+ imag(xc(F)), xcm, θcm, σ, ϕinf, Γ, a, b)

        outfile = joinpath(outdir, @sprintf("tandem_snapshot_%s_t%d.dat", ab, round(Int, t_snap)))
        open(outfile, "w") do io
            @printf(io, "# Nx=%d Ny=%d t=%.4f\n", length(xg), length(yg), F[1,1])
            println(io, join([@sprintf("%.4f", v) for v in xg], " "))
            println(io, join([@sprintf("%.4f", v) for v in yg], " "))
            for iy in eachindex(yg)
                println(io, join([@sprintf("%.5f", ϕ[iy, ix]) for ix in eachindex(xg)], " "))
            end
        end
        println("Wrote ", relpath(outfile, @__DIR__))
    end
end
