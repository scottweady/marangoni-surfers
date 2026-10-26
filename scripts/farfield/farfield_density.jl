# Coarse-grained number density of a far-field expansion run (farfield_expansion.jl), from a
# Gaussian kernel of width `ℓ` about each particle center, normalized by the initial lattice
# density 1/d², and its azimuthal average about the center of mass.
#
# Output: <framedir>/density_<n>.dat for each frame
#   Line 1: # Nx=... Ny=... t=...
#   Line 2: x values (Nx floats)
#   Line 3: y values (Ny floats)
#   Lines 4..(3+Ny): ρ matrix, Ny rows × Nx columns
# and <framedir>/radial_density.dat, one row per frame
#   t  ρ(r₁) ρ(r₂) ...   (first row: NaN followed by the radii r)

using DelimitedFiles
using Printf

include("../../src/script_params.jl")

aspect_ratio = script_param(:aspect_ratio, 1.0)
d            = script_param(:d, 8.0)
ℓ            = script_param(:ℓ, 2 * d)      # kernel width
h            = script_param(:h, d / 2)      # grid spacing
framedir     = script_param(:framedir, joinpath(@__DIR__, "output", "expansion_$(aspect_ratio)_d$(d)"))

files = sort(filter(f -> startswith(f, "frame_") && endswith(f, ".dat"), readdir(framedir)))
frames = [readdlm(joinpath(framedir, f)) for f in files]    # columns: t x y θ ux uy ω

# the grid covers the last (largest) frame, so all frames share it
last = frames[end]
pad  = 4ℓ
xg   = range(minimum(last[:, 2]) - pad, maximum(last[:, 2]) + pad; step=h)
yg   = range(minimum(last[:, 3]) - pad, maximum(last[:, 3]) + pad; step=h)
rg   = range(0, maximum(hypot.(last[:, 2] .- sum(last[:, 2]) / size(last, 1),
                                last[:, 3] .- sum(last[:, 3]) / size(last, 1))) + pad; step=h)

function density(x, y)
    ρ = zeros(length(yg), length(xg))
    w = ceil(Int, 4ℓ / h)
    for (xi, yi) in zip(x, y)
        ix = round(Int, (xi - first(xg)) / h) + 1
        iy = round(Int, (yi - first(yg)) / h) + 1
        for jy in max(1, iy - w):min(length(yg), iy + w), jx in max(1, ix - w):min(length(xg), ix + w)
            ρ[jy, jx] += exp(-((xg[jx] - xi)^2 + (yg[jy] - yi)^2) / (2ℓ^2))
        end
    end
    return ρ .* d^2 ./ (2π * ℓ^2)
end

# azimuthal average of ρ in annuli of width h about (xc, yc)
function radial(ρ, xc, yc)
    s, c = zeros(length(rg)), zeros(length(rg))
    for jy in eachindex(yg), jx in eachindex(xg)
        k = round(Int, hypot(xg[jx] - xc, yg[jy] - yc) / h) + 1
        k <= length(rg) || continue
        s[k] += ρ[jy, jx]
        c[k] += 1
    end
    return s ./ max.(c, 1)
end

radial_rows = [vcat(NaN, collect(rg))']
for (n, F) in enumerate(frames)
    t = F[1, 1]
    ρ = density(F[:, 2], F[:, 3])
    push!(radial_rows, vcat(t, radial(ρ, sum(F[:, 2]) / size(F, 1), sum(F[:, 3]) / size(F, 1)))')

    open(joinpath(framedir, "density_$(lpad(n - 1, 6, '0')).dat"), "w") do io
        @printf(io, "# Nx=%d Ny=%d t=%.4f\n", length(xg), length(yg), t)
        println(io, join([@sprintf("%.6e", v) for v in xg], " "))
        println(io, join([@sprintf("%.6e", v) for v in yg], " "))
        for iy in eachindex(yg)
            println(io, join([@sprintf("%.6e", ρ[iy, ix]) for ix in eachindex(xg)], " "))
        end
    end
end
writedlm(joinpath(framedir, "radial_density.dat"), reduce(vcat, radial_rows))
println("Wrote ", length(frames), " density frames to ", framedir)
