# Surfactant concentration for a snapshot of the patch-expansion run (data/patch), for
# plotting in MATLAB (fig_expansion.m).
#
# Output: output/expansion_snapshot.dat
#   Line 1: # Nx=... Ny=... t=...
#   Line 2: x values (Nx floats)
#   Line 3: y values (Ny floats)
#   Lines 4..(3+Ny): ϕ matrix, Ny rows × Nx columns (NaN inside particles)

using DelimitedFiles

include("../src/marangoni.jl")

t_snap = 250.0

# -- Geometry ------------------------------------------------------------------

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

Ns = 128

Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns)

# -- Load frame nearest t_snap -------------------------------------------------

datadir = joinpath(@__DIR__, "..", "data", "patch")
files   = sort(filter(f -> startswith(f, "frame_") && endswith(f, ".dat"), readdir(datadir)))
frames  = [readdlm(joinpath(datadir, f)) for f in files]    # columns: t x y θ vx vy ω
k       = argmin([abs(F[1,1] - t_snap) for F in frames])
A       = frames[k]
t_snap  = A[1,1]

xcm = A[:,2] .+ im .* A[:,3]
θcm = A[:,4]
Nparticles = length(xcm)
@printf("Snapshot at t = %.2f, %d particles\n", t_snap, Nparticles)

# -- Solve surfactant and evaluate on a grid -----------------------------------

σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

pad  = 5.0
Nplt = 400
xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

for k in 1:Nparticles
    z_body = exp(-im * θcm[k]) .* (z .- xcm[k])
    ϕ[(real.(z_body) / a).^2 .+ (imag.(z_body) / b).^2 .< 1.0] .= NaN
end

# -- Write results -------------------------------------------------------------

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
outfile = joinpath(outdir, "expansion_snapshot.dat")
open(outfile, "w") do io
    @printf(io, "# Nx=%d Ny=%d t=%.4f\n", Nplt, Nplt, t_snap)
    println(io, join([@sprintf("%.10e", v) for v in xplt], " "))
    println(io, join([@sprintf("%.10e", v) for v in yplt], " "))
    for iy in 1:Nplt
        println(io, join([@sprintf("%.10e", real(ϕ[iy, ix])) for ix in 1:Nplt], " "))
    end
end

println("Wrote ", relpath(outfile, @__DIR__))
