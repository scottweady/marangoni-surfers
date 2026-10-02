"""
Evaluate the surfactant concentration field for a snapshot of the patch-expansion
run (data/patch) and write it to a dat file for plotting in MATLAB (fig_expansion).

Output: data/patch/expansion_snapshot.dat
  Line 1: # Nx=... Ny=... t=...
  Line 2: x values (Nx floats)
  Line 3: y values (Ny floats)
  Lines 4..(3+Ny): phi matrix, Ny rows × Nx columns (NaN inside particles)
"""

using DelimitedFiles

include("../src/marangoni.jl")

t_snap = 250.0

# ── Geometry (unit disks, as in data/patch) ──────────────────────────────────

a, b  = 1.0, 1.0
α_arc = π/2
g(θ)  = -cos(θ)

Ns = 64

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Γ_D = openarc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = openarc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Load frame nearest t_snap ────────────────────────────────────────────────

dataDir = joinpath(@__DIR__, "..", "data", "patch")
files   = sort(filter(f -> startswith(f, "frame_") && endswith(f, ".dat"), readdir(dataDir)))
frames  = [readdlm(joinpath(dataDir, f)) for f in files]    # columns: t x y θ vx vy ω
k       = argmin([abs(F[1,1] - t_snap) for F in frames])
A       = frames[k]
t_snap  = A[1,1]

xcm = A[:,2] .+ im .* A[:,3]
θcm = A[:,4]
Nparticles = length(xcm)
@printf("Snapshot at t = %.2f (%s), %d particles\n", t_snap, files[k], Nparticles)

# ── Solve surfactant and evaluate on a grid ──────────────────────────────────

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

# ── Write ────────────────────────────────────────────────────────────────────

outfile = joinpath(dataDir, "expansion_snapshot.dat")
open(outfile, "w") do io
    @printf(io, "# Nx=%d Ny=%d t=%.4f\n", Nplt, Nplt, t_snap)
    println(io, join([@sprintf("%.10e", v) for v in xplt], " "))
    println(io, join([@sprintf("%.10e", v) for v in yplt], " "))
    for iy in 1:Nplt
        println(io, join([@sprintf("%.10e", real(ϕ[iy, ix])) for ix in 1:Nplt], " "))
    end
end

println("Written to $outfile")
