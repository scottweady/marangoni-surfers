# Surfactant concentration and surface velocity for a snapshot of a pair trajectory,
# for plotting in MATLAB.
#
# Output: output/concentration_snapshot.dat
#   Line 1: # Nx=... Ny=... t=... ucm=... vcm=...
#   Line 2: x values (Nx floats)
#   Line 3: y values (Ny floats)
#   Lines 4..(3+Ny):        ϕ matrix, Ny rows × Nx columns
#   Lines (4+Ny)..(3+2Ny):  u matrix (x-velocity)
#   Lines (4+2Ny)..(3+3Ny): v matrix (y-velocity)

using DelimitedFiles
using Statistics

include("../src/marangoni.jl")

# -- Geometry ------------------------------------------------------------------

aspect_ratio = 1.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

M  = 32
Ns = 128

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Load snapshot from trajectory dat file ------------------------------------

datfile = joinpath(@__DIR__, "../data/disk", "a1.0000_b1.0000_d16.0.dat")
A = readdlm(datfile)

# Snapshot shortly before the end of the trajectory
row = A[end-15, :]
t_snap = row[1]
xcm = [row[2] + im*row[3], row[5] + im*row[6]]
θcm = [row[4], row[7]]

@printf("Snapshot at t = %.2f\n", t_snap)
@printf("  xcm = [%.4f%+.4fi, %.4f%+.4fi]\n", real(xcm[1]), imag(xcm[1]), real(xcm[2]), imag(xcm[2]))
@printf("  θcm = [%.4f, %.4f]\n", θcm[1], θcm[2])

# -- Solve ---------------------------------------------------------------------

U, σ, Γ, ϕinf, f = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
println("U[1] = $(U[1])")

# -- Evaluate concentration and velocity on grid -------------------------------

pad  = 5.0
Nplt = 200
xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

inside(z) = (real(z) / a)^2 + (imag(z) / b)^2 < 1
ϕ, u, v = surface_fields(z, xcm, θcm, Ω, σ, Γ, ϕinf, f; inside=inside)

# -- Write results -------------------------------------------------------------

function write_matrix(io, M, Nplt)
    for iy in 1:Nplt
        println(io, join([@sprintf("%.10e", M[iy, ix]) for ix in 1:Nplt], " "))
    end
end

ucm_vel = mean(Ui[1] for Ui in U)
vcm_vel = mean(Ui[2] for Ui in U)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
outfile = joinpath(outdir, "concentration_snapshot.dat")
open(outfile, "w") do io
    @printf(io, "# Nx=%d Ny=%d t=%.4f ucm=%.10e vcm=%.10e\n", Nplt, Nplt, t_snap, ucm_vel, vcm_vel)
    println(io, join([@sprintf("%.10e", v) for v in xplt], " "))
    println(io, join([@sprintf("%.10e", v) for v in yplt], " "))
    write_matrix(io, ϕ, Nplt)
    write_matrix(io, u, Nplt)
    write_matrix(io, v, Nplt)
end

println("Written to $outfile")
