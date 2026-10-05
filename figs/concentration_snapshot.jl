# Surfactant concentration and surface velocity for a snapshot of a pair trajectory,
# for plotting in MATLAB.
#
# Output: concentration_snapshot.dat
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

M  = 16
Ns = 128

Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Load snapshot from trajectory dat file ------------------------------------

datfile = joinpath(@__DIR__, "../data/disk", "a1.0000_b1.0000_d16.0.dat")
A = readdlm(datfile)

# Use midpoint of trajectory as representative snapshot
row = A[end-15, :]
t_snap = row[1]
xcm = [row[2] + im*row[3], row[5] + im*row[6]]
θcm = [row[4], row[7]]

@printf("Snapshot at t = %.2f\n", t_snap)
@printf("  xcm = [%.4f%+.4fi, %.4f%+.4fi]\n", real(xcm[1]), imag(xcm[1]), real(xcm[2]), imag(xcm[2]))
@printf("  θcm = [%.4f, %.4f]\n", θcm[1], θcm[2])

# -- Solve: surfactant + Stokes mobility ---------------------------------------

println("Solving surfactant...")
σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

Nparticles = length(xcm)
Npts       = length(Ω.z)
offs       = [0; cumsum(fill(Npts, Nparticles))]

# ∇ℬ at Stokes quadrature nodes (uinf for the mobility problem)
pts = [exp(im * θcm[i]) .* Ω.z .+ xcm[i] for i in 1:Nparticles]
u_pts, v_pts = ∇ℬ(vcat(vec.(pts)...), σ, Γ)

uinf = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Nparticles)
for i in 1:Nparticles
    ri = offs[i]+1:offs[i+1]
    shp = size(pts[i])
    uinf[i] = (reshape(copy(u_pts[ri]), shp), reshape(copy(v_pts[ri]), shp))
end

println("Solving Stokes mobility...")
F = [zeros(3) for _ in 1:Nparticles]
f_body, U = suspension_mobility_solve(xcm, θcm, uinf, nothing, F, Ω)
println("Done. U[1] = $(U[1])")

# -- Evaluate concentration and velocity on grid -------------------------------

pad  = 5.0
Nplt = 200
xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

# Concentration
ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

# Full velocity: u = S[f](z) + ∇ℬ[σ](z)
println("Evaluating velocity field...")
u, v = ∇ℬ(z, σ, Γ)   # Marangoni component

for k in 1:Nparticles
    # Transform grid to body frame of particle k
    z_body = ComplexF64.(exp(-im * θcm[k]) .* (z .- xcm[k]))
    # Evaluate Stokes single-layer at body-frame exterior points
    vel = 𝒮_st(z_body, f_body[k], Ω)   # Array of (u_body, v_body) tuples
    u_b = real.(first.(vel))
    v_b = real.(last.(vel))
    # Rotate body-frame velocity back to world frame
    ct = cos(θcm[k]); st = sin(θcm[k])
    u .+= ct .* u_b .- st .* v_b
    v .+= st .* u_b .+ ct .* v_b
end

# Mask points inside either particle
for k in 1:Nparticles
    z_body = exp(-im * θcm[k]) .* (z .- xcm[k])
    inside = (real.(z_body) / a).^2 .+ (imag.(z_body) / b).^2 .< 1.0
    ϕ[inside] .= NaN
    u[inside] .= NaN
    v[inside] .= NaN
end

# -- Write results -------------------------------------------------------------

function write_matrix(io, M, Nplt)
    for iy in 1:Nplt
        println(io, join([@sprintf("%.10e", M[iy, ix]) for ix in 1:Nplt], " "))
    end
end

ucm_vel = mean(U[i][1] for i in 1:Nparticles)
vcm_vel = mean(U[i][2] for i in 1:Nparticles)

outfile = joinpath(@__DIR__, "concentration_snapshot.dat")
open(outfile, "w") do io
    @printf(io, "# Nx=%d Ny=%d t=%.4f ucm=%.10e vcm=%.10e\n", Nplt, Nplt, t_snap, ucm_vel, vcm_vel)
    println(io, join([@sprintf("%.10e", v) for v in xplt], " "))
    println(io, join([@sprintf("%.10e", v) for v in yplt], " "))
    write_matrix(io, ϕ, Nplt)
    write_matrix(io, u, Nplt)
    write_matrix(io, v, Nplt)
end

println("Written to $outfile")
