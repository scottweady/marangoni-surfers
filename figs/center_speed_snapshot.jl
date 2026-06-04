"""
Evaluate surfactant concentration and surface velocity for N=13 circles in a line,
for θ=0 and θ=π/2.  Writes one dat file per orientation.

Output files: center_speed_snapshot_th0.dat, center_speed_snapshot_thpi2.dat
  Line 1: # Nx=... Ny=... theta=...
  Line 2: x values (Nx floats)
  Line 3: y values (Ny floats)
  Lines 4..(3+Ny):         phi matrix, Ny rows × Nx columns
  Lines (4+Ny)..(3+2Ny):   u matrix
  Lines (4+2Ny)..(3+3Ny):  v matrix
"""


include("../src/marangoni.jl")

# ── Geometry: unit circles ────────────────────────────────────────────────────

a, b  = 1.0, 1.0
α_arc = π / 2
g(θ)  = -cos(θ)

M  = 8
Ns = 128

γ(s)   =  a*cos(s) + im*b*sin(s)
γs(s)  = -a*sin(s) + im*b*cos(s)
γss(s) = -a*cos(s) - im*b*sin(s)

Ω   = ellipse(a, b, M)
Γ_D = arc(γ, γs, γss, Ns; dom=[α_arc, 2π - α_arc])
Γ_N = arc(γ, γs, γss, Ns; dom=[-α_arc, α_arc])

# ── Configuration ─────────────────────────────────────────────────────────────

Nparticles = 13
d          = 20.0
xcm_base   = ComplexF64[d * (k - (Nparticles + 1) / 2) for k in 1:Nparticles]

# Grid: 4:1 aspect ratio (x-range / 4 = y-range) for collective flow visualization
pad    = 3.0
Nx     = 400
Ny     = 100
xmin   = real(xcm_base[1])   - pad
xmax   = real(xcm_base[end]) + pad
xrange = xmax - xmin
xplt   = range(xmin, xmax,        length=Nx)
yplt   = range(-xrange/8, xrange/8, length=Ny)   # half-height = xrange/8 → 4:1

# ── Helper: write a Ny×Nx matrix row by row ───────────────────────────────────

function write_matrix(io, M, Nx, Ny)
    for iy in 1:Ny
        println(io, join([@sprintf("%.10e", M[iy, ix]) for ix in 1:Nx], " "))
    end
end

# ── Loop over orientations ────────────────────────────────────────────────────

for (θ0, tag) in [(0.0, "th0"), (π/2, "thpi2")]

    xcm = xcm_base
    θcm = fill(θ0, Nparticles)

    @printf("\n─── θ = %.4f ───\n", θ0)

    # Surfactant solve
    println("  Solving surfactant...")
    σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

    # ∇ℬ at Stokes quadrature nodes → uinf for mobility problem
    Npts = length(Ω.z)
    offs = [0; cumsum(fill(Npts, Nparticles))]
    pts  = [exp(im * θcm[i]) .* Ω.z .+ xcm[i] for i in 1:Nparticles]
    u_pts, v_pts = ∇ℬ(vcat(vec.(pts)...), σ, Γ)

    uinf = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Nparticles)
    for i in 1:Nparticles
        ri      = offs[i]+1:offs[i+1]
        shp     = size(pts[i])
        uinf[i] = (reshape(copy(u_pts[ri]), shp), reshape(copy(v_pts[ri]), shp))
    end

    # Stokes mobility
    println("  Solving Stokes mobility...")
    F = [zeros(3) for _ in 1:Nparticles]
    f_body, U = suspension_mobility_solve(xcm, θcm, uinf, nothing, F, Ω)
    @printf("  U[center] = (%.6f, %.6f, %.6f)\n", U[(Nparticles+1)÷2]...)

    # Grid
    z = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    # Concentration
    println("  Evaluating concentration...")
    ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

    # Full velocity: u = ∇ℬ[σ] + S[f]
    println("  Evaluating velocity field...")
    u, v = ∇ℬ(z, σ, Γ)

    for k in 1:Nparticles
        z_body = ComplexF64.(exp(-im * θcm[k]) .* (z .- xcm[k]))
        vel    = 𝒮_st(z_body, f_body[k], Ω)
        u_b    = real.(first.(vel))
        v_b    = real.(last.(vel))
        ct     = cos(θcm[k]); st = sin(θcm[k])
        u .+= ct .* u_b .- st .* v_b
        v .+= st .* u_b .+ ct .* v_b
    end

    # Mask interior points
    for k in 1:Nparticles
        z_body  = exp(-im * θcm[k]) .* (z .- xcm[k])
        inside  = (real.(z_body) / a).^2 .+ (imag.(z_body) / b).^2 .< 1.0
        ϕ[inside] .= NaN
        u[inside]  .= NaN
        v[inside]  .= NaN
    end

    # Write
    outfile = joinpath(@__DIR__, "center_speed_snapshot_$(tag).dat")
    open(outfile, "w") do io
        @printf(io, "# Nx=%d Ny=%d theta=%.6f\n", Nx, Ny, θ0)
        println(io, join([@sprintf("%.10e", x) for x in xplt], " "))
        println(io, join([@sprintf("%.10e", y) for y in yplt], " "))
        write_matrix(io, ϕ, Nx, Ny)
        write_matrix(io, u, Nx, Ny)
        write_matrix(io, v, Nx, Ny)
    end
    println("  Written to $outfile")
end
