# Surfactant concentration and surface velocity fields for a parallel pair of elliptical
# surfers separated by `Δx` along the x-axis, sweeping `aspect_ratios` at fixed area.
#
# Output: output/fields_a<a>_b<b>.dat, one file per aspect ratio
#   Line 1: # Nx=... Ny=... a=... b=... d=...
#   Line 2: x values (Nx floats)
#   Line 3: y values (Ny floats)
#   Lines 4..(3+Ny):        ϕ matrix (concentration), Ny rows × Nx columns
#   Lines (4+Ny)..(3+2Ny):  u matrix (x-velocity)
#   Lines (4+2Ny)..(3+3Ny): v matrix (y-velocity)

include("../../src/marangoni.jl")

# -- Parameters ----------------------------------------------------------------

s₀   = π/2
g(s) = -cos(s)

aspect_ratios = [0.25, 1.0, 4.0]
Δx = 8.0  # initial center-to-center separation

# -- Discretization ------------------------------------------------------------

Ns = 128      # arc quadrature points
M  = 16       # spectral degree for ellipse

# -- Field grid parameters -----------------------------------------------------

pad  = 10.0
Nplt = 512

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

# -- Main loop over aspect ratio -----------------------------------------------

for aspect_ratio in aspect_ratios

    a = sqrt(aspect_ratio)
    b = 1.0 / a

    @printf("(a, b) = (%.4f, %.4f)\n", a, b)

    # Define ellipse and arc boundaries (body frame, θ=0)
    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

    xcm = [-Δx/2 + 0im, Δx/2 + 0im]
    θcm = [π/2, π/2]

    # -- Solve: surfactant + Stokes mobility -----------------------------------

    println("  Solving surfactant...")
    σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

    Nparticles = length(xcm)
    Npts       = length(Ω.z)
    offs       = [0; cumsum(fill(Npts, Nparticles))]

    # ∇ℬ at Stokes quadrature nodes (uinf for the mobility problem)
    pts = [exp(im * θcm[i]) .* Ω.z .+ xcm[i] for i in 1:Nparticles]
    u_pts, v_pts = ∇ℬ(vcat(vec.(pts)...), σ, Γ)

    uinf = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Nparticles)
    for i in 1:Nparticles
        ri  = offs[i]+1:offs[i+1]
        shp = size(pts[i])
        uinf[i] = (reshape(copy(u_pts[ri]), shp), reshape(copy(v_pts[ri]), shp))
    end

    println("  Solving Stokes mobility...")
    F = [zeros(3) for _ in 1:Nparticles]
    f_body, U = suspension_mobility_solve(xcm, θcm, uinf, nothing, F, Ω)
    @printf("  U[1] = %s\n", U[1])

    # -- Evaluate concentration and velocity on grid ---------------------------

    xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
    yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
    z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    # Concentration
    ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in 1:Nparticles) .+ ϕinf

    # Full velocity: u = S[f](z) + ∇ℬ[σ](z)
    println("  Evaluating velocity field...")
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

    # -- Write results ---------------------------------------------------------

    function write_matrix(io, A, Nplt)
        for iy in 1:Nplt
            println(io, join([@sprintf("%.10e", A[iy, ix]) for ix in 1:Nplt], " "))
        end
    end

    outfile = joinpath(outdir, @sprintf("fields_a%.4f_b%.4f.dat", a, b))
    open(outfile, "w") do io
        @printf(io, "# Nx=%d Ny=%d a=%.6f b=%.6f d=%.1f\n", Nplt, Nplt, a, b, Δx)
        println(io, join([@sprintf("%.10e", xv) for xv in xplt], " "))
        println(io, join([@sprintf("%.10e", yv) for yv in yplt], " "))
        write_matrix(io, ϕ, Nplt)
        write_matrix(io, u, Nplt)
        write_matrix(io, v, Nplt)
    end

    println("  Wrote $outfile")

end

println("\nDone.")
