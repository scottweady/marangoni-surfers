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
M  = 32       # Stokes expansion order

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

    Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)
    inside(z) = (real(z) / a)^2 + (imag(z) / b)^2 < 1

    xcm = [-Δx/2 + 0im, Δx/2 + 0im]
    θcm = [π/2, π/2]

    # -- Solve ---------------------------------------------------------------

    U, σ, Γ, ϕinf, f = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)
    @printf("  U[1] = %s\n", U[1])

    # -- Evaluate concentration and velocity on grid ---------------------------

    xplt = range(minimum(real.(xcm)) - pad, maximum(real.(xcm)) + pad, length=Nplt)
    yplt = range(minimum(imag.(xcm)) - pad, maximum(imag.(xcm)) + pad, length=Nplt)
    z    = ComplexF64[xn + im*yn for yn in yplt, xn in xplt]

    ϕ, u, v = surface_fields(z, xcm, θcm, Ω, σ, Γ, ϕinf, f; inside=inside)

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
