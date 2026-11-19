# Example: a single elliptical surfer.
#
#     julia --project=scripts scripts/single/single.jl [aspect_ratio=4.0] [plot=true]
#
# Solves for the surfer's rigid-body velocity, checks it against the exact propulsion speed,
# and prints its multipole moments. With plot=true, the surfactant concentration around it
# is saved to scripts/single/output/.

using MarangoniSurfers
using Printf

# -- Parameters ----------------------------------------------------------------

aspect_ratio = script_param(:aspect_ratio, 4.0)   # a/b, with a along the propulsion axis
show_plots   = script_param(:plot, false)
show_plots && @eval using Plots                   # plotting comes from the Plots extension

a = sqrt(aspect_ratio)    # semi-axes of a unit-area ellipse, ab = 1
b = 1 / a
s₀   = π/2                # the Dirichlet arc is s ∈ (s₀, 2π - s₀), the rear of the particle
g(s) = -cos(s)            # surfactant concentration imposed on the Dirichlet arc

# -- Discretization ------------------------------------------------------------

Ns = 64    # quadrature nodes per surfactant arc
M  = 16    # Stokes expansion order on the particle

# Ω: the particle's Stokes domain; Γ_D, Γ_N: its Dirichlet and Neumann arcs, in the body frame
Ω, Γ_D, Γ_N = surfer(a, b, s₀, Ns, M)

# -- Solve ---------------------------------------------------------------------

println("Single surfer, a/b = $aspect_ratio")

U, σ, Γ, ϕinf = marangoni_solve(Ω, Γ_D, Γ_N, g)
U_exact, = farfield_constants(aspect_ratio; s₀=s₀)
@printf("  velocity (ux, uy, ω) = (%.8f, %.8f, %.8f)\n", U...)
@printf("  exact speed          = %.8f   (relative error %.1e)\n", U_exact,
        abs(hypot(U[1], U[2]) - U_exact) / U_exact)

mp = marangoni_multipoles(Ω, Γ_D, Γ_N, g)
@printf("  force dipole         = [%9.5f %9.5f; %9.5f %9.5f]\n", mp.force_dipole'...)
@printf("  surfactant monopole  = (%.5f, %.5f)\n", mp.surf_monopole...)
@printf("  surfactant dipole    = [%9.5f %9.5f; %9.5f %9.5f]\n", mp.surf_dipole'...)

# -- Concentration plot --------------------------------------------------------

if show_plots
    outdir = mkpath(joinpath(@__DIR__, "output"))
    p = plot_concentration([0.0im], σ, ϕinf, Γ; title="single surfer, a/b = $aspect_ratio")
    savefig(p, joinpath(outdir, "single.png"))
    println("\nPlot saved to $(joinpath(outdir, "single.png"))")
end
