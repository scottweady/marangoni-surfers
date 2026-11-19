"""
    MarangoniSurfers

Boundary integral and far-field simulation of Marangoni surfers: rigid particles on a fluid
interface, self-propelled by the surface tension gradient of the surfactant they release.

The module is organized in layers:
- `bie/`: 2D boundary integral tools — quadrature, boundary representation, kernel matrix
  assembly, layer potentials, and close evaluation.
- `solvers.jl`: surfactant and Marangoni flow solves for a suspension of elliptical surfers.
- `simulation.jl`: adaptive time integration of the full (boundary integral) model.
- `farfield.jl`: far-field (multipole) model of the suspension.
- `plotting.jl`: surfactant concentration plots.
- `script_params.jl`: command-line parameter overrides for the scripts in `scripts/`.
"""
module MarangoniSurfers

using FastGaussQuadrature
using FFTW
using FMM3D
using LinearAlgebra
using Plots
using Printf
using ProjectedSphericalHarmonics
using SpecialFunctions

# Boundary integral tools. `dot` and `cross` (complex numbers as 2-vectors) are not exported,
# to avoid clashing with LinearAlgebra.
export Boundary, closedcurve, openarc, circle, segment, ellipse_arc, with_coupling
export integral, laplace_solve
export 𝒱, 𝒦, ℬ, ℋ, ∇𝒱, ∇𝒦, ∇ℬ, ∇ℋ, 𝒱_close, K′_close

# Marangoni surfers: geometry, solvers, and full simulation
export ellipse, surfer, surfer_boundary
export surfactant_solve, surfactant_response, surface_fields
export marangoni_solve, marangoni_self_operator, marangoni_multipoles
export rigid_velocities, suspension_simulate, midpoint_step, step_size, min_gap, boundary_gap

# Far-field model
export farfield_constants, farfield_moments, farfield_closure, farfield_response
export farfield_velocities, farfield_simulate, repulsion_velocities!, pair_bounds, ellipse_outline

# Plotting and scripts
export plot_concentration, concentration_grid
export script_param

include("script_params.jl")

# Boundary integral tools
include(joinpath("bie", "boundary.jl"))
include(joinpath("bie", "quadrature.jl"))
include(joinpath("bie", "matrices.jl"))
include(joinpath("bie", "operators.jl"))
include(joinpath("bie", "solvers.jl"))
include(joinpath("bie", "closeeval.jl"))

# Marangoni surfers
include("solvers.jl")
include("simulation.jl")
include("farfield.jl")
include("plotting.jl")

end
