# MarangoniSurfers.jl

Boundary integral simulation of Marangoni surfers: rigid particles on a fluid interface,
self-propelled by the surface tension gradient of the surfactant they release. This is the
code for *A boundary integral formulation for microscale Marangoni surfing* (S. Weady).

The package solves, for a suspension of elliptical surfers,

- the surfactant problem: a mixed Dirichlet/Neumann Laplace problem on each particle's
  boundary arcs, by a boundary integral method (dense for a few particles; GMRES with
  Ewald-summed interactions and a far-field preconditioner for many);
- the Marangoni flow and the particles' rigid-body velocities, through the Stokes mobility
  solver of [ProjectedSphericalHarmonics.jl](https://github.com/scottweady/ProjectedSphericalHarmonics.jl);
- the time evolution of the suspension, by an adaptive midpoint method;

along with the far-field (multipole) model of the suspension.

## Installation

Requires Julia 1.12. The package is not registered; from a clone of this repository,

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()   # also fetches ProjectedSphericalHarmonics from GitHub ([sources])
```

## Usage

```julia
using MarangoniSurfers

g(s) = -cos(s)                                  # concentration on the Dirichlet arc
Ω, Γ_D, Γ_N = surfer(2.0, 0.5, π/2, 64, 16)     # a/b = 4 ellipse, arc junction at s₀ = π/2

# a single surfer: rigid-body velocity (ux, uy, ω)
U, = marangoni_solve(Ω, Γ_D, Γ_N, g)

# a suspension: velocities of each particle, for centers xcm and orientations θcm
xcm, θcm = [0.0im, 8.0 + 0.0im], [π/2, π/2]
U, σ, Γ, ϕ∞, f = marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g)

# time integration
traj_xcm, traj_θcm, traj_t = suspension_simulate(xcm, θcm, Ω, Γ_D, Γ_N, g;
                                                 t_end=10.0, render=false)
```

Plotting functions (`plot_concentration`, and `suspension_simulate(...; render=true)`) are
provided by a package extension and need `using Plots`.

The solvers are multithreaded; start Julia with threads, e.g. `julia -t 8`.

A worked single-surfer example, which checks the velocity against the exact propulsion speed,
prints the multipole moments and plots the concentration, is `scripts/single/single.jl`:

```sh
julia -t 8 --project=scripts scripts/single/single.jl plot=true
```

## Tests

```julia
using Pkg
Pkg.test("MarangoniSurfers")
```

or, from `test/`, `julia --project=. runtests.jl [--verbose]` in the test environment
(`test/Project.toml`), where `--verbose` runs and prints the full convergence study.

## Reproducing the figures

The scripts in `scripts/` and the figure data generators in `figs/` run in their own
environment, which adds Plots and the scripts' other dependencies, and whose
`scripts/Manifest.toml` pins every package version:

```sh
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts scripts/makedata.jl            # all data, into data/
julia --project=scripts scripts/makedata.jl tandem     # or one group: single, pairs, suspension, tandem
```

`makedata.jl` lists the groups and options (`resume`, `data=<dir>`). The figures are then
drawn in MATLAB with `figs/makefigs.m`, which reads `data/` and writes PDFs to
`figs/output/`.

The results in the paper were computed with the code tagged `paper-v1`.

## Layout

```
src/
├── MarangoniSurfers.jl   # module: dependencies, exports, includes
├── bie/                  # 2D boundary integral tools: boundaries, quadrature, layer
│                         # potentials, Laplace solver, close evaluation, Ewald sums
├── solvers.jl            # surfer geometry, surfactant and Marangoni solves
├── simulation.jl         # time integration
├── farfield.jl           # far-field (multipole) model
├── plotting.jl           # concentration on a grid; plotting stubs
└── script_params.jl      # name=value command-line overrides for the scripts
ext/                      # Plots.jl extension
test/                     # convergence tests, with their own environment
scripts/                  # data generators and the single/ example, with their own environment
figs/                     # MATLAB figure scripts and *_snapshot.jl data generators
```

## License

MIT; see [LICENSE](LICENSE).
