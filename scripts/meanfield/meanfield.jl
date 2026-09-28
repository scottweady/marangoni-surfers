
# Pseudo-spectral solver for the mean-field (Vlasov) limit of the far-field model.
#
# The distribution ψ(x, θ, t) of particle positions and orientations satisfies
#
#     ∂ψ/∂t + ∇·[(U₀ p + u) ψ] + ∂/∂θ[θ̇ ψ] = d_T Δψ + d_R ∂²ψ/∂θ²,
#
#     p(θ) = (cos θ, sin θ),   θ̇ = λ p⊥·E·p,
#
# closed by the polarization field m(x) = ∫ p ψ dθ through
#
#     u = ∇B[σ_eff],   σ_eff = β ∇·m,   B[f](x) = (1/4π)∫ |x-y| f(y) dy,
#
# and E = ∇u.  Since u is a gradient the flow is irrotational, so Jeffery's
# equation retains only its strain-rate contribution.
#
# All three coordinates are periodic, so everything is Fourier. Writing
# S(k) = β(k·m̂)/(2k³), the operators are diagonal:
#
#     û_i = k_i S,      Ê_ij = i k_j û_i,
#
# with the k = 0 mode set to zero (no mean flow). Note E carries no 1/k factor
# even though u does. In terms of E the rotation rate is
#
#     θ̇ = λ[(E₂₂ - E₁₁) sin(2θ)/2 + E₁₂ cos(2θ)],
#
# which needs only three scalar fields.
#
# Caveat: the interaction kernel has symbol ~1/k and therefore no intrinsic
# length scale, so results depend on the box size Lx, Ly — the lowest nonzero
# wavenumber is the most strongly amplified. Convergence in box size must be
# checked separately from convergence in resolution.

using FFTW
using Printf
using Random

"""
    MeanFieldGrid(Nx, Ny, Nθ; Lx=2π, Ly=2π)

Periodic collocation grid on `[0,Lx) × [0,Ly) × [0,2π)` for the mean-field solver.

Fields `ψ` are stored as real `Nx × Ny × Nθ` arrays, with the first index running over
`x`, the second over `y`, and the third over `θ`.

# Arguments
- `Nx`, `Ny`: number of grid points in `x` and `y`.
- `Nθ`: number of orientation grid points.
- `Lx`, `Ly`: box dimensions (default `2π`).

# Returns
- A `MeanFieldGrid` holding the collocation points, wavenumbers, and the 2/3-rule
  dealiasing masks.
"""
struct MeanFieldGrid
    Nx::Int
    Ny::Int
    Nθ::Int
    Lx::Float64
    Ly::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    θ::Vector{Float64}
    kx::Vector{Float64}
    ky::Vector{Float64}
    kθ::Vector{Float64}
    dx::Float64
    dy::Float64
    dθ::Float64
    keep_x::Vector{Bool}
    keep_y::Vector{Bool}
    keep_θ::Vector{Bool}
end

function MeanFieldGrid(Nx::Int, Ny::Int, Nθ::Int; Lx=2π, Ly=2π)

    x = collect((0:Nx-1) .* (Lx / Nx))
    y = collect((0:Ny-1) .* (Ly / Ny))
    θ = collect((0:Nθ-1) .* (2π / Nθ))

    nx = collect(fftfreq(Nx, Nx))       # integer mode numbers, FFT order
    ny = collect(fftfreq(Ny, Ny))
    nθ = collect(fftfreq(Nθ, Nθ))

    kx = (2π / Lx) .* nx
    ky = (2π / Ly) .* ny

    # 2/3 rule: the quadratic fluxes are truncated above N/3
    keep_x = abs.(nx) .<= Nx ÷ 3
    keep_y = abs.(ny) .<= Ny ÷ 3
    keep_θ = abs.(nθ) .<= Nθ ÷ 3

    return MeanFieldGrid(Nx, Ny, Nθ, Lx, Ly, x, y, θ, kx, ky, nθ,
                         Lx/Nx, Ly/Ny, 2π/Nθ, keep_x, keep_y, keep_θ)
end

"""
    meanfield_moments(ψ, G)

Zeroth and first orientational moments of `ψ`.

# Arguments
- `ψ`: distribution, an `Nx × Ny × Nθ` array.
- `G::MeanFieldGrid`: the grid.

# Returns
- `(ρ, m1, m2)`: concentration `ρ = ∫ψ dθ` and polarization `m = ∫pψ dθ`, each
  `Nx × Ny`.
"""
function meanfield_moments(ψ, G::MeanFieldGrid)
    cθ = reshape(cos.(G.θ), 1, 1, :)
    sθ = reshape(sin.(G.θ), 1, 1, :)
    ρ  = G.dθ .* dropdims(sum(ψ,        dims=3), dims=3)
    m1 = G.dθ .* dropdims(sum(ψ .* cθ,  dims=3), dims=3)
    m2 = G.dθ .* dropdims(sum(ψ .* sθ,  dims=3), dims=3)
    return ρ, m1, m2
end

"""
    meanfield_velocity(ψ, G; β=π/2)

Mean-field surface velocity and strain rate induced by `ψ`.

Applies the diagonal symbols `û_i = k_i β(k·m̂)/(2k³)` and `Ê_ij = i k_j û_i`, with the
`k = 0` mode set to zero.

# Arguments
- `ψ`: distribution, an `Nx × Ny × Nθ` array.
- `G::MeanFieldGrid`: the grid.
- `β`: surfactant dipole strength (default `π/2`, the value for a unit disk with
  `g(s) = -cos s`).

# Returns
- `(u1, u2, E11, E12, E22)`: velocity components and the independent components of the
  symmetric strain rate `E = ∇u`, each `Nx × Ny`.
"""
function meanfield_velocity(ψ, G::MeanFieldGrid; β=π/2)

    _, m1, m2 = meanfield_moments(ψ, G)
    m̂1, m̂2 = fft(m1), fft(m2)

    û1  = zeros(ComplexF64, G.Nx, G.Ny)
    û2  = zeros(ComplexF64, G.Nx, G.Ny)
    Ê11 = zeros(ComplexF64, G.Nx, G.Ny)
    Ê12 = zeros(ComplexF64, G.Nx, G.Ny)
    Ê22 = zeros(ComplexF64, G.Nx, G.Ny)

    @inbounds for j in 1:G.Ny, i in 1:G.Nx
        kx, ky = G.kx[i], G.ky[j]
        k2 = kx^2 + ky^2
        k2 == 0 && continue
        S = β * (kx * m̂1[i,j] + ky * m̂2[i,j]) / (2 * k2 * sqrt(k2))
        û1[i,j] = kx * S
        û2[i,j] = ky * S
        Ê11[i,j] = im * kx * û1[i,j]
        Ê12[i,j] = im * ky * û1[i,j]
        Ê22[i,j] = im * ky * û2[i,j]
    end

    return real.(ifft(û1)),  real.(ifft(û2)),
           real.(ifft(Ê11)), real.(ifft(Ê12)), real.(ifft(Ê22))

end

"""
    meanfield_rhs(ψ, G; U0=0.098175, λ=0.0, β=π/2, dT=0.0, dR=0.0, dealias=true)

Right-hand side `∂ψ/∂t` of the mean-field transport equation.

Fluxes are formed in physical space and differentiated spectrally, so the scheme is
conservative: `∫ψ` is preserved to roundoff. Diffusion is applied as an exact Fourier
multiplier and is not dealiased.

# Arguments
- `ψ`: distribution, an `Nx × Ny × Nθ` array.
- `G::MeanFieldGrid`: the grid.
- `U0`: single-particle propulsion speed (default `0.098175`, the unit disk).
- `λ`: Jeffery coefficient `(a²-b²)/(a²+b²)`; `0` for disks (default `0.0`).
- `β`: surfactant dipole strength (default `π/2`).
- `dT`, `dR`: translational and rotational diffusivities (default `0.0`).
- `dealias`: apply the 2/3-rule truncation to the nonlinear fluxes (default `true`).

# Returns
- `dψ`: the time derivative, an `Nx × Ny × Nθ` real array.
"""
function meanfield_rhs(ψ, G::MeanFieldGrid; U0=0.098175, λ=0.0, β=π/2,
                       dT=0.0, dR=0.0, dealias=true)

    u1, u2, E11, E12, E22 = meanfield_velocity(ψ, G; β=β)

    cθ  = reshape(cos.(G.θ),      1, 1, :)
    sθ  = reshape(sin.(G.θ),      1, 1, :)
    c2θ = reshape(cos.(2 .* G.θ), 1, 1, :)
    s2θ = reshape(sin.(2 .* G.θ), 1, 1, :)

    u1r = reshape(u1, G.Nx, G.Ny, 1)
    u2r = reshape(u2, G.Nx, G.Ny, 1)

    # conservative fluxes
    Fx = (U0 .* cθ .+ u1r) .* ψ
    Fy = (U0 .* sθ .+ u2r) .* ψ

    θdot = λ .* (0.5 .* reshape(E22 .- E11, G.Nx, G.Ny, 1) .* s2θ .+
                        reshape(E12,        G.Nx, G.Ny, 1) .* c2θ)
    Fθ = θdot .* ψ

    # spatial divergence (+ translational diffusion)
    F̂x, F̂y = fft(Fx, (1,2)), fft(Fy, (1,2))
    ψ̂xy = dT > 0 ? fft(ψ, (1,2)) : nothing
    A = similar(F̂x)
    @inbounds for m in 1:G.Nθ, j in 1:G.Ny, i in 1:G.Nx
        kx, ky = G.kx[i], G.ky[j]
        nl = -im * (kx * F̂x[i,j,m] + ky * F̂y[i,j,m])
        if dealias && !(G.keep_x[i] && G.keep_y[j])
            nl = 0.0 + 0.0im
        end
        A[i,j,m] = dT > 0 ? nl - dT * (kx^2 + ky^2) * ψ̂xy[i,j,m] : nl
    end
    dψ = real.(ifft(A, (1,2)))

    # orientational divergence (+ rotational diffusion)
    F̂θ = fft(Fθ, 3)
    ψ̂θ = dR > 0 ? fft(ψ, 3) : nothing
    B = similar(F̂θ)
    @inbounds for m in 1:G.Nθ, j in 1:G.Ny, i in 1:G.Nx
        n = G.kθ[m]
        nl = -im * n * F̂θ[i,j,m]
        if dealias && !G.keep_θ[m]
            nl = 0.0 + 0.0im
        end
        B[i,j,m] = dR > 0 ? nl - dR * n^2 * ψ̂θ[i,j,m] : nl
    end
    dψ .+= real.(ifft(B, 3))

    return dψ

end

"""
    meanfield_timestep(ψ, G; U0, λ, β, cfl=0.2, dt_max=Inf)

Time step from the advective CFL condition in `(x, y, θ)`. Diffusion is treated
implicitly and imposes no restriction.

# Returns
- `(dt, umax, θdotmax)`: the step and the maximum advection speeds, for reporting.
"""
function meanfield_timestep(ψ, G::MeanFieldGrid; U0=0.098175, λ=0.0, β=π/2, cfl=0.2, dt_max=Inf)

    u1, u2, E11, E12, E22 = meanfield_velocity(ψ, G; β=β)

    ux = U0 + maximum(abs, u1)
    uy = U0 + maximum(abs, u2)
    # |θ̇| ≤ λ sqrt(((E22-E11)/2)² + E12²), maximized over θ
    θdot = λ == 0 ? 0.0 :
           abs(λ) * maximum(sqrt.(((E22 .- E11) ./ 2).^2 .+ E12.^2))

    rate = ux / G.dx + uy / G.dy + θdot / G.dθ
    dt = rate > 0 ? cfl / rate : dt_max
    return min(dt, dt_max), max(ux, uy) - U0, θdot

end

"""
    meanfield_simulate(ψ0, G; t_end=10.0, ...)

Integrate the mean-field equation with the variable-step IMEX scheme SBDF2.

Advection in `(x, y, θ)` is extrapolated explicitly and diffusion `dT Δ + dR ∂²_θ`,
diagonal in Fourier space, is treated implicitly, so the step is limited only by the
advective CFL condition. With `ω = dtₙ/dtₙ₋₁` the update is
`a₀ψⁿ⁺¹ + a₁ψⁿ + a₂ψⁿ⁻¹ = dt[(1+ω)Nⁿ - ωNⁿ⁻¹] + dt Dψⁿ⁺¹` with
`a₀ = (1+2ω)/(1+ω)`, `a₁ = -(1+ω)`, `a₂ = ω²/(1+ω)`; the first step takes `ω = 0`, which
is IMEX Euler. The step ratio is capped at `ω_max` (zero-stability requires `ω < 1+√2`).

The explicit part is weakly unstable on the imaginary axis (`|ζ| - 1 ≈ 0.75(dt λ)⁴` per
step, as for AB2), so undamped spectral advection needs a smaller CFL factor than
SSPRK3; diffusion damps the high modes where this matters.

Mass `∫ψ dx dy dθ` is conserved to roundoff; `min ψ` is reported so that loss of
positivity — which signals under-resolution — is visible.

# Arguments
- `ψ0`: initial distribution, an `Nx × Ny × Nθ` array.
- `G::MeanFieldGrid`: the grid.
- `t_end`: final time (default `10.0`).
- `U0`, `λ`, `β`, `dT`, `dR`: model parameters, as in [`meanfield_rhs`](@ref).
- `cfl`: CFL safety factor (default `0.2`).
- `dt_max`: cap on the time step (default `Inf`).
- `ω_max`: cap on the step ratio `dtₙ/dtₙ₋₁` (default `1.25`).
- `tplt`: time between recorded snapshots (default `t_end/20`).
- `output_dir`: if non-empty, write `frame_<n>.dat` files of `(x, y, ρ, m1, m2, u1, u2)`
  here (default `""`).
- `callback`: if given, called as `callback(snap, G)` on every recorded snapshot,
  including the initial one. Used for live plotting (default `nothing`).
- `verbose`: print per-step diagnostics (default `true`).

# Returns
- `(ψ, snapshots)` where `snapshots` is a vector of named tuples
  `(t, ρ, m1, m2, u1, u2, mass, minψ)`.
"""
function meanfield_simulate(ψ0, G::MeanFieldGrid; t_end=10.0,
                            U0=0.098175, λ=0.0, β=π/2, dT=0.0, dR=0.0,
                            cfl=0.2, dt_max=Inf, ω_max=1.25, tplt=t_end/20,
                            output_dir="", callback=nothing,
                            verbose=true, dealias=true)

    ψ = copy(ψ0)
    cell = G.dx * G.dy * G.dθ

    # explicit advection and implicit diffusion symbol D(k, n) = -(dT k² + dR n²)
    Nexp(f) = meanfield_rhs(f, G; U0=U0, λ=λ, β=β, dT=0.0, dR=0.0, dealias=dealias)
    D = [-(dT * (G.kx[i]^2 + G.ky[j]^2) + dR * G.kθ[m]^2)
         for i in 1:G.Nx, j in 1:G.Ny, m in 1:G.Nθ]

    ψ_prev = similar(ψ)
    N_prev = similar(ψ)
    dt_prev = 0.0

    snapshots = NamedTuple[]

    function record!(t)
        ρ, m1, m2 = meanfield_moments(ψ, G)
        u1, u2, _, _, _ = meanfield_velocity(ψ, G; β=β)
        snap = (t=t, ρ=copy(ρ), m1=copy(m1), m2=copy(m2), u1=copy(u1), u2=copy(u2),
                mass=sum(ψ)*cell, minψ=minimum(ψ))
        push!(snapshots, snap)
        if !isempty(output_dir)
            mkpath(output_dir)
            path = joinpath(output_dir, "frame_$(lpad(length(snapshots)-1, 6, '0')).dat")
            open(path, "w") do io
                @printf(io, "# t = %.8f   Nx = %d   Ny = %d\n", t, G.Nx, G.Ny)
                println(io, "# x y rho m1 m2 u1 u2")
                for j in 1:G.Ny, i in 1:G.Nx
                    @printf(io, "%.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
                            G.x[i], G.y[j], ρ[i,j], m1[i,j], m2[i,j], u1[i,j], u2[i,j])
                end
            end
        end
        callback === nothing || callback(snap, G)
        return snap
    end

    t = 0.0
    next_save = tplt
    s0 = record!(t)
    mass0 = s0.mass
    verbose && @printf("t = %8.4f  mass = %.12e  min ψ = %+.3e\n", t, s0.mass, s0.minψ)

    nstep = 0
    while t < t_end

        dt, umax, θdotmax = meanfield_timestep(ψ, G; U0=U0, λ=λ, β=β,
                                               cfl=cfl, dt_max=dt_max)
        nstep > 0 && (dt = min(dt, ω_max * dt_prev))

        # land on the next save time without leaving a tiny final step, which
        # would make the following step ratio large
        t_stop = min(t_end, next_save)
        if t + dt >= t_stop - 1e-12
            dt = t_stop - t
        elseif t + 2dt > t_stop
            dt = (t_stop - t) / 2
        end

        # SBDF2 (variable step); ω = 0 on the first step gives IMEX Euler
        ω  = nstep == 0 ? 0.0 : dt / dt_prev
        a0 = (1 + 2ω) / (1 + ω)
        a1 = -(1 + ω)
        a2 = ω^2 / (1 + ω)

        Nn  = Nexp(ψ)
        rhs = nstep == 0 ? ψ .+ dt .* Nn :
              -a1 .* ψ .- a2 .* ψ_prev .+ dt .* ((1 + ω) .* Nn .- ω .* N_prev)
        ψ̂ = fft(rhs) ./ (a0 .- dt .* D)

        ψ_prev .= ψ
        N_prev .= Nn
        dt_prev = dt
        ψ = real.(ifft(ψ̂))

        t += dt
        nstep += 1

        if t >= next_save - 1e-12
            s = record!(t)
            verbose && @printf("t = %8.4f  dt = %.2e  |u|max = %.3e  |θ̇|max = %.3e  Δmass/mass = %+.2e  min ψ = %+.3e\n",
                               t, dt, umax, θdotmax, (s.mass - mass0)/mass0, s.minψ)
            next_save += tplt
        end

    end

    verbose && @printf("done: %d steps\n", nstep)
    return ψ, snapshots

end

# ── Initial conditions ────────────────────────────────────────────────────────

"""
    ic_vonmises(G, ρ, θ0; κ=4.0)

Build `ψ(x,θ) = ρ(x) f(θ - θ0(x))` with `f` the von Mises density of concentration `κ`,
normalized so that `∫ψ dθ = ρ(x)`. Positive by construction.

# Arguments
- `G::MeanFieldGrid`: the grid.
- `ρ`: function `(x,y) -> ρ` giving the local concentration.
- `θ0`: function `(x,y) -> θ` giving the local preferred orientation.
- `κ`: orientational concentration; `κ = 0` is isotropic, large `κ` is strongly aligned
  (default `4.0`).

# Returns
- `ψ`: an `Nx × Ny × Nθ` array.
"""
function ic_vonmises(G::MeanFieldGrid, ρ, θ0; κ=4.0)
    ψ = zeros(G.Nx, G.Ny, G.Nθ)
    for j in 1:G.Ny, i in 1:G.Nx
        w = exp.(κ .* cos.(G.θ .- θ0(G.x[i], G.y[j])))
        w ./= sum(w) * G.dθ                      # ∫ w dθ = 1
        @views ψ[i,j,:] .= ρ(G.x[i], G.y[j]) .* w
    end
    return ψ
end

"""
    ic_isotropic(G; ρ0=1.0)

Uniform isotropic state `ψ = ρ₀/2π`. An exact steady state of the mean-field equation.
"""
ic_isotropic(G::MeanFieldGrid; ρ0=1.0) = fill(ρ0 / (2π), G.Nx, G.Ny, G.Nθ)

"""
    ic_polar(G; ρ0=1.0, θ0=0.0, κ=4.0, ε=0.0, nx=1, ny=0)

Uniform polar state, optionally seeded with a density perturbation of relative amplitude
`ε` at integer wavenumbers `(nx, ny)`.

The unperturbed state has `∇·m = 0` and hence `u ≡ 0`, so it is an exact steady state;
`ε > 0` seeds a perturbation for linear stability measurements.
"""
function ic_polar(G::MeanFieldGrid; ρ0=1.0, θ0=0.0, κ=4.0, ε=0.0, nx=1, ny=0)
    kx = 2π * nx / G.Lx
    ky = 2π * ny / G.Ly
    ρfun(x, y) = ρ0 * (1 + ε * cos(kx * x + ky * y))
    return ic_vonmises(G, ρfun, (x, y) -> θ0; κ=κ)
end

"""
    ic_isotropic_perturbed(G; ρ0=1.0, modes)

Uniform isotropic state seeded with a sum of plane waves in the concentration and/or the
polarization,

    ψ = (ρ₀/2π) [1 + Σⱼ ερⱼ cos(kⱼ·x + φⱼ) + 2 Σⱼ εmⱼ cos(θ - θmⱼ) cos(kⱼ·x + φⱼ)],

with `kⱼ = 2π(nxⱼ/Lx, nyⱼ/Ly)`. The moments are then exactly

    ρ = ρ₀[1 + Σⱼ ερⱼ cos(kⱼ·x + φⱼ)],
    m = ρ₀ Σⱼ εmⱼ cos(kⱼ·x + φⱼ) (cos θmⱼ, sin θmⱼ),

so concentration and polarization can be dialled independently mode by mode. With no
modes this is the uniform isotropic state, an exact steady state.

Since `u` responds only to `∇·m`, each mode drives a flow proportional to
`kⱼ·(cos θmⱼ, sin θmⱼ)`: `θmⱼ` along `kⱼ` is the longitudinal perturbation that drives
the strongest flow, while `θmⱼ` perpendicular to `kⱼ` drives none at `t = 0`. A pure
concentration perturbation (`εm = 0`) still evolves, because self-propulsion advects each
orientation class in its own direction and so generates polarization from `∇ρ`.

# Arguments
- `G::MeanFieldGrid`: the grid.
- `ρ0`: background concentration (default `1.0`).
- `modes`: iterable of named tuples, one per plane wave. `nx` and `ny` (integer
  wavenumbers, not both zero) are required; `ερ` (concentration amplitude, default `0`),
  `εm` (polarization amplitude, default `0`), `phase` (default `0`), and `θm`
  (polarization direction, default: along `k`, i.e. longitudinal) are optional. Build a
  random set with [`random_modes`](@ref).

# Returns
- `ψ`: an `Nx × Ny × Nθ` array. Positivity requires `Σⱼ (ερⱼ + 2|εmⱼ|) < 1`.
"""
function ic_isotropic_perturbed(G::MeanFieldGrid; ρ0=1.0, modes=())

    budget = sum(mode -> get(mode, :ερ, 0.0) + 2abs(get(mode, :εm, 0.0)), modes; init=0.0)
    budget < 1 || error("perturbation too large: need Σ(ερ + 2|εm|) < 1 for ψ ≥ 0 " *
                        "(got $budget over $(length(modes)) modes)")

    ψ = fill(ρ0 / (2π), G.Nx, G.Ny, G.Nθ)

    for mode in modes
        nx, ny = mode.nx, mode.ny
        (nx == 0 && ny == 0) && error("mode (0,0) is a uniform shift, not a perturbation")

        ερ    = get(mode, :ερ, 0.0)
        εm    = get(mode, :εm, 0.0)
        phase = get(mode, :phase, 0.0)

        kx = 2π * nx / G.Lx
        ky = 2π * ny / G.Ly

        θmj = get(mode, :θm, nothing)
        θm  = θmj === nothing ? atan(ky, kx) : θmj    # default: longitudinal

        for m in 1:G.Nθ, j in 1:G.Ny, i in 1:G.Nx
            c = cos(kx * G.x[i] + ky * G.y[j] + phase)
            ψ[i,j,m] += (ρ0 / (2π)) * (ερ * c + 2 * εm * cos(G.θ[m] - θm) * c)
        end
    end

    return ψ

end

"""
    random_modes(nmodes; nmin=1, nmax=4, ερ=0.0, εm=1e-3, θm=nothing, seed=0)

Build a set of `nmodes` distinct plane waves with random integer wavevectors in the
band `nmin ≤ |n| ≤ nmax` and random phases, for use as the `modes` argument of
[`ic_isotropic_perturbed`](@ref).

Modes are drawn without replacement, and `(nx, ny)` and `(-nx, -ny)` are treated as the
same wave (the perturbation is real), so at most one of each pair appears. If the band
holds fewer than `nmodes` distinct waves, all of them are returned.

# Arguments
- `nmodes`: number of plane waves requested.
- `nmin`, `nmax`: inclusive bounds on `|n| = hypot(nx, ny)` (default `1`, `4`).
- `ερ`, `εm`: per-mode concentration and polarization amplitudes (default `0.0`, `1e-3`).
  Positivity requires `nmodes * (ερ + 2|εm|) < 1`.
- `θm`: polarization direction for every mode; `nothing` leaves each mode longitudinal,
  i.e. aligned with its own `k` (default `nothing`).
- `seed`: RNG seed, so a run is reproducible (default `0`).

# Returns
- A vector of named tuples suitable for `ic_isotropic_perturbed(G; modes=...)`.
"""
function random_modes(nmodes::Int; nmin=1, nmax=4, ερ=0.0, εm=1e-3, θm=nothing, seed=0)

    rng = MersenneTwister(seed)

    # enumerate the band, keeping one of each ±(nx,ny) pair
    candidates = NTuple{2,Int}[]
    for ny in -nmax:nmax, nx in -nmax:nmax
        (nx == 0 && ny == 0) && continue
        n = hypot(nx, ny)
        (nmin - 1e-9 <= n <= nmax + 1e-9) || continue
        (-nx, -ny) in candidates && continue
        push!(candidates, (nx, ny))
    end

    nkeep = min(nmodes, length(candidates))
    nkeep < nmodes && @warn "band [$nmin, $nmax] holds only $nkeep distinct waves; using all of them"

    picked = candidates[randperm(rng, length(candidates))[1:nkeep]]

    return [(nx = nx, ny = ny, ερ = ερ, εm = εm, θm = θm,
             phase = 2π * rand(rng)) for (nx, ny) in picked]

end

"""
    ic_patch(G; ρ0=1.0, R=1.0, center=nothing, width=0.25, κ=4.0)

Radially-oriented circular patch of surfers of radius `R`, with a `tanh` edge of width
`width`. Reproduces the initial condition of the patch-expansion simulations.
"""
function ic_patch(G::MeanFieldGrid; ρ0=1.0, R=1.0, center=nothing, width=0.25, κ=1.0)
    xc, yc = something(center, (G.Lx/2, G.Ly/2))
    ρfun(x, y) = ρ0 * 0.5 * (1 - tanh((hypot(x - xc, y - yc) - R) / width))
    θfun(x, y) = atan(y - yc, x - xc)
    return ic_vonmises(G, ρfun, θfun; κ=κ)
end

# ── Diagnostics ───────────────────────────────────────────────────────────────

"""
    radial_profile(f, G; center=nothing, nbins=64, Rmax=nothing)

Azimuthally averaged radial profile of a scalar field `f` on the grid.

# Returns
- `(r, f̄, counts)`: bin centers, bin means, and the number of grid points per bin.
"""
function radial_profile(f, G::MeanFieldGrid; center=nothing, nbins=64, Rmax=nothing)
    xc, yc = something(center, (G.Lx/2, G.Ly/2))
    Rmax = something(Rmax, min(G.Lx, G.Ly)/2)
    edges = range(0, Rmax, length=nbins+1)
    sums = zeros(nbins); counts = zeros(Int, nbins)
    for j in 1:G.Ny, i in 1:G.Nx
        r = hypot(G.x[i] - xc, G.y[j] - yc)
        r >= Rmax && continue
        b = max(1, min(nbins, Int(floor(r / (Rmax/nbins))) + 1))
        sums[b] += f[i,j]; counts[b] += 1
    end
    means = [counts[b] > 0 ? sums[b]/counts[b] : NaN for b in 1:nbins]
    return collect(edges[1:end-1] .+ step(edges)/2), means, counts
end

"""
    mode_amplitude(f, G, nx, ny)

Amplitude `|f̂|` of the integer Fourier mode `(nx, ny)` of a scalar field, normalized by
the number of grid points. Used to measure growth rates in stability studies.
"""
function mode_amplitude(f, G::MeanFieldGrid, nx::Int, ny::Int)
    f̂ = fft(f)
    i = mod(nx, G.Nx) + 1
    j = mod(ny, G.Ny) + 1
    return abs(f̂[i,j]) / (G.Nx * G.Ny)
end
