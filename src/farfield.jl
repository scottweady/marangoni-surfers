using LinearAlgebra
using Printf
using SpecialFunctions

"""
    farfield_constants(aspect_ratio; s₀=π/2)

Exact single-particle constants entering the far-field model, for an elliptical surfer of
unit area (`ab = 1`) with boundary concentration `g(s) = -cos(s)` on the Dirichlet arc
`s ∈ (s₀, 2π - s₀)`.

The propulsion speed follows from the disk solution of Lauga & Davis mapped to the ellipse
through elliptic coordinates; for `a < b` the roles of the principal axes are swapped. The
dipole strength `η` is defined by `∫ (y - x_cm) σ dΓ = -η p`.

# Arguments
- `aspect_ratio`: `a/b`, with `a` the half-axis along the propulsion direction.
- `s₀`: Dirichlet-Neumann junction parameter (default `π/2`).

# Returns
- `(U, η, β)`: propulsion speed, dipole strength, and Bretherton constant
  `(a² - b²)/(a² + b²)`.
"""
function farfield_constants(aspect_ratio; s₀=π/2)
    a = sqrt(aspect_ratio)
    b = 1 / a
    lo = min(a, b)
    m  = 1 - (lo / max(a, b))^2

    # k_i - k_0 along the propulsion axis, with K and E the complete elliptic integrals
    dk = if a == b
        -π / 8 * lo
    elseif a > b
        lo / (2m) * (ellipk(m) - ellipe(m)) - lo / 2 * ellipk(m)
    else
        lo / (2m) * (ellipe(m) - (1 - m) * ellipk(m)) - lo / 2 * ellipk(m)
    end

    U = -(a + b) * dk * (1 + cos(s₀))^2 / (8 * a * b)
    η = π * (a + b) * (1 + cos(s₀))^2 / 4
    β = (a^2 - b^2) / (a^2 + b^2)
    return U, η, β
end

"""
    farfield_response(aspect_ratio)

Moments `(τ₀, τ₁)` of the response of an isolated surfer to a uniform shift in ambient
concentration, `∫ τ dΓ = τ₀` and `∫ (y - x_cm) τ dΓ = τ₁ p`, for the unit-area ellipse with
the Dirichlet arc `s ∈ (π/2, 3π/2)`.

For the unit disk `τ₀ = 2π/log 2`, so that the effective capacity radius
`ρₑ = exp(-2π/τ₀)` is `1/2`, and `τ₁ = -τ₀`. The ellipse follows from the disk through
elliptic coordinates, which rescale lengths at infinity by `(a + b)/2`, giving
`τ₀ = 2π/log(4/(a + b))` and `τ₁ = -τ₀ (a + b)/2`. These agree with
[`surfactant_response`](@ref) to machine precision.

# Arguments
- `aspect_ratio`: `a/b`, with `a` the half-axis along the propulsion direction.

# Returns
- `(τ₀, τ₁)`: the response moments.
"""
function farfield_response(aspect_ratio)
    a = sqrt(aspect_ratio)
    b = 1 / a
    τ₀ = 2π / log(4 / (a + b))
    return τ₀, -τ₀ * (a + b) / 2
end

"""
    farfield_moments(xcm, θcm, η, τ₀, τ₁)

Surfactant monopoles and dipole strengths of `N` surfers from the monopole closure.

Each particle responds to the concentration produced at its center by all others, through
the amplitude `αᵢ` of the response `τ`: `m₀ᵢ = τ₀ αᵢ` and `ηᵢ = η - τ₁ αᵢ`. The amplitudes and
the shift in the far-field constant `c∞′` solve the `(N+1) × (N+1)` system

    αᵢ - (τ₀/2π) Σⱼ αⱼ log rᵢⱼ + (τ₁/2π) Σⱼ αⱼ (r̂ᵢⱼ·pⱼ)/rᵢⱼ + c∞′ = (η/2π) Σⱼ (r̂ᵢⱼ·pⱼ)/rᵢⱼ,
    Σᵢ αᵢ = 0,

with `rᵢⱼ = xᵢ - xⱼ`, which is solved directly.

# Arguments
- `xcm`: complex particle centers.
- `θcm`: particle orientations, in radians.
- `η`: single-particle dipole strength from [`farfield_constants`](@ref).
- `τ₀`, `τ₁`: response moments from [`farfield_response`](@ref).

# Returns
- `(m₀, ηᵢ)`: per-particle monopoles and dipole strengths.
"""
function farfield_moments(xcm, θcm, η, τ₀, τ₁; cache=nothing)
    α, = farfield_closure(xcm, θcm, η, τ₀, τ₁; cache=cache)
    return τ₀ .* α, η .- τ₁ .* α
end

"""
    farfield_closure(xcm, θcm, η, τ₀, τ₁; cache=nothing, rtol=1e-12, maxiter=10)

Solve the monopole closure system of [`farfield_moments`](@ref) for the response
amplitudes `αᵢ` and the shift `c∞′` of the far-field constant relative to an isolated
particle.

With `cache` a `Ref` holding an LU factorization from an earlier, nearby configuration (or
`nothing`), the system is solved by iterative refinement with that factorization to relative
residual `rtol`, and refactored, updating `cache`, if this takes more than `maxiter` steps.
Since the configuration changes little between evaluations, this avoids an `O(N³)`
factorization at almost every step.

# Returns
- `(α, c∞′)`: per-particle amplitudes and the far-field shift.
"""
function farfield_closure(xcm, θcm, η, τ₀, τ₁; cache=nothing, rtol=1e-12, maxiter=10)
    N = length(xcm)
    p = cis.(θcm)
    A = zeros(N + 1, N + 1)
    b = zeros(N + 1)
    Threads.@threads for i in 1:N
        A[i, i] = 1.0
        A[i, N+1] = 1.0
        A[N+1, i] = 1.0
        for j in 1:N
            j == i && continue
            r  = xcm[i] - xcm[j]
            d  = abs(r)
            rp = real(conj(r) * p[j]) / d^2     # (r̂ᵢⱼ · pⱼ) / rᵢⱼ
            A[i, j] = -τ₀ / 2π * log(d) + τ₁ / 2π * rp
            b[i]   += η / 2π * rp
        end
    end
    cache === nothing && (s = A \ b; return s[1:N], s[N+1])

    # iterative refinement with a factorization from an earlier configuration, refactoring
    # when it no longer converges
    if cache[] !== nothing
        s = cache[] \ b
        for _ in 1:maxiter
            r = b - A * s
            norm(r) <= rtol * norm(b) && return s[1:N], s[N+1]
            s .+= cache[] \ r
        end
    end
    cache[] = lu(A)
    s = cache[] \ b
    return s[1:N], s[N+1]
end

"""
    farfield_velocities(xcm, θcm, U, η, β, m₀=zeros(length(xcm)), ηᵢ=fill(η, length(xcm)))

Translational and angular velocities of `N` surfers in the far-field model.

Each particle self-propels along `p = (cos θ, sin θ)` with speed `U ηᵢ/η`, since the
self-induced velocity of an ellipse depends only on its own dipole, and is advected and
rotated (Jeffery, with Bretherton constant `β`) by the monopolar and dipolar Marangoni flows
`m₀ⱼ r̂/(4π) + ηⱼ (I - r̂r̂)·pⱼ / (4π r)` of all other particles, with `r = xᵢ - xⱼ`. With the
default moments this is the dipole-only model; with moments from
[`farfield_moments`](@ref) it is the monopole closure. The pairwise sum is evaluated
directly at `O(N²)` cost.

# Arguments
- `xcm`: complex particle centers.
- `θcm`: particle orientations, in radians.
- `U`, `η`, `β`: constants from [`farfield_constants`](@ref).
- `m₀`: per-particle monopoles (default zero).
- `ηᵢ`: per-particle dipole strengths (default `η`).

# Returns
- `(u, ω)`: complex translational velocities and angular velocities.
"""
function farfield_velocities(xcm, θcm, U, η, β, m₀=zeros(length(xcm)), ηᵢ=fill(η, length(xcm)))
    N = length(xcm)
    p = cis.(θcm)
    u = (U / η) .* ηᵢ .* p
    ω = zeros(N)

    Threads.@threads for i in 1:N
        ui, ωi = zero(ComplexF64), 0.0
        q = im * p[i]
        for j in 1:N
            j == i && continue
            r  = xcm[i] - xcm[j]
            d  = abs(r)
            e  = r / d
            c  = ηᵢ[j] / 4π
            ep  = real(conj(e) * p[j])     # r̂ · p_j
            qe  = real(conj(q) * e)
            epi = real(conj(e) * p[i])
            ui += m₀[j] / 4π * e + c * (p[j] - e * ep) / d
            ωi += -β * m₀[j] / 4π * qe * epi / d -
                  c * β * (-3 * qe * epi * ep + qe * real(conj(p[j]) * p[i]) +
                           real(conj(q) * p[j]) * epi) / d^2
        end
        u[i] += ui
        ω[i]  = ωi
    end
    return u, ω
end

"""
    repulsion_velocities!(u, xcm, A, r₀)

Add to `u` the velocities from a soft pairwise repulsion between particle centers,
`A (1 - r/r₀)² r̂` for `r < r₀`, with `r = xᵢ - xⱼ`. The pairwise sum is evaluated directly at
`O(N²)` cost.
"""
function repulsion_velocities!(u, xcm, A, r₀)
    Threads.@threads for i in eachindex(xcm)
        ui = zero(ComplexF64)
        for j in eachindex(xcm)
            j == i && continue
            r = xcm[i] - xcm[j]
            d = abs(r)
            d < r₀ && (ui += A * (1 - d / r₀)^2 * r / d)
        end
        u[i] += ui
    end
    return u
end

"""
    pair_bounds(xcm, u, r_c)

Smallest center-to-center distance over all pairs of particles, and the largest relative
speed `|uᵢ - uⱼ|` over pairs closer than `r_c`.
"""
function pair_bounds(xcm, u, r_c)
    N = length(xcm)
    dmin, vmax = fill(Inf, N), zeros(N)
    Threads.@threads for i in 1:N
        for j in i+1:N
            d = abs(xcm[i] - xcm[j])
            dmin[i] = min(dmin[i], d)
            d < r_c && (vmax[i] = max(vmax[i], abs(u[i] - u[j])))
        end
    end
    return minimum(dmin), maximum(vmax)
end

"""
    ellipse_outline(aspect_ratio, n)

Body-frame boundary points of a unit-area ellipse with half-axis `a = √aspect_ratio` along
the real axis, at `n` equispaced parameter values.
"""
function ellipse_outline(aspect_ratio, n)
    a = sqrt(aspect_ratio)
    s = range(0, 2π; length=n + 1)[1:n]
    return @. a * cos(s) + im * sin(s) / a
end

"""
    min_gap(xcm, θcm, z_body)

Smallest boundary-to-boundary distance over all pairs of particles, measured between the
placed boundary points `z_body`. Pairs whose centers are farther apart than twice the
largest body radius plus the current minimum are skipped.

# Returns
- `(gap, i, j)`: the minimum gap and the pair attaining it.
"""
function min_gap(xcm, θcm, z_body)
    N = length(xcm)
    R = maximum(abs, z_body)
    gmin, imin, jmin = Inf, 0, 0
    for i in 1:N, j in i+1:N
        abs(xcm[i] - xcm[j]) - 2R < gmin || continue
        zi = cis(θcm[i]) .* z_body .+ xcm[i]
        zj = cis(θcm[j]) .* z_body .+ xcm[j]
        g  = sqrt(minimum(abs2.(zi .- transpose(zj))))
        g < gmin && ((gmin, imin, jmin) = (g, i, j))
    end
    return gmin, imin, jmin
end

"""
    farfield_simulate(xcm0, θcm0, aspect_ratio; t_end=100.0, dt_max=0.5, dt_min=1e-4,
                      safety=0.25, save_every=5.0, stop_gap=1e-2, outdir="",
                      n_outline=64, closure=false, repulsion=0.0, r_rep=2.0, U_factor=1.0,
                      verbose=true)

Integrate the far-field model with classical RK4 from the configuration `(xcm0, θcm0)`.

With `closure = true`, the monopoles and dipole strengths of all particles are recomputed
from [`farfield_moments`](@ref) at every stage; otherwise every particle carries the
single-particle dipole and no monopole.

The step is limited so no pair can close more than a fraction `safety` of the current
minimum gap within a step, using twice the largest particle speed as a bound on the
closing speed. Integration stops at `t_end` or on contact, when the minimum gap falls
below `stop_gap`.

With `repulsion > 0`, particles are treated as points with the soft repulsion of
[`repulsion_velocities!`](@ref) and integration does not stop on contact. The gap is then
the minimum center-to-center distance, and the closing speed is bounded by the largest
relative speed among pairs closer than `4 r_rep`, while pairs farther apart may close by at
most half that cutoff, at twice the largest speed relative to the mean. This leaves out the
collective motion that dominates the particle speeds in a large patch.

Configurations are saved every `save_every` time units, and at the final time, as
`frame_<n>.dat` files in the same format as [`suspension_simulate`](@ref): one row
`t x y θ ux uy ω` per particle.

# Arguments
- `xcm0`: initial complex particle centers.
- `θcm0`: initial orientations, in radians.
- `aspect_ratio`: `a/b` of the (unit-area) elliptical particles.
- `t_end`: final time (default `100.0`).
- `dt_max`, `dt_min`: bounds on the time step (default `0.5`, `1e-4`).
- `safety`: fraction of the minimum gap that may close in one step (default `0.25`).
- `save_every`: interval between saved frames (default `5.0`).
- `stop_gap`: contact threshold on the minimum gap (default `1e-2`).
- `outdir`: directory for frame files; nothing is written if empty (default `""`).
- `n_outline`: boundary points per particle used to measure gaps (default `64`).
- `closure`: use the monopole closure (default `false`).
- `repulsion`, `r_rep`: strength and range of the soft repulsion (default `0.0`, `2.0`).
- `U_factor`: factor multiplying the self-propulsion speed, leaving `η` and `β` unchanged
  (default `1.0`).
- `verbose`: print progress at each saved frame (default `true`).

# Returns
- `(xcm, θcm, t, contact)`: final configuration, final time, and whether the run stopped
  on contact.
"""
function farfield_simulate(xcm0, θcm0, aspect_ratio;
                           t_end=100.0, dt_max=1e-2, dt_min=1e-4, safety=0.25,
                           save_every=5.0, stop_gap=1e-2, outdir="", n_outline=64,
                           closure=false, repulsion=0.0, r_rep=2.0, U_factor=1.0,
                           verbose=true)

    U, η, β = farfield_constants(aspect_ratio)
    U *= U_factor
    τ₀, τ₁  = farfield_response(aspect_ratio)
    z_body  = ellipse_outline(aspect_ratio, n_outline)
    point   = repulsion > 0
    lu_cache = Ref{Any}(nothing)
    function rhs(x, θ)
        u, ω = closure ? farfield_velocities(x, θ, U, η, β, farfield_moments(x, θ, η, τ₀, τ₁; cache=lu_cache)...) :
                         farfield_velocities(x, θ, U, η, β)
        point && repulsion_velocities!(u, x, repulsion, r_rep)
        return u, ω
    end
    r_c     = 4 * r_rep

    # time step from the gap and the velocities at the start of the step
    function step_size(x, u, gap)
        point || return safety * gap / (2 * maximum(abs, u))
        d, v = pair_bounds(x, u, r_c)
        return min(safety * d / v, r_c / (4 * maximum(abs, u .- sum(u) / length(u))))
    end
    gap_of(x, θ) = point ? first(pair_bounds(x, zeros(ComplexF64, length(x)), 0.0)) :
                           first(min_gap(x, θ, z_body))

    xcm, θcm = ComplexF64.(xcm0), Float64.(θcm0)
    t, frame = 0.0, 0
    next_save = save_every
    isempty(outdir) || mkpath(outdir)

    function save(t, x, θ)
        isempty(outdir) && return
        u, ω = rhs(x, θ)
        open(joinpath(outdir, "frame_$(lpad(frame, 6, '0')).dat"), "w") do io
            for i in eachindex(x)
                @printf(io, "%.6f %.10f %.10f %.10f %.10f %.10f %.10f\n",
                        t, real(x[i]), imag(x[i]), θ[i], real(u[i]), imag(u[i]), ω[i])
            end
        end
        frame += 1
    end

    save(t, xcm, θcm)
    gap = gap_of(xcm, θcm)
    contact = false

    while t < t_end - 1e-12
        if !point && gap < stop_gap
            contact = true
            break
        end

        u1, ω1 = rhs(xcm, θcm)
        dt = clamp(step_size(xcm, u1, gap), dt_min, dt_max)
        dt = min(dt, t_end - t, next_save - t)

        u2, ω2 = rhs(xcm .+ dt/2 .* u1, θcm .+ dt/2 .* ω1)
        u3, ω3 = rhs(xcm .+ dt/2 .* u2, θcm .+ dt/2 .* ω2)
        u4, ω4 = rhs(xcm .+ dt   .* u3, θcm .+ dt   .* ω3)
        xcm .+= dt/6 .* (u1 .+ 2u2 .+ 2u3 .+ u4)
        θcm .+= dt/6 .* (ω1 .+ 2ω2 .+ 2ω3 .+ ω4)
        t   += dt

        point || (gap = gap_of(xcm, θcm))
        if t >= next_save - 1e-12
            point && (gap = gap_of(xcm, θcm))
            save(t, xcm, θcm)
            next_save += save_every
            verbose && @printf("t = %8.2f   min gap = %.3e\n", t, gap)
        end
    end

    # the final configuration is saved unless it coincides with a save time
    abs(t - (next_save - save_every)) > 1e-12 && save(t, xcm, θcm)
    verbose && contact && @printf("contact at t = %.4f (gap %.2e)\n", t, gap)
    return xcm, θcm, t, contact
end
