
"""
    surfer_boundary(a, b, s₀, Ns)

Dirichlet and Neumann arcs of an elliptical surfer in the body frame.

The ellipse `γ(s) = a cos s + i b sin s` is split into the Neumann arc `s ∈ [-s₀, s₀]` and
the Dirichlet arc `s ∈ [s₀, 2π - s₀]`. The Dirichlet arc uses `:open_weight` quadrature
and the Neumann arc uses `:open`.

# Arguments
- `a`, `b`: semi-axes along `x` and `y`.
- `s₀`: half-angle of the Neumann arc, in the ellipse parameter `s`.
- `Ns`: number of quadrature nodes on each arc.

# Returns
- `(Γ_D, Γ_N)`: the Dirichlet and Neumann arc `Boundary`s.
"""
function surfer_boundary(a, b, s₀, Ns)

    γ(s)   =  a * cos(s) + im * b * sin(s)
    γs(s)  = -a * sin(s) + im * b * cos(s)
    γss(s) = -a * cos(s) - im * b * sin(s)

    Γ_D = openarc(γ, γs, γss, Ns; dom=[s₀, 2π - s₀], quad=:open_weight)
    Γ_N = openarc(γ, γs, γss, Ns; dom=[-s₀, s₀], quad=:open)
    return Γ_D, Γ_N

end

"""
    surfer(a, b, s₀, Ns, M)

Stokes domain and surfactant arcs of an elliptical surfer in the body frame.

# Arguments
- `a`, `b`: semi-axes along `x` and `y`.
- `s₀`: half-angle of the Neumann arc, in the ellipse parameter `s`.
- `Ns`: number of quadrature nodes on each arc.
- `M`: expansion order of the Stokes domain.

# Returns
- `(Ω, Γ_D, Γ_N)`: the Stokes domain and the arcs from [`surfer_boundary`](@ref).
"""
surfer(a, b, s₀, Ns, M) = (ellipse(a, b, M), surfer_boundary(a, b, s₀, Ns)...)

"""
    integral(f, Ω)

Domain integral `∫_Ω f dA = ∑ fᵢ dzᵢ`.

# Arguments
- `f`: field sampled at the quadrature nodes of `Ω`.
- `Ω::Domain`: Stokes domain to integrate over.

# Returns
- The scalar value of the integral.
"""
integral(f, Ω::Domain) = sum(f .* Ω.dz)

"""
    _boundary_eval(f, Ω, s)

Evaluate a domain field `f` at arc parameter values `s` by tracing to the boundary and
FFT-interpolating.

# Arguments
- `f`: field sampled at the quadrature nodes of `Ω`.
- `Ω`: Stokes domain whose boundary `f` is traced to.
- `s`: arc parameter values at which to evaluate the traced field.

# Returns
- The traced field values at `s`.
"""
function _boundary_eval(f, Ω, s)
    vals  = trace(f, Ω)
    N     = length(vals)
    f̂     = fft(vals) ./ N
    freqs = fftfreq(N) * N
    return real.(sum(f̂[k] .* cis.(freqs[k] .* s) for k in 1:N))
end

"""
    surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

Solve the surfactant Laplace problem given reference arc boundaries.

`Γ_D` and `Γ_N` are reference (body-frame) Dirichlet and Neumann arc boundaries. They
are rigidly placed at each `(xcm[k], θcm[k])` — the self-interaction
matrices `V`, `K′` are reused since the Laplace kernel is translation/rotation invariant.

# Arguments
- `xcm`: vector of complex particle centers.
- `θcm`: vector of particle orientations, in radians.
- `Γ_D::Boundary`: reference Dirichlet arc boundary.
- `Γ_N::Boundary`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.

# Returns
- `(σ, ϕ∞, Γ)`: layer densities per placed boundary, the far-field constant, and the
  vector of placed `Boundary`s (alternating `Γ_D`, `Γ_N` per particle).
"""
function surfactant_solve(xcm, θcm, Γ_D::Boundary, Γ_N::Boundary, g)

    Γ = Boundary[]
    bc = Symbol[]
    rhs = Vector{Float64}[]

    for (xcm_i, θcm_i) in zip(xcm, θcm)
        Γ_D_i = _place(Γ_D, xcm_i, θcm_i)
        Γ_N_i = _place(Γ_N, xcm_i, θcm_i)
        push!(Γ, Γ_D_i, Γ_N_i)
        push!(bc, :D, :N)
        push!(rhs, g.(Γ_D_i.s), zeros(length(Γ_N_i.s)))
    end

    σ, ϕinf = laplace_solve(Γ, bc, rhs)
    return σ, ϕinf, Γ

end

"""
    marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=1e-6)

Solve for the translational velocity `U` and angular velocity `ω` of a collection of
Marangoni surfers via the mobility formulation.

`Ω` is the Stokes domain boundary, `Γ_D` and `Γ_N` are the Dirichlet and Neumann arc
boundaries for the surfactant solve, and `g` is the Dirichlet BC as a function of arc
parameter `s`.

# Arguments
- `xcm`: vector of complex particle centers.
- `θcm`: vector of particle orientations, in radians.
- `Ω`: Stokes domain boundary (shared by all particles).
- `Γ_D::Boundary`: reference Dirichlet arc boundary.
- `Γ_N::Boundary`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.
- `tol`: FMM/GMRES accuracy tolerance (default `1e-6`).

# Returns
- `(U, σ, Γ, ϕ∞, f)`: per-particle rigid-body velocities `U[i] = (ux, uy, ω)`, layer
  densities, placed boundaries, the far-field constant, and the Stokes force densities
  `f[i]` on each particle's domain, in its body frame.
"""
function marangoni_solve(xcm, θcm, Ω, Γ_D::Boundary, Γ_N::Boundary, g; tol=1e-6)

    Nparticles = length(xcm)

    σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

    pts  = [exp(im * θcm[i]) .* Ω.z .+ xcm[i] for i in 1:Nparticles]
    Npts     = length(Ω.z)
    offs     = [0; cumsum(fill(Npts, Nparticles))]
    bx, by = ∇ℬ(vcat(vec.(pts)...), σ, Γ; eps=tol)

    uinf = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Nparticles)
    for i in 1:Nparticles
        ri = offs[i]+1:offs[i+1]
        shp = size(pts[i])
        uinf[i] = (reshape(copy(bx[ri]), shp), reshape(copy(by[ri]), shp))
    end

    F = [zeros(3) for _ in 1:Nparticles]
    f, U = suspension_mobility_solve(xcm, θcm, uinf, nothing, F, Ω; eps=tol, gmrestol=tol)

    return U, σ, Γ, ϕinf, f

end

"""
    surface_fields(z, xcm, θcm, Ω, σ, Γ, ϕ∞, f; inside=nothing, tol=1e-6)

Surfactant concentration and surface velocity at points `z`, from the solution of the
multi-particle [`marangoni_solve`](@ref).

The concentration is the single-layer potential of `σ` plus `ϕ∞`. The velocity is the
Marangoni flow `∇ℬ[σ]` plus, for each particle, the Stokes single layer of its force
density, evaluated in the body frame and rotated back.

# Arguments
- `z`: array of complex evaluation points.
- `xcm`: vector of complex particle centers.
- `θcm`: vector of particle orientations, in radians.
- `Ω`: Stokes domain boundary (shared by all particles).
- `σ`, `Γ`, `ϕ∞`, `f`: layer densities, placed boundaries, far-field constant, and force
  densities returned by [`marangoni_solve`](@ref).
- `inside`: predicate on body-frame points marking the particle interior, where the
  fields are set to `NaN` (default `nothing`, no masking).
- `tol`: FMM accuracy tolerance for `∇ℬ` (default `1e-6`).

# Returns
- `(ϕ, u, v)`: concentration and velocity components, with the same shape as `z`.
"""
function surface_fields(z, xcm, θcm, Ω, σ, Γ, ϕinf, f; inside=nothing, tol=1e-6)

    z = ComplexF64.(z)
    ϕ = sum(𝒱(z, σ[2k-1], Γ[2k-1]) .+ 𝒱(z, σ[2k], Γ[2k]) for k in eachindex(xcm)) .+ ϕinf
    u, v = ∇ℬ(z, σ, Γ; eps=tol)

    for k in eachindex(xcm)
        z_body = exp(-im * θcm[k]) .* (z .- xcm[k])
        vel = 𝒮_st(z_body, f[k], Ω)
        u_b, v_b = real.(first.(vel)), real.(last.(vel))
        c, s = cos(θcm[k]), sin(θcm[k])
        u .+= c .* u_b .- s .* v_b
        v .+= s .* u_b .+ c .* v_b
        if inside !== nothing
            mask = inside.(z_body)
            ϕ[mask] .= NaN
            u[mask] .= NaN
            v[mask] .= NaN
        end
    end

    return ϕ, u, v

end

"""
    marangoni_solve(Ω, Γ_D, Γ_N, g)

Single-particle Marangoni solve via the mobility formulation.

# Arguments
- `Ω`: Stokes domain boundary.
- `Γ_D`: reference Dirichlet arc boundary.
- `Γ_N`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.

# Returns
- `(U, σ, Γ, ϕ∞)`: rigid-body velocity `U = (ux, uy, ω)`, layer densities `(σ_D, σ_N)`,
  placed boundaries, and the far-field constant.
"""
function marangoni_solve(Ω, Γ_D, Γ_N, g)

    σ, ϕinf, Γ = surfactant_solve([0.0im], [0.0], Γ_D, Γ_N, g)
    σ_D, σ_N = σ[1], σ[2]
    x, y = real.(Ω.z), imag.(Ω.z)

    function _marangoni_force(e)
        μ = 𝒮_st⁻¹(e, Ω)
        h = 𝒮(μ[1] .* x + μ[2] .* y, Ω) - (x .* 𝒮(μ[1], Ω) + y .* 𝒮(μ[2], Ω))
        h_D = _boundary_eval(h, Ω, Γ_D.s)
        h_N = _boundary_eval(h, Ω, Γ_N.s)
        return integral(h_D .* σ_D, Γ_D) + integral(h_N .* σ_N, Γ_N)
    end

    F = [_marangoni_force((1.0, 0.0)), _marangoni_force((0.0, 1.0)), _marangoni_force((-y, x))]
    M = stokes_mobility_matrix(Ω)
    U = real.(M \ F)

    return U, σ, Γ, ϕinf

end

"""
    marangoni_multipoles(Ω, Γ_D, Γ_N, g)

Compute the Marangoni multipole moments for a single surfer.

# Arguments
- `Ω`: Stokes domain boundary.
- `Γ_D`: reference Dirichlet arc boundary.
- `Γ_N`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.

# Returns
- `NamedTuple` with fields `velocity`, `force_monopole`, `force_dipole`,
  `surf_monopole`, `surf_dipole`.
"""
function marangoni_multipoles(Ω, Γ_D, Γ_N, g)

    U, σ, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)

    ucm, ωcm = U[1:2], U[3]
    α_D, α_N = σ

    z    = Ω.z
    x, y = real.(z), imag.(z)

    bx = ∇ℬ(z, α_D, Γ_D)[1] .+ ∇ℬ(z, α_N, Γ_N)[1]
    by = ∇ℬ(z, α_D, Γ_D)[2] .+ ∇ℬ(z, α_N, Γ_N)[2]
    f  = 𝒮_st⁻¹((ucm[1] .- y*ωcm .- bx, ucm[2] .+ x*ωcm .- by), Ω)

    fm = [real.(sum(f[1] .* Ω.dz)), real.(sum(f[2] .* Ω.dz))]
    fd = [real.(sum(x .* f[1] .* Ω.dz))  real.(sum(x .* f[2] .* Ω.dz));
          real.(sum(y .* f[1] .* Ω.dz))  real.(sum(y .* f[2] .* Ω.dz))]

    xD, xN = real.(Γ_D.γ), real.(Γ_N.γ)
    yD, yN = imag.(Γ_D.γ), imag.(Γ_N.γ)

    sm = [integral(xD .* α_D, Γ_D) + integral(xN .* α_N, Γ_N),
          integral(yD .* α_D, Γ_D) + integral(yN .* α_N, Γ_N)]
    sd = [integral(xD.*xD .* α_D, Γ_D) + integral(xN.*xN .* α_N, Γ_N)  integral(xD.*yD .* α_D, Γ_D) + integral(xN.*yN .* α_N, Γ_N);
          integral(yD.*xD .* α_D, Γ_D) + integral(yN.*xN .* α_N, Γ_N)  integral(yD.*yD .* α_D, Γ_D) + integral(yN.*yN .* α_N, Γ_N)]

    return (velocity=U, force_monopole=fm, force_dipole=fd, surf_monopole=sm, surf_dipole=sd)

end
