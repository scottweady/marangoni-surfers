
"""
    surfer_boundary(a, b, s₀, Ns; quad=:phi)

Dirichlet and Neumann arcs of an elliptical surfer in the body frame.

The ellipse `γ(s) = a cos s + i b sin s` is split into the Neumann arc `s ∈ [-s₀, s₀]` and
the Dirichlet arc `s ∈ [s₀, 2π - s₀]`.

With `quad = :phi` (default) both arcs are [`ellipse_arc`](@ref)s, on which the density
near the Dirichlet-Neumann junctions (`u/√(1-t²) + v` on the Dirichlet arc, half-integer
powers on the Neumann arc) is smooth, and each arc carries its coupling to the other. With
`quad = :mixed` the Dirichlet arc uses `:weighted` and the Neumann arc `:smooth`, which
resolves only the leading junction singularity.

# Arguments
- `a`, `b`: semi-axes along `x` and `y`.
- `s₀`: half-angle of the Neumann arc, in the ellipse parameter `s`.
- `Ns`: number of quadrature nodes on each arc.
- `quad`: arc discretization, `:phi` (default) or `:mixed`.

# Returns
- `(Γ_D, Γ_N)`: the Dirichlet and Neumann arc `Boundary`s.
"""
function surfer_boundary(a, b, s₀, Ns; quad=:phi)

    if quad == :phi
        Γ_D = ellipse_arc(a, b, s₀, 2π - s₀, Ns)
        Γ_N = ellipse_arc(a, b, -s₀, s₀, Ns)
        V_DN = _laplace2d_single_layer_phi_cross(a, b, -s₀, s₀, Γ_N.t, Γ_N.w, Γ_D.s)
        K_ND = _laplace2d_double_layer_adjoint_phi(a, b, s₀, 2π - s₀, Γ_D.t, Γ_D.w, Γ_N.s)
        return with_coupling(Γ_D, V_DN), with_coupling(Γ_N, K_ND)
    elseif quad == :mixed
        γ(s)   =  a * cos(s) + im * b * sin(s)
        γs(s)  = -a * sin(s) + im * b * cos(s)
        γss(s) = -a * cos(s) - im * b * sin(s)
        Γ_D = openarc(γ, γs, γss, Ns; dom=[s₀, 2π - s₀], quad=:weighted)
        Γ_N = openarc(γ, γs, γss, Ns; dom=[-s₀, s₀], quad=:smooth)
        return Γ_D, Γ_N
    else
        error("Unknown surfer quadrature: $quad")
    end

end

"""
    surfer(a, b, s₀, Ns, M; quad=:phi)

Stokes domain and surfactant arcs of an elliptical surfer in the body frame.

# Arguments
- `a`, `b`: semi-axes along `x` and `y`.
- `s₀`: half-angle of the Neumann arc, in the ellipse parameter `s`.
- `Ns`: number of quadrature nodes on each arc.
- `M`: expansion order of the Stokes domain.
- `quad`: arc discretization, passed to [`surfer_boundary`](@ref) (default `:phi`).

# Returns
- `(Ω, Γ_D, Γ_N)`: the Stokes domain and the arcs from [`surfer_boundary`](@ref).
"""
surfer(a, b, s₀, Ns, M; quad=:phi) = (ellipse(a, b, M), surfer_boundary(a, b, s₀, Ns; quad=quad)...)

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
    _boundary_interp(Ω, s)

Matrix taking the boundary trace of a domain field to its values at arc parameters `s`, by
trigonometric interpolation.

# Arguments
- `Ω`: Stokes domain whose boundary the field is traced to.
- `s`: arc parameter values at which to evaluate the traced field.

# Returns
- Matrix `A` with `A * trace(f, Ω)` the traced field values at `s`.
"""
function _boundary_interp(Ω, s)
    N     = size(Ω.z, 2)
    freqs = fftfreq(N) * N
    return cis.(s .* transpose(freqs)) * (fft(Matrix{ComplexF64}(I, N, N), 1) ./ N)
end

"""
    _boundary_eval(f, Ω, s)

Evaluate a domain field `f` at arc parameter values `s` by tracing to the boundary and
interpolating with [`_boundary_interp`](@ref).

# Arguments
- `f`: field sampled at the quadrature nodes of `Ω`.
- `Ω`: Stokes domain whose boundary `f` is traced to.
- `s`: arc parameter values at which to evaluate the traced field.

# Returns
- The traced field values at `s`.
"""
_boundary_eval(f, Ω, s) = _boundary_interp(Ω, s) * vec(trace(f, Ω))

"""
    _reciprocal_kernel(μ, Ω, Γ_D, Γ_N)

Boundary kernel that moves a domain integral of the Marangoni flow onto the surfactant arcs,

    ∫_Ω μ · ∇ℬ[σ] dA = ∫_{Γ_D} h σ_D ds + ∫_{Γ_N} h σ_N ds,

with `h = 𝒮[μ₁x + μ₂y] - (x 𝒮[μ₁] + y 𝒮[μ₂])` traced to the arcs. The arcs lie on ∂Ω, where
`∇ℬ[σ]` is not smooth, while `h` is.

# Arguments
- `μ`: tuple `(μ₁, μ₂)` of fields sampled at the quadrature nodes of `Ω`, with the `1/w`
  edge behaviour of a force density.
- `Ω`: Stokes domain.
- `Γ_D`, `Γ_N`: reference Dirichlet and Neumann arc boundaries, on ∂Ω.

# Returns
- `(h_D, h_N)`: the kernel at the nodes of `Γ_D` and `Γ_N`.
"""
function _reciprocal_kernel(μ, Ω, Γ_D, Γ_N)
    x, y = real.(Ω.z), imag.(Ω.z)
    h = 𝒮(μ[1] .* x + μ[2] .* y, Ω) - (x .* 𝒮(μ[1], Ω) + y .* 𝒮(μ[2], Ω))
    return _boundary_eval(h, Ω, Γ_D.s), _boundary_eval(h, Ω, Γ_N.s)
end

"""
    marangoni_self_operator(Ω, Γ_D, Γ_N)

Operator giving the PSH coefficients of a surfer's own Marangoni flow `∇ℬ[σ]` on its Stokes
domain, in the body frame.

`psh` approximates the inner products `∫ conj(Y_lm) u dA_ζ / w` by quadrature, which
converges slowly for the self flow since it is not smooth at ∂Ω. The inner products are
instead evaluated exactly on the arcs with the reciprocal kernel of
[`_reciprocal_kernel`](@ref), taking `μ = conj(Y_lm) / (w J)`, where `J` is the area
Jacobian of the map from the unit disk.

The kernel is built for all modes at once in coefficient space: `psh` on the even modes is
the matrix `Y' diag(dw)`, `𝒮` is `K̂_S`, and an even expansion on the rim is evaluated
directly at the arc parameters. The rows for the top degree `l = M` would need `x Y_lm` of
degree `M + 1`, which `Ω` cannot represent, so they are left zero.

# Arguments
- `Ω`: Stokes domain.
- `Γ_D`, `Γ_N`: reference Dirichlet and Neumann arc boundaries, on ∂Ω.

# Returns
- `(Hx, Hy)`: matrices with `vec(psh(∇ℬ[σ]ᵢ, Ω)) ≈ Hᵢ * [σ_D; σ_N]` on the even modes
  below the top degree.
"""
function marangoni_self_operator(Ω, Γ_D, Γ_N)

    modes = findall(vec(Ω.even))
    lm    = [(Ω.Lspan[I[1]], Ω.Mspan[I[2]]) for I in CartesianIndices(Ω.shp)[modes]]

    # Even basis at the Stokes nodes, and psh on the even modes as a matrix
    x, y = real.(vec(Ω.z)), imag.(vec(Ω.z))
    J    = vec(Ω.dz ./ Ω.dζ)
    dw   = vec(repeat(vec(Ω.dw), 1, size(Ω.ζ, 2)))
    Y    = reduce(hcat, [vec(ylm(l, m, Ω.ζ)) for (l, m) in lm])
    P    = Y' .* transpose(dw)

    # Coefficients of 𝒮[φ], 𝒮[x φ] and 𝒮[y φ] for φ = conj(Y_lm) / (w J), one column per mode
    φw  = conj.(Y) ./ J
    Sφ  = Ω.K̂_S * (P * φw)
    Sxφ = Ω.K̂_S * (P * (x .* φw))
    Syφ = Ω.K̂_S * (P * (y .* φw))

    # Reciprocal kernel on the arcs, evaluating each even mode at the rim point e^{is}
    s  = [Γ_D.s; Γ_N.s]
    zb = Ω.f.(cis.(s))
    Yb = reduce(hcat, [ylm(l, m, 1.0) .* cis.(m .* s) for (l, m) in lm])
    YS = Yb * Sφ
    Tx = Yb * Sxφ .- real.(zb) .* YS
    Ty = Yb * Syφ .- imag.(zb) .* YS

    keep = [k for (k, (l, _)) in enumerate(lm) if l < Ω.Mℓ]
    dΓ   = [Γ_D.dΓ; Γ_N.dΓ]
    Hx   = zeros(ComplexF64, length(Ω.z), length(s))
    Hy   = zeros(ComplexF64, length(Ω.z), length(s))
    Hx[modes[keep], :] = transpose(Tx[:, keep] .* dΓ)
    Hy[modes[keep], :] = transpose(Ty[:, keep] .* dΓ)
    return Hx, Hy

end

"""
    _marangoni_self_flow(self_op, σ_D, σ_N, Ω)

PSH coefficients of a surfer's own Marangoni flow, in the body frame.

# Arguments
- `self_op`: operator from [`marangoni_self_operator`](@ref).
- `σ_D`, `σ_N`: layer densities on the surfer's Dirichlet and Neumann arcs.
- `Ω`: Stokes domain.

# Returns
- `(b̂x, b̂y)`: coefficient arrays with the shape of `Ω.z`.
"""
function _marangoni_self_flow(self_op, σ_D, σ_N, Ω)
    Hx, Hy = self_op
    σ = [σ_D; σ_N]
    return reshape(Hx * σ, size(Ω.z)), reshape(Hy * σ, size(Ω.z))
end

"""
    surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

Solve the surfactant Laplace problem given reference arc boundaries.

`Γ_D` and `Γ_N` are reference (body-frame) Dirichlet and Neumann arc boundaries. They
are rigidly placed at each `(xcm[k], θcm[k])` — the self-interaction
matrices `V`, `K′` and the couplings between a particle's own arcs are reused since the
Laplace kernel is translation/rotation invariant.

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

    # each particle's Dirichlet and Neumann arcs meet at its junctions
    partner = [isodd(k) ? k + 1 : k - 1 for k in eachindex(Γ)]
    σ, ϕinf = laplace_solve(Γ, bc, rhs; partner=partner)
    return σ, ϕinf, Γ

end

"""
    marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=1e-6, self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

Solve for the translational velocity `U` and angular velocity `ω` of a collection of
Marangoni surfers via the mobility formulation.

`Ω` is the Stokes domain boundary, `Γ_D` and `Γ_N` are the Dirichlet and Neumann arc
boundaries for the surfactant solve, and `g` is the Dirichlet BC as a function of arc
parameter `s`.

The Marangoni flow on each particle is split into the flow from the other particles, which
is smooth there and sampled at the Stokes nodes, and the particle's own flow, which is not
smooth at its edge and is projected with `self_op`.

# Arguments
- `xcm`: vector of complex particle centers.
- `θcm`: vector of particle orientations, in radians.
- `Ω`: Stokes domain boundary (shared by all particles).
- `Γ_D::Boundary`: reference Dirichlet arc boundary.
- `Γ_N::Boundary`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.
- `tol`: FMM/GMRES accuracy tolerance (default `1e-6`).
- `self_op`: self-flow operator from [`marangoni_self_operator`](@ref); pass a prebuilt one
  when solving repeatedly with the same `Ω`, `Γ_D`, `Γ_N`.

# Returns
- `(U, σ, Γ, ϕ∞, f)`: per-particle rigid-body velocities `U[i] = (ux, uy, ω)`, layer
  densities, placed boundaries, the far-field constant, and the Stokes force densities
  `f[i]` on each particle's domain, in its body frame.
"""
function marangoni_solve(xcm, θcm, Ω, Γ_D::Boundary, Γ_N::Boundary, g; tol=1e-6,
                         self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

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
        own = 2i-1:2i

        # Flow from the other particles: remove the particle's own arcs, sampled with the
        # same quadrature so the near-edge error cancels
        sx, sy = ∇ℬ(pts[i], σ[own], Γ[own]; eps=tol)
        cx = reshape(bx[ri], shp) .- sx
        cy = reshape(by[ri], shp) .- sy

        # Own flow from its exact PSH coefficients, rotated to the world frame
        ux, uy = real.(ipsh(_marangoni_self_flow(self_op, σ[own]..., Ω), Ω))
        c, s = cos(θcm[i]), sin(θcm[i])
        uinf[i] = (cx .+ c .* ux .- s .* uy, cy .+ s .* ux .+ c .* uy)
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
        h_D, h_N = _reciprocal_kernel(𝒮_st⁻¹(e, Ω), Ω, Γ_D, Γ_N)
        return integral(h_D .* σ_D, Γ_D) + integral(h_N .* σ_N, Γ_N)
    end

    F = [_marangoni_force((1.0, 0.0)), _marangoni_force((0.0, 1.0)), _marangoni_force((-y, x))]
    M = stokes_mobility_matrix(Ω)
    U = real.(M \ F)

    return U, σ, Γ, ϕinf

end

"""
    marangoni_multipoles(Ω, Γ_D, Γ_N, g; self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

Compute the Marangoni multipole moments for a single surfer.

The force density solves `𝒮_st f = U + ω × r - ∇ℬ[σ]` in PSH coefficient space, with the
coefficients of the Marangoni flow taken from `self_op`.

# Arguments
- `Ω`: Stokes domain boundary.
- `Γ_D`: reference Dirichlet arc boundary.
- `Γ_N`: reference Neumann arc boundary.
- `g`: Dirichlet boundary data, as a function of arc parameter `s`.
- `self_op`: self-flow operator from [`marangoni_self_operator`](@ref).

# Returns
- `NamedTuple` with fields `velocity`, `force_monopole`, `force_dipole`,
  `surf_monopole`, `surf_dipole`.
"""
function marangoni_multipoles(Ω, Γ_D, Γ_N, g; self_op=marangoni_self_operator(Ω, Γ_D, Γ_N))

    U, σ, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)

    ucm, ωcm = U[1:2], U[3]
    α_D, α_N = σ

    z    = Ω.z
    x, y = real.(z), imag.(z)

    b̂x, b̂y = _marangoni_self_flow(self_op, α_D, α_N, Ω)
    r̂x, r̂y = psh((ucm[1] .- y*ωcm, ucm[2] .+ x*ωcm), Ω)
    f̂  = solve(Ω.K̂_G, (r̂x .- b̂x, r̂y .- b̂y), Ω, parity=:even)
    f  = ipsh(f̂, Ω, parity=:even) ./ Ω.w

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
