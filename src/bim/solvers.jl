
# +1/2 jump term for the Neumann diagonal block.
# Open arcs (Chebyshev nodes, |s| ≤ 1) use the 1/√(1-t²) weighted jump; closed curves use 1/2 I.
function _jump_term(Γ::Boundary)
    if all(abs.(Γ.s) .<= 1.0 + 1e-10)
        return Diagonal(1 ./ (2 .* sqrt.(1 .- Γ.t.^2)))
    else
        return 0.5 * I(length(Γ))
    end
end

"""
Assemble and solve the mixed Dirichlet/Neumann Laplace BIM system for a collection of boundaries.

`types[i]` is `:D` (Dirichlet, single-layer ansatz) or `:N` (Neumann, double-layer ansatz).
An extra unknown ϕ∞ (the far-field constant) is appended and fixed by enforcing ∫_all σ ds = 0.
Returns `(densities, ϕ∞)` where `densities[i]` is the layer density on boundary i.
"""
function laplace2d(boundaries::Vector{Boundary}, types::Vector{Symbol}, g::Vector)
    n     = length(boundaries)
    Ns    = length.(boundaries)
    offs  = [0; cumsum(Ns)]
    N_tot = offs[end]

    A   = zeros(N_tot + 1, N_tot + 1)
    rhs = zeros(N_tot + 1)

    for i in 1:n
        ri = offs[i]+1 : offs[i+1]
        rhs[ri] = g[i]
        for j in 1:n
            cj = offs[j]+1 : offs[j+1]
            if i == j
                A[ri, cj] = types[i] == :D ? boundaries[i].V :
                                              Matrix(_jump_term(boundaries[i])) + boundaries[i].V′
            else
                A[ri, cj] = types[i] == :D ? _laplace2d_single_layer(boundaries[j], boundaries[i]) :
                                              _laplace2d_double_layer_adjoint(boundaries[j], boundaries[i])
            end
        end
        if types[i] == :D
            A[ri, end] .= 1.0
        end
        A[end, ri] = abs.(boundaries[i].γt) .* boundaries[i].w
    end

    σ_aug     = A \ rhs
    densities = [σ_aug[offs[i]+1:offs[i+1]] for i in 1:n]
    return densities, σ_aug[end]
end

"Two-boundary convenience wrapper: Dirichlet arc Γ_D and Neumann arc Γ_N."
function laplace2d(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)
    densities, ϕinf = laplace2d([Γ_D, Γ_N], [:D, :N], [g_D, g_N])
    return densities[1], densities[2], ϕinf
end

_place(Γ::Boundary, xcm, θcm) = (eθ = exp(1im*θcm);
    Boundary(Γ.t, Γ.s, xcm .+ eθ.*Γ.γ, eθ.*Γ.γt, eθ.*Γ.γtt, Γ.w, Γ.V, Γ.V′))

"""
Solve the surfactant Laplace problem given reference arc boundaries.

`Γ_D` and `Γ_N` are reference (body-frame) Dirichlet and Neumann arc boundaries.
They are rigidly placed at each `(centers[k], orientations[k])` — the self-interaction
matrices `V`, `V′` are reused since the Laplace kernel is translation/rotation invariant.
`g` is the Dirichlet BC as a function of arc parameter s.
Returns `(densities, boundaries, ϕ∞)`.
"""
function surfactant_solve(centers, orientations, Γ_D::Boundary, Γ_N::Boundary, g)
    boundaries = Boundary[]
    types      = Symbol[]
    rhs        = Vector{Float64}[]

    for (xcm, θcm) in zip(centers, orientations)
        γd = _place(Γ_D, xcm, θcm)
        γn = _place(Γ_N, xcm, θcm)
        push!(boundaries, γd, γn)
        push!(types, :D, :N)
        push!(rhs, g.(γd.s), zeros(length(γn)))
    end

    densities, ϕinf = laplace2d(boundaries, types, rhs)
    return densities, boundaries, ϕinf
end

"""
Solve for the translational velocity U and angular velocity ω of a single ellipse surfer.

`Ns` is the number of arc quadrature nodes for the surfactant solve.
`Nq` is the number of quadrature nodes for the Stokes domain discretization.
Returns `(U, ω, α_D, α_N, Γ_D, Γ_N, Ω)`. Pass the last five to `marangoni_multipoles` if needed.
"""
function marangoni_mobility(Ns, Nq, a, b, α, g)
    Ω      = ellipse(a, b, Nq)
    z      = Ω.z
    x, y   = real.(z), imag.(z)

    γ(s)   =  a * cos.(s) + im * b * sin.(s)
    γs(s)  = -a * sin.(s) + im * b * cos.(s)
    γss(s) = -a * cos.(s) - im * b * sin.(s)

    # ── Surfactant solve ──────────────────────────────────────────────────────

    Γ_D  = arc(γ, γs, γss, Ns; dom=[α, 2π-α])
    Γ_N  = arc(γ, γs, γss, Ns; dom=[-α, α])
    α_D, α_N, _ = laplace2d(Γ_D, Γ_N, g.(Γ_D.s), zeros(Ns))

    # ── Stokes mobility solve ─────────────────────────────────────────────────

    # reciprocal_vel(μ): boundary trace of S[μ₁x + μ₂y] - xS[μ₁] - yS[μ₂]
    function reciprocal_vel(μ)
        tmp = 𝒮(μ[1] .* x + μ[2] .* y, Ω) - (x .* 𝒮(μ[1], Ω) + y .* 𝒮(μ[2], Ω))
        return boundary_eval(tmp, Ω, Γ_D.s), boundary_eval(tmp, Ω, Γ_N.s)
    end

    M_mob = stokes_mobility_matrix(Ω)
    rv1   = reciprocal_vel(𝒮_st⁻¹((1, 0), Ω))
    rv2   = reciprocal_vel(𝒮_st⁻¹((0, 1), Ω))
    rv3   = reciprocal_vel(𝒮_st⁻¹((-y, x), Ω))

    F = [integral(rv1[1] .* α_D, Γ_D) + integral(rv1[2] .* α_N, Γ_N),
         integral(rv2[1] .* α_D, Γ_D) + integral(rv2[2] .* α_N, Γ_N),
         integral(rv3[1] .* α_D, Γ_D) + integral(rv3[2] .* α_N, Γ_N)]

    result = M_mob \ F
    U  = real(result[1:2])
    ω  = real(result[3])

    return U, ω, α_D, α_N, Γ_D, Γ_N, Ω
end

"""
Compute force and surfactant multipoles from the output of `marangoni_mobility`.

Returns a NamedTuple with fields `force_monopole`, `force_dipole`, `surf_monopole`, `surf_dipole`.
Each is a vector or matrix of the corresponding moments.
"""
function marangoni_multipoles(U, ω, α_D, α_N, Γ_D, Γ_N, Ω)
    z    = Ω.z
    x, y = real.(z), imag.(z)

    bx = ∇ℬ(z, α_D, Γ_D)[1] .+ ∇ℬ(z, α_N, Γ_N)[1]
    by = ∇ℬ(z, α_D, Γ_D)[2] .+ ∇ℬ(z, α_N, Γ_N)[2]
    f  = 𝒮_st⁻¹((U[1] .- y*ω .- bx, U[2] .+ x*ω .- by), Ω)

    fm = [real.(sum(f[1] .* Ω.dz)), real.(sum(f[2] .* Ω.dz))]
    fd = [real.(sum(x .* f[1] .* Ω.dz))  real.(sum(x .* f[2] .* Ω.dz));
          real.(sum(y .* f[1] .* Ω.dz))  real.(sum(y .* f[2] .* Ω.dz))]

    xD, xN = cos.(Γ_D.s), cos.(Γ_N.s)
    yD, yN = sin.(Γ_D.s), sin.(Γ_N.s)

    sm = [integral(xD .* α_D, Γ_D) + integral(xN .* α_N, Γ_N),
          integral(yD .* α_D, Γ_D) + integral(yN .* α_N, Γ_N)]
    sd = [integral(xD.*xD .* α_D, Γ_D) + integral(xN.*xN .* α_N, Γ_N)  integral(xD.*yD .* α_D, Γ_D) + integral(xN.*yN .* α_N, Γ_N);
          integral(yD.*xD .* α_D, Γ_D) + integral(yN.*xN .* α_N, Γ_N)  integral(yD.*yD .* α_D, Γ_D) + integral(yN.*yN .* α_N, Γ_N)]

    return (force_monopole=fm, force_dipole=fd, surf_monopole=sm, surf_dipole=sd)
end

function marangoni_solve(centers, orientations, Ω, Γ_D::Boundary, Γ_N::Boundary, g; eps=1e-6)
    Np = length(centers)

    σs, boundaries, ϕinf = surfactant_solve(centers, orientations, Γ_D, Γ_N, g)

    uslip = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Np)
    for i in 1:Np
        pts = exp(im * orientations[i]) .* Ω.z .+ centers[i]
        bx  = zeros(size(pts))
        by  = zeros(size(pts))
        for j in eachindex(σs)
            bxj, byj = ∇ℬ(pts, σs[j], boundaries[j])
            bx .+= bxj
            by .+= byj
        end
        uslip[i] = (-bx, -by)
    end

    F = [zeros(3) for _ in 1:Np]
    _, U = suspension_mobility_solve(centers, orientations, nothing, uslip, F, Ω; eps=eps)
    return U, σs, boundaries, ϕinf
end
