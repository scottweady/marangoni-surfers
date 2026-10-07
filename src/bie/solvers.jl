
"""
    _place(Γ, xcm, θcm)

Rigidly translate and rotate a reference `Boundary`, reusing its precomputed quadrature,
self matrices and coupling (`V`, `K′` and `coupling` are frame-invariant under rigid motion).

# Arguments
- `Γ::Boundary`: reference boundary, centered at the origin.
- `xcm`: complex translation (new center of mass).
- `θcm`: rotation angle, in radians.

# Returns
- `Boundary` placed at `xcm` with orientation `θcm`.
"""
function _place(Γ::Boundary, xcm, θcm)
    eθ = exp(1im*θcm)
    return Boundary(Γ.quad, Γ.t, Γ.s, xcm .+ eθ.*Γ.γ, eθ.*Γ.γt, eθ.*Γ.γtt, Γ.w, Γ.dΓ, Γ.V, Γ.K′, Γ.coupling)
end

"""
    _jump_term(Γ)

Diagonal jump term for the double-layer boundary integral equation, accounting for the
endpoint density weighting of the open-arc quadratures.

# Arguments
- `Γ::Boundary`: boundary whose `quad` field selects the jump-term form.

# Returns
- Diagonal operator (`Diagonal` or `UniformScaling`-backed matrix) of size `length(Γ)`:
  `1/(2√(1-t²))` for `:weighted`, `1/2` otherwise.
"""
function _jump_term(Γ::Boundary)
    if Γ.quad == :weighted
        return Diagonal(-1 ./ (2 .* sqrt.(1 .- Γ.t.^2)))
    else
        return -0.5 * I(length(Γ))
    end
end

"""
    laplace_solve(Γ, bc, g; decay=true, partner=nothing)

Assemble and solve the mixed Dirichlet/Neumann Laplace BIM system for a collection of
boundaries.

An extra unknown `ϕ∞` (the far-field constant) is appended and fixed by enforcing
`∫_all σ ds = 0`.

# Arguments
- `Γ::Vector{Boundary}`: boundaries in the system.
- `bc::Vector{Symbol}`: boundary condition per boundary — `:D` (Dirichlet, single-layer
  ansatz) or `:N` (Neumann, double-layer ansatz).
- `g::Vector`: boundary data per boundary, matching the length of each `Γ[i]`.
- `decay`: if `true` (default), enforce the far-field decay constraint and solve for
  `ϕ∞`; if `false`, drop the constraint and solve the unaugmented system with `ϕ∞ = 0`.
- `partner`: optional vector giving, for each boundary, the index of the arc it meets at a
  junction (default `nothing`). For such a pair the block uses `Γ[i].coupling`, when set,
  in place of plain quadrature.

# Returns
- `(σ, ϕ∞)`: `σ[i]` is the layer density on boundary `i`; `ϕ∞` is the far-field constant.
"""
function laplace_solve(Γ::Vector{Boundary}, bc::Vector{Symbol}, g::Vector; decay=true, partner=nothing)

    Nparticles = length(Γ)
    offs  = [0; cumsum(length.(Γ))]
    Ntotal = offs[end]

    A   = zeros(Ntotal + 1, Ntotal + 1)
    rhs = zeros(Ntotal + 1)

    for i in 1:Nparticles
        ri = offs[i]+1 : offs[i+1]
        rhs[ri] = g[i]
        for j in 1:Nparticles
            cj = offs[j]+1 : offs[j+1]
            if i == j
                A[ri, cj] = bc[i] == :D ? Γ[i].V : Matrix(_jump_term(Γ[i])) + Γ[i].K′
            elseif partner !== nothing && partner[i] == j && Γ[i].coupling !== nothing
                A[ri, cj] = Γ[i].coupling
            else
                A[ri, cj] = bc[i] == :D ? _laplace2d_single_layer_cross(Γ[j], Γ[i]) : _laplace2d_double_layer_adjoint(Γ[j], Γ[i])
            end
        end

        if bc[i] == :D
            A[ri, end] .= 1.0
        end

        A[end, ri] = Γ[i].dΓ

    end

    if decay
      σ_aug     = A \ rhs
      σ = [σ_aug[offs[i]+1:offs[i+1]] for i in 1:Nparticles]
      λ = σ_aug[end]
    else
      σ = A[1:Ntotal, 1:Ntotal] \ rhs[1:Ntotal]
      σ = [σ[offs[i]+1:offs[i+1]] for i in 1:Nparticles]
      λ = 0.0
    end

    return σ, λ

end

"""
    laplace_solve(Γ_D, Γ_N, g_D, g_N)

Two-boundary convenience wrapper for [`laplace_solve`](@ref): Dirichlet arc `Γ_D` and
Neumann arc `Γ_N`.

# Arguments
- `Γ_D::Boundary`: Dirichlet boundary.
- `Γ_N::Boundary`: Neumann boundary.
- `g_D`: Dirichlet boundary data on `Γ_D`.
- `g_N`: Neumann boundary data on `Γ_N`.

# Returns
- `(σ_D, σ_N, ϕ∞)`: layer densities on `Γ_D` and `Γ_N`, and the far-field constant.
"""
function laplace_solve(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)
    σ, ϕinf = laplace_solve([Γ_D, Γ_N], [:D, :N], [g_D, g_N])
    return σ[1], σ[2], ϕinf
end
