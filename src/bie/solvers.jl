
function _place(Γ::Boundary, xcm, θcm)
    eθ = exp(1im*θcm)
    return Boundary(Γ.quad, Γ.θ, Γ.t, Γ.s, xcm .+ eθ.*Γ.γ, eθ.*Γ.γt, eθ.*Γ.γtt, Γ.w, Γ.dΓ, Γ.V, Γ.K′)
end

function _jump_term(Γ::Boundary)
    if Γ.quad == :bruno
        return Diagonal(1 ./ (2 .* sqrt.(1 .- Γ.t.^2)))
    elseif Γ.quad == :theta
        return Diagonal(1 ./ (2 .* sin.(Γ.θ)))
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
function laplace_solve(Γ::Vector{Boundary}, bc::Vector{Symbol}, g::Vector; decay=true)

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
            else
                A[ri, cj] = bc[i] == :D ? _laplace2d_single_layer(Γ[j], Γ[i]) : _laplace2d_double_layer_adjoint(Γ[j], Γ[i])
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

"Two-boundary convenience wrapper: Dirichlet arc Γ_D and Neumann arc Γ_N."
function laplace_solve(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)
    σ, ϕinf = laplace_solve([Γ_D, Γ_N], [:D, :N], [g_D, g_N])
    return σ[1], σ[2], ϕinf
end