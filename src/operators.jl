
# Workhorse: ∑ⱼ kernel(ζ - γⱼ) σⱼ |γ′ⱼ| wⱼ for scalar or array targets.
function _smooth_evaluation(ζ::Number, σ, kernel, Γ::Boundary)
    δz = ζ .- Γ.γ
    return sum(kernel.(δz) .* σ .* Γ.dΓ)
end

function _smooth_evaluation(ζ, σ, kernel, Γ::Boundary)
    shp = size(ζ)
    δz  = vec(ζ) .- transpose(Γ.γ)
    A   = kernel.(δz) .* transpose(σ .* Γ.dΓ)
    return reshape(sum(A, dims=2), shp)
end

"Laplace single-layer potential: 𝒱(ζ) = -1/(2π) ∫_Γ log|ζ-γ'| σ ds'."
function 𝒱(ζ, σ, Γ::Boundary)
    return _smooth_evaluation(ζ, σ, r -> -(1/(2π)) * log(abs(r)), Γ)
end

"Boundary-to-boundary single-layer: evaluate 𝒱 at the nodes of Γ_tgt."
function 𝒱(Γ_tgt::Boundary, σ, Γ_src::Boundary)
    return 𝒱(Γ_tgt.γ, σ, Γ_src)
end

"Biharmonic single-layer potential: ℬ(ζ) = 1/(4π) ∫_Γ |ζ-γ'| σ ds'."
function ℬ(ζ, σ, Γ::Boundary)
    return _smooth_evaluation(ζ, σ, r -> abs(r)/(4π), Γ)
end

"""
Gradient of the biharmonic potential: ∇ℬ(ζ) = 1/(4π) ∫_Γ (ζ-γ')/|ζ-γ'| σ ds'.

Accepts a single Boundary or a Vector{Boundary}.  Dispatches to direct
`_smooth_evaluation` when `length(ζ) * total_sources ≤ 4000`, otherwise
uses three FMM3D calls.  Returns (bx, by) arrays matching the shape of ζ.
Valid when ζ is well-separated from Γ; use `∇ℬ_arc` near open arcs.
"""
function ∇ℬ(ζ, σs, Γs::Vector{Boundary}; eps=1e-6)

    Nsrc = sum(length(Γ) for Γ in Γs)

    if length(ζ) * Nsrc <= 10000
        bx = _smooth_evaluation(ζ, σs[1], r -> real(r)/(4π*abs(r)), Γs[1])
        by = _smooth_evaluation(ζ, σs[1], r -> imag(r)/(4π*abs(r)), Γs[1])
        for k in 2:length(Γs)
            bx = bx .+ _smooth_evaluation(ζ, σs[k], r -> real(r)/(4π*abs(r)), Γs[k])
            by = by .+ _smooth_evaluation(ζ, σs[k], r -> imag(r)/(4π*abs(r)), Γs[k])
        end
        return bx, by
    end

    q  = vcat([σs[j] .* Γs[j].dΓ for j in eachindex(Γs)]...)
    γs = vcat([Γs[j].γ for j in eachindex(Γs)]...)
    Ns   = length(γs)
    srcs = vcat(real.(γs)', imag.(γs)', zeros(1, Ns))
    ζv   = vec(ζ)
    Nt   = length(ζv)
    tgts = vcat(real.(ζv)', imag.(ζv)', zeros(1, Nt))

    P1  = real.(lfmm3d(eps, srcs; charges=q,              targets=tgts, pgt=1).pottarg)
    P2x = real.(lfmm3d(eps, srcs; charges=real.(γs) .* q, targets=tgts, pgt=1).pottarg)
    P2y = real.(lfmm3d(eps, srcs; charges=imag.(γs) .* q, targets=tgts, pgt=1).pottarg)

    bx = reshape(real.(ζv) .* P1 .- P2x, size(ζ))
    by = reshape(imag.(ζv) .* P1 .- P2y, size(ζ))
    return bx, by
end

function ∇ℬ(ζ, σ, Γ::Boundary; eps=1e-6)
    return ∇ℬ(ζ, [σ], [Γ]; eps=eps)
end
