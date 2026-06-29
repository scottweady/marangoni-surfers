
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


"Laplace double-layer potential: 𝒦(ζ) = 1/(2π) ∫_Γ n̂(y)·(ζ-y)/|ζ-y|² μ(y) ds(y).  Valid for ζ off Γ."
function 𝒦(ζ::Number, μ, Γ::Boundary)
    n̂ = (imag.(Γ.γt) .- im .* real.(Γ.γt)) ./ abs.(Γ.γt)
    r  = ζ .- Γ.γ
    return (1/(2π)) * sum(real.(conj.(n̂) ./ conj.(r)) .* μ .* Γ.dΓ)
end

function 𝒦(ζ, μ, Γ::Boundary)
    n̂  = (imag.(Γ.γt) .- im .* real.(Γ.γt)) ./ abs.(Γ.γt)
    shp = size(ζ)
    r   = vec(ζ) .- transpose(Γ.γ)
    A   = real.(conj.(transpose(n̂)) ./ conj.(r)) .* transpose(μ .* Γ.dΓ)
    return reshape((1/(2π)) * sum(A, dims=2), shp)
end

"Helmholtz single-layer potential: ℋ(ζ) = 1/(2π) ∫_Γ K₀(k|ζ-γ′|) τ ds′.  Valid for ζ off Γ."
function ℋ(ζ, τ, Γ::Boundary; k=1.0)
    return _smooth_evaluation(ζ, τ, r -> besselk(0, k*abs(r)) / (2π), Γ)
end

"Biharmonic single-layer potential: ℬ(ζ) = 1/(4π) ∫_Γ |ζ-γ'| σ ds'."
function ℬ(ζ, σ, Γ::Boundary)
    return _smooth_evaluation(ζ, σ, r -> abs(r)/(4π), Γ)
end


"Gradient of Laplace single-layer: ∇𝒱(ζ) = -1/(2π) ∫_Γ (ζ-γ′)/|ζ-γ′|² σ ds′.  Valid for ζ off Γ."
function ∇𝒱(ζ, σ, Γ::Boundary)
    Vx = _smooth_evaluation(ζ, σ, r -> -real(r) / (2π * abs2(r)), Γ)
    Vy = _smooth_evaluation(ζ, σ, r -> -imag(r) / (2π * abs2(r)), Γ)
    return Vx, Vy
end

"""
Gradient of Laplace double-layer: ∇𝒦(ζ) = 1/(2π) ∫_Γ (I/|r|² - 2rr̂/|r|⁴)·n̂(y) μ(y) ds(y),  r = ζ-y.

Returns (Kx, Ky) arrays matching the shape of ζ.  Valid for ζ off Γ.
"""
function ∇𝒦(ζ::Number, μ, Γ::Boundary)
    n̂     = (imag.(Γ.γt) .- im .* real.(Γ.γt)) ./ abs.(Γ.γt)
    r     = ζ .- Γ.γ
    r2    = abs2.(r)
    ndotr = real.(conj.(n̂) .* r)
    w     = μ .* Γ.dΓ
    Kx = (1/(2π)) * sum((real.(n̂) ./ r2 .- 2 .* real.(r) .* ndotr ./ r2.^2) .* w)
    Ky = (1/(2π)) * sum((imag.(n̂) ./ r2 .- 2 .* imag.(r) .* ndotr ./ r2.^2) .* w)
    return Kx, Ky
end

function ∇𝒦(ζ, μ, Γ::Boundary)
    n̂     = (imag.(Γ.γt) .- im .* real.(Γ.γt)) ./ abs.(Γ.γt)
    shp   = size(ζ)
    r     = vec(ζ) .- transpose(Γ.γ)
    r2    = abs2.(r)
    ndotr = real.(conj.(transpose(n̂)) .* r)
    w     = transpose(μ .* Γ.dΓ)
    Kx = reshape((1/(2π)) * sum((real.(transpose(n̂)) ./ r2 .- 2 .* real.(r) .* ndotr ./ r2.^2) .* w, dims=2), shp)
    Ky = reshape((1/(2π)) * sum((imag.(transpose(n̂)) ./ r2 .- 2 .* imag.(r) .* ndotr ./ r2.^2) .* w, dims=2), shp)
    return Kx, Ky
end

"Gradient of Helmholtz single-layer: ∇ℋ(ζ) = -k/(2π) ∫_Γ K₁(k|ζ-γ′|)/|ζ-γ′| (ζ-γ′) τ ds′.  Valid for ζ off Γ."
function ∇ℋ(ζ, τ, Γ::Boundary; k=1.0)
    function kernel_x(r)
        ρ = abs(r)
        ρ < eps() && return 0.0
        return -k * real(r) * besselk(1, k*ρ) / (2π * ρ)
    end
    function kernel_y(r)
        ρ = abs(r)
        ρ < eps() && return 0.0
        return -k * imag(r) * besselk(1, k*ρ) / (2π * ρ)
    end
    Hx = _smooth_evaluation(ζ, τ, kernel_x, Γ)
    Hy = _smooth_evaluation(ζ, τ, kernel_y, Γ)
    return Hx, Hy
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
