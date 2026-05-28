
using FFTW
using FMM3D
using LinearAlgebra
using QuadGK

# Real inner product and cross product treating complex numbers as 2-vectors.
dot(a::Complex, b::Complex)   = real(a)*real(b) + imag(a)*imag(b)
cross(a::Complex, b::Complex) = real(a)*imag(b) - imag(a)*real(b)

"Apply the pre-built single-layer matrix: (Vσ)[i] = ∫_Γ G(γᵢ, γ') σ(γ') ds'."
function V(σ, Γ::Boundary)
    return Γ.V * σ
end

"""
Solve Vσ = g on a closed curve, augmented to remove the log-capacity null mode.

With `decay=true` solves the (N+1)×(N+1) system [V e; dt' 0][σ; λ] = [g; 0]
enforcing ∫σ ds = 0.  Returns (σ, λ) where λ is the additive constant.
"""
function V⁻¹(g, Γ::Boundary; decay=true)
    N  = length(Γ)
    e  = ones(N)
    dt = abs.(Γ.γt) .* Γ.w
    if decay
        A   = [Γ.V  e; dt' 0.0]
        sol = A \ [g; 0.0]
        return sol[1:N], sol[N+1]
    else
        return Γ.V \ g, 0.0
    end
end

"Apply the pre-built adjoint double-layer matrix: (K′σ)[i] = ∫_Γ ∂ₙG(γᵢ, γ') σ(γ') ds'."
function V′(σ, Γ::Boundary)
    return Γ.V′ * σ
end

"Double layer"
function K′(σ, Γ::Boundary)
    return transpose(Γ.V') * σ
end

# Workhorse: ∑ⱼ kernel(ζ - γⱼ) σⱼ |γ′ⱼ| wⱼ for scalar or array targets.
function _smooth_evaluation(ζ::Number, σ, kernel, Γ::Boundary)
    δz = ζ .- Γ.γ
    return sum(kernel.(δz) .* σ .* abs.(Γ.γt) .* Γ.w)
end

function _smooth_evaluation(ζ, σ, kernel, Γ::Boundary)
    shp = size(ζ)
    δz  = vec(ζ) .- transpose(Γ.γ)
    A   = kernel.(δz) .* transpose(σ .* abs.(Γ.γt) .* Γ.w)
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

    if length(ζ) * Nsrc <= Inf
        bx = _smooth_evaluation(ζ, σs[1], r -> real(r)/(4π*abs(r)), Γs[1])
        by = _smooth_evaluation(ζ, σs[1], r -> imag(r)/(4π*abs(r)), Γs[1])
        for k in 2:length(Γs)
            bx = bx .+ _smooth_evaluation(ζ, σs[k], r -> real(r)/(4π*abs(r)), Γs[k])
            by = by .+ _smooth_evaluation(ζ, σs[k], r -> imag(r)/(4π*abs(r)), Γs[k])
        end
        return bx, by
    end

    q  = vcat([σs[j] .* abs.(Γs[j].γt) .* Γs[j].w for j in eachindex(Γs)]...)
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

"Evaluate a domain field f at arc parameter values θ by tracing to the boundary and FFT-interpolating."
function boundary_eval(f, Ω, s)
    vals  = trace(f, Ω)
    N     = length(vals)
    f̂     = fft(vals) ./ N
    freqs = fftfreq(N) * N
    return real.(sum(f̂[k] .* cis.(freqs[k] .* s) for k in 1:N))
end

"Boundary integral: ∫_Γ σ ds = ∑ⱼ σⱼ |γ′ⱼ| wⱼ."
function integral(σ, Γ::Boundary)
    return sum(σ .* abs.(Γ.γt) .* Γ.w)
end

function integral(μ, Ω::Domain)
    return sum(μ .* Ω.dz)
end