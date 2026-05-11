
using FFTW
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

Returns (ux, uy).  Used to compute the Marangoni slip velocity on the particle surface.
Valid when ζ is well-separated from Γ; use `∇ℬ_arc` near open arcs.
"""
function ∇ℬ(ζ, σ, Γ::Boundary)
    return (_smooth_evaluation(ζ, σ, r -> real(r)/(4π*abs(r)), Γ),
            _smooth_evaluation(ζ, σ, r -> imag(r)/(4π*abs(r)), Γ))
end

"""
Near-field accurate ∇ℬ for open arc boundaries with Bruno (1/√) endpoint densities.

Far from Γ (dist > 5π/N) uses the direct Chebyshev-Gauss sum.
Near Γ reconstructs the integrand via DCT cosine series and calls QuadGK adaptively.
Returns (ux, uy).
"""
function ∇ℬ_arc(ζ::Number, α, Γ::Boundary; rtol=1e-8)
    N    = length(Γ)
    dist = minimum(abs.(ζ .- Γ.γ))

    if dist > 5π / N
        δz  = ζ .- Γ.γ
        val = sum(@. (δz / abs(δz)) * α * abs(Γ.γt) * Γ.w) / (4π)
        return real(val), imag(val)
    end

    function dct_coeffs(v)
        raw = FFTW.r2r(collect(Float64, v), FFTW.REDFT10)
        ĉ = raw ./ N;  ĉ[1] /= 2
        return ĉ
    end
    function coseval(ĉ, φ)
        s = ĉ[1]
        for k in 1:N-1;  s += ĉ[k+1] * cos(k * φ);  end
        return s
    end

    α̂  = dct_coeffs(α)
    γ̂r = dct_coeffs(real.(Γ.γ))
    γ̂i = dct_coeffs(imag.(Γ.γ))
    speed = abs(Γ.γt[1])

    function integrand(φ)
        γi = coseval(γ̂r, φ) + im * coseval(γ̂i, φ)
        δz = ζ - γi
        return (δz / abs(δz)) * coseval(α̂, φ) * speed / (4π)
    end

    val, _ = quadgk(integrand, 0.0, π; rtol=rtol)
    return real(val), imag(val)
end

function ∇ℬ_arc(ζ, α, Γ::Boundary)
    shp     = size(ζ)
    results = [∇ℬ_arc(ζi, α, Γ) for ζi in vec(ζ)]
    return reshape(first.(results), shp), reshape(last.(results), shp)
end

"Evaluate a domain field f at arc parameter values θ by tracing to the boundary and FFT-interpolating."
function boundary_eval(f, Ω, θ)
    vals  = trace(f, Ω)
    N     = length(vals)
    f̂     = fft(vals) ./ N
    freqs = fftfreq(N) * N
    return real.(sum(f̂[k] .* cis.(freqs[k] .* θ) for k in 1:N))
end

"Quadrature weight integral: ∫_Γ σ ds = ∑ⱼ σⱼ |γ′ⱼ| wⱼ."
function integral(σ, Γ::Boundary)
    return sum(σ .* abs.(Γ.γt) .* Γ.w)
end
