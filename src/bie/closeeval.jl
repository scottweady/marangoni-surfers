
# ── Close evaluation ──────────────────────────────────────────────────────────

"Find the complex preimage t of ζ under γ via Halley's method on a truncated Chebyshev interpolant."
function _find_preimage(ζ, Γ::Boundary; maxiter=20, tol=1e-15, verbose=false)

    # Truncate to cap Clenshaw amplification outside [-1,1].
    # T_k(t) ~ cosh(k·cosh⁻¹|t|); 32 terms limits amplification to ~87 at |t|=1.013.
    Nroots = min(length(Γ), 32)
    γ̂   = _cheb_coeffs(Γ.γ;   Nmax=Nroots)
    γ̂ₜ  = _cheb_coeffs(Γ.γt;  Nmax=Nroots)
    γ̂ₜₜ = _cheb_coeffs(Γ.γtt; Nmax=Nroots)

    # Affine starting guess: inverse of chord map s(τ) = (τ - τ₀)/s₀ (Klinteberg 2020, eq. 30).
    τ₀ = (Γ.γ[end] + Γ.γ[1]) / 2
    s₀ = (Γ.γ[end] - Γ.γ[1]) / 2
    t  = (ζ - τ₀) / s₀
    Q  = _cheb_eval(γ̂, t) - ζ

    verbose && println("  iter   0:  t = ", real(t), " + ", imag(t), "i  |Q| = ", abs(Q))

    for k in 1:maxiter

        abs(Q) < tol && break
        Qt  = _cheb_eval(γ̂ₜ,  t)
        Qtt = _cheb_eval(γ̂ₜₜ, t)
        iszero(Qt) && break

        # Halley step: cubic convergence, robust near near-double roots.
        denom = 2Qt^2 - Q * Qtt
        iszero(denom) && break
        δt = 2Q * Qt / denom
        α  = 1.0

        for _ in 1:8
            abs(_cheb_eval(γ̂, t - α*δt) - ζ) < abs(Q) && break
            α /= 2
        end

        t -= α * δt
        Q  = _cheb_eval(γ̂, t) - ζ

        verbose && println("  iter ", lpad(k, 3), ":  t = ", real(t), " + ", imag(t), "i  |Q| = ", abs(Q), "  α = ", α)

    end

    isconverged = abs(Q) < tol

    return t, isconverged

end

"Outer Joukowski map Φ(z) = z + √(z²-1) choosing the branch with |Φ| ≥ 1."
function _Φ(z)
    s  = sqrt(z - 1) * sqrt(z + 1)   # avoids cancellation vs sqrt(z^2-1)
    w1 = z + s
    w2 = z - s
    return abs(w1) >= abs(w2) ? w1 : w2
end

"""
Close-evaluation of the single-layer potential via singularity swap (af Klinteberg 2020).

Splits log|γ(t)-ζ| = log|Q(t)/(t-t₀)| + log|t-t₀| where t₀ = γ⁻¹(ζ).
The first part is smooth → direct Gauss-Chebyshev.
The second part → Chebyshev expansion of g = |γ'|σ, evaluated exactly via Iₖ(t₀,k).
"""
function 𝒱_close(ζ_tgt, σ, Γ::Boundary; t_tgt=nothing)

    _Iₖ(z, k)  = k == 0 ? -0.5 * log(abs(_Φ(z))/2) : (1/(2k)) * real(_Φ(z)^(-k))

    if t_tgt === nothing
        t_tgt, _ = _find_preimage(ζ_tgt, Γ)
    end

    g = σ .* abs.(Γ.γt)
    ĝ = _cheb_coeffs(g)

    smooth   = -(1/(2π)) * sum(log.(abs.((ζ_tgt .- Γ.γ) ./ (t_tgt .- Γ.t))) .* σ .* Γ.dΓ)
    singular = sum(ĝ[k+1] * _Iₖ(t_tgt, k) for k in 0:length(ĝ)-1)

    return smooth + singular

end

"""
Close-evaluation of the adjoint double-layer potential via singularity swap.

Uses the Chebyshev-Cauchy identity (af Klinteberg 2020, image formula):
  ∫₋₁¹ Tₖ(t)/(z-t) dt/√(1-t²) = π/√(z²-1) · (z - √(z²-1))^k

Writes the kernel as Re(n̂/(ζ-γ(t))) and introduces the preimage τ = γ⁻¹(ζ):
  ∫ σ|γ'|/(ζ-γ) dt/√(1-t²) = ∫ f(t)/(τ-t) dt/√(1-t²)
where f(t) = σ(t)|γ'(t)|(τ-t)/(ζ-γ(t)) is smooth (removable singularity at t=τ).
Chebyshev-expand f and apply the identity exactly.
"""
function K′_close(ζ_tgt, n̂_tgt, σ, Γ::Boundary; t_tgt=nothing)

    if t_tgt === nothing
        t_tgt, _ = _find_preimage(ζ_tgt, Γ)
    end

    g = σ .* abs.(Γ.γt)
    f = g .* (t_tgt .- Γ.t) ./ (ζ_tgt .- Γ.γ)

    ĉ = _cheb_coeffs(f)

    Φ = _Φ(t_tgt)
    s = Φ - t_tgt          # √(τ²-1), outer branch, so |Ψ| ≤ 1
    Ψ = inv(Φ)             # τ - √(τ²-1)

    cauchy = (π / s) * sum(ĉ[k+1] * Ψ^k for k in 0:length(ĉ)-1)

    return -(1/(2π)) * dot(n̂_tgt, conj(cauchy))

end
