
# ── Close evaluation ──────────────────────────────────────────────────────────

"""
    _find_preimage(ζ, Γ; maxiter=20, tol=1e-15, verbose=false, t_init=nothing)

Find the complex preimage `t` of `ζ` under `γ` via Halley's method on a truncated
Chebyshev interpolant.

Truncates the interpolant to `min(length(Γ), 32)` Chebyshev modes to cap Clenshaw
amplification outside `[-1,1]` (`T_k(t) ~ cosh(k·cosh⁻¹|t|)`; 32 terms limits
amplification to ~87 at `|t|=1.013`), and by default starts from the affine chord-map
inverse (Klinteberg 2020, eq. 30). For strongly curved arcs (e.g. a half ellipse) the
chord-map guess can converge to a spurious preimage; pass `t_init` (e.g. the parameter
of the nearest node) to avoid this.

# Arguments
- `ζ`: target point in the complex plane.
- `Γ::Boundary`: boundary whose parametrization is being inverted.
- `maxiter`: maximum Halley iterations (default `20`).
- `tol`: convergence tolerance on `|γ(t) - ζ|` (default `1e-15`).
- `verbose`: print iteration diagnostics if `true` (default `false`).
- `t_init`: initial guess for `t` (default `nothing`, use the chord-map inverse).

# Returns
- `(t, isconverged)`: the preimage and whether `|γ(t) - ζ| < tol` was reached.
"""
function _find_preimage(ζ, Γ::Boundary; maxiter=20, tol=1e-15, verbose=false, t_init=nothing)

    Nroots = min(length(Γ), 32)
    γ̂   = _cheb_coeffs(Γ.γ;   Nmax=Nroots)
    γ̂ₜ  = _cheb_coeffs(Γ.γt;  Nmax=Nroots)
    γ̂ₜₜ = _cheb_coeffs(Γ.γtt; Nmax=Nroots)

    τ₀ = (Γ.γ[end] + Γ.γ[1]) / 2
    s₀ = (Γ.γ[end] - Γ.γ[1]) / 2
    t  = t_init === nothing ? (ζ - τ₀) / s₀ : complex(t_init)
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

"""
    _Φ(z)

Outer Joukowski map `Φ(z) = z + √(z²-1)`, choosing the branch with `|Φ| ≥ 1`.

# Arguments
- `z`: (possibly complex) evaluation point.

# Returns
- `Φ(z)` on the outer branch.
"""
function _Φ(z)
    s  = sqrt(z - 1) * sqrt(z + 1)   # avoids cancellation vs sqrt(z^2-1)
    w1 = z + s
    w2 = z - s
    return abs(w1) >= abs(w2) ? w1 : w2
end

"""
    𝒱_close(ζ_tgt, σ, Γ; t_tgt=nothing)

Close evaluation of the single-layer potential via singularity swap (af Klinteberg 2020).

Splits `log|γ(t)-ζ| = log|Q(t)/(t-t₀)| + log|t-t₀|` where `t₀ = γ⁻¹(ζ)`. The first part
is smooth and evaluated by direct Gauss-Chebyshev quadrature; the second part is the
Chebyshev expansion of `g = |γ'|σ`, evaluated exactly via `Iₖ(t₀,k)`.

# Arguments
- `ζ_tgt`: target point, close to or on `Γ`.
- `σ`: density sampled at the quadrature nodes of `Γ`.
- `Γ::Boundary`: source boundary.
- `t_tgt`: precomputed preimage of `ζ_tgt` under `Γ.γ`, from [`_find_preimage`](@ref)
  (computed internally if `nothing`).

# Returns
- The single-layer potential at `ζ_tgt`.
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
    K′_close(ζ_tgt, n̂_tgt, σ, Γ; t_tgt=nothing)

Close evaluation of the adjoint double-layer potential via singularity swap.

Uses the Chebyshev-Cauchy identity (af Klinteberg 2020, image formula)
`∫₋₁¹ Tₖ(t)/(z-t) dt/√(1-t²) = π/√(z²-1) · (z - √(z²-1))^k`.
Writes the kernel as `Re(n̂/(ζ-γ(t)))` and introduces the preimage `τ = γ⁻¹(ζ)` so that
`∫ σ|γ'|/(ζ-γ) dt/√(1-t²) = ∫ f(t)/(τ-t) dt/√(1-t²)`, where
`f(t) = σ(t)|γ'(t)|(τ-t)/(ζ-γ(t))` is smooth (removable singularity at `t=τ`).
Chebyshev-expands `f` and applies the identity exactly.

# Arguments
- `ζ_tgt`: target point, close to or on `Γ`.
- `n̂_tgt`: unit normal at the target, used to project the (vector) kernel onto a scalar.
- `σ`: density sampled at the quadrature nodes of `Γ`.
- `Γ::Boundary`: source boundary.
- `t_tgt`: precomputed preimage of `ζ_tgt` under `Γ.γ`, from [`_find_preimage`](@ref)
  (computed internally if `nothing`).

# Returns
- The adjoint double-layer potential at `ζ_tgt`.
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
