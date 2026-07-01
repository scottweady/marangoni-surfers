
# ── Kernel matrix assembly ────────────────────────────────────────────────────

"""
Cross single-layer matrix: sources on Γ_src evaluated at targets on Γ_tgt.

Off-diagonal block (no self-interaction); kernel is smooth so plain quadrature suffices.
A[i,j] = -1/(2π) log|γ_tgt[i] - γ_src[j]| |γ′_src[j]| w_src[j].
"""
function _laplace2d_single_layer(Γ_src::Boundary, Γ_tgt::Boundary)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(N_tgt, N_src)

    for i in 1:N_tgt, j in 1:N_src
        A[i,j] = -(1/(2π)) * log(abs(Γ_tgt.γ[i] - Γ_src.γ[j])) * Γ_src.dΓ[j]
    end

    return A

end

function _laplace2d_single_layer_smooth(t, γ, γt, w)

    N = length(t)
    I1 = zeros(Float64, N, N)
    I2 = zeros(Float64, N, N)
    I3 = zeros(Float64, N, N)

    for i in 1:N, j in 1:N

        if i == j
            I1[i, j] = 0.0
            I2[i, j] = abs(γt[i]) * log(abs(γt[j])) * w[j]
            I3[i, j] = abs(γt[i]) * ((1.0 - t[i]) * log(abs(1.0 - t[i])) + (1.0 + t[i]) * log(abs(1.0 + t[i])) - 2.0)
        else
            I1[i, j] = log(abs(γ[i] - γ[j])) * abs(γt[j]) * w[j]
            I2[i, j] = abs(γt[i]) * log(abs(γ[i] - γ[j]) / abs(t[i] - t[j])) * w[j]
            I3[i, j] = 0.0
        end
    end

    A = -1.0 / (2π) * ( (I1 - abs.(γt) .* sum(I1 ./ abs.(γt)', dims=2) .* I(N)) + sum(I2, dims=2) .* I(N) + I3)

    return A

end

function _laplace2d_single_layer_theta(θ, γ, γt, w, t̃, γ̃, γ̃t, w̃)

    N = length(θ)
    I1 = zeros(Float64, N, N)
    I2 = zeros(Float64, N, N)

    for i in 1:N, j in 1:N

        if i == j
            I1[i, j] = 0.0
        else
            I1[i, j] = log(abs(γ[i] - γ[j])) * abs(γt[j]) * w[j]
        end

        I2[i,j] = log(abs(γ[i] - γ̃[j])) * abs(γ̃t[j]) * w̃[j]

    end

    A = -1.0 / (2π) * ( (I1 - sum(I1, dims=2) .* I(N)) + sum(I2, dims=2) .* I(N) )

    return A

end

"""
Single-layer matrix for the Laplace log-kernel on a smooth closed curve.

Uses Kress splitting: A[i,j] = (R[|i-j|] + h * regKernel[i,j]) * speed[j], where
the singular part is absorbed into the Kress weights and the remainder is smooth.
Arguments are the N uniform parameter nodes `t`, curve points `γ`, and tangents `γt`.
"""
function _laplace2d_single_layer_kress(t, γ, γt)

    N = length(t)
    A = zeros(Float64, N, N)
    _, w = _kresspts(N)
    h = 2π / N

    for i in 1:N, j in 1:N
        k = mod(i - j, N) + 1
        A_reg = i == j ? -(1/(2π)) * log(abs(γt[j])) : -(1/(2π)) * log(abs(γ[i] - γ[j]) / (2 * abs(sin((t[i] - t[j]) / 2))))
        A[i, j] = (w[k] + h * A_reg) * abs(γt[j])
    end

    return A

end

"""
Single-layer matrix for the Laplace log-kernel on an open arc with 1/√(1-t²) endpoint density.

Uses Bruno (2021) quadrature: the density is split as σ(t) = α(t)/√(1-t²) where α is smooth,
and the log singularity is absorbed into the precomputed coefficients R^(N) (eq. 4.4 of Bruno 2021).
Arguments are Gauss-Chebyshev nodes `t` ∈ (-1,1), curve points `γ`, tangents `γt`, and weights `w`.
"""
function _laplace2d_single_layer_bruno(t, γ, γt, w)

    N = length(t)

    # R^(N)(l) = -∑_{m=0}^{N-1} (2-δ_{0m}) λ_m cos(lmπ/N), Bruno 2021 eq. 4.4
    R = zeros(2N)

    for l in 0:2N-1
        for m in 0:N-1
            coeff  = m == 0 ? 1.0 : 2.0
            λm     = m == 0 ? log(2) / 2 : 1 / (2m)
            R[l+1] -= coeff * λm * cos(l * m * π / N)
        end
    end

    A = zeros(Float64, N, N)

    for i in 1:N, j in 1:N
      Rij = R[abs(i - j) + 1] + R[(i + j)]
      E2  = i == j ? log(abs(γt[j])) : log(abs(γ[i] - γ[j]) / abs(t[i] - t[j]))
      A[i, j] = -(1/(2π)) * (Rij + E2) * abs(γt[j]) * w[j]
    end

    return A

end

"""
Single-layer matrix for the Helmholtz K₀-kernel on a smooth closed curve (Kress quadrature).

Splits K₀(r)/(2π) = -log(r)/(2π) + [K₀(r)+log(r)]/(2π).
The log part uses the same Kress weights as the Laplace matrix.
The smooth correction [K₀(r)+log(r)]/(2π) → (log 2 - γ_E)/(2π) as r → 0 is integrated
by the trapezoidal rule (h * smooth_value * |γ′ⱼ|).
"""
function _helmholtz2d_single_layer_kress(t, γ, γt; k=1.0)

    N = length(t)
    H = _laplace2d_single_layer_kress(t, γ, γt)  # Laplace SL matrix (handles log singularity)
    h = 2π / N
    γ_E = Base.MathConstants.eulergamma

    for i in 1:N, j in 1:N
        if i == j
            # K₀(kr) + log(r) → log(2) - γ_E - log(k) as r → 0
            H[i, j] += h / (2π) * (log(2) - γ_E - log(k)) * abs(γt[j])
        else
            r = abs(γ[i] - γ[j])
            H[i, j] += h / (2π) * (besselk(0, k*r) + log(r)) * abs(γt[j])
        end
    end

    return H

end

"""
Adjoint double-layer (K′) matrix for the Laplace kernel on a single boundary.

K′[i,j] = -1/(2π) (n(γᵢ)·(γᵢ-γⱼ)) / |γᵢ-γⱼ|² |γ′ⱼ| wⱼ  (i ≠ j)
Diagonal entry uses the curvature limit.  Normal convention: outward from the bounded domain.
"""
function _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)

    N = length(t)
    A = zeros(Float64, N, N)
    n̂ = (imag.(γt) - im * real.(γt)) ./ abs.(γt)

    for i in 1:N, j in 1:N

      if i == j
        A[i, j] = (1/(4π)) * dot(n̂[i], γtt[i]) / abs(γt[i]) * w[i]
      else
        δz = γ[i] - γ[j]
        A[i, j] = -(1/(2π)) * (dot(n̂[i], δz) / abs2(δz)) * abs(γt[j]) * w[j]
      end

    end

    return A

end

"""
Adjoint double-layer (K_k′) matrix for the Helmholtz kernel on a smooth closed curve.

K_k′[i,j] = -k/(2π) K₁(k|γᵢ-γⱼ|) (n(γᵢ)·(γᵢ-γⱼ))/|γᵢ-γⱼ| |γ′ⱼ| wⱼ  (i ≠ j)
Diagonal uses the curvature limit  κ(γᵢ)/(4π).
The kernel has a removable singularity (n·r = O(|r|²), K₁(kr)/r ~ 1/(kr²)), so
plain trapezoidal quadrature is spectrally accurate.
"""
function _helmholtz2d_double_layer_adjoint(t, γ, γt, γtt, w; k=1.0)

    N = length(t)
    A = zeros(Float64, N, N)
    n̂ = (imag.(γt) .- im .* real.(γt)) ./ abs.(γt)

    for i in 1:N, j in 1:N

        if i == j
            A[i, j] = (1/(4π)) * real(conj(n̂[i]) * γtt[i]) / abs(γt[i]) * w[i]
        else
            δz = γ[i] - γ[j]
            r  = abs(δz)
            A[i, j] = -(k/(2π)) * besselk(1, k*r) * real(conj(n̂[i]) * δz) / r * abs(γt[j]) * w[j]
        end

    end

    return A

end

"""
Self double-layer matrix K on a closed curve: normal at SOURCE j.

K[i,j] = +1/(2π) (n̂(γⱼ)·(γᵢ-γⱼ)) / |γᵢ-γⱼ|² |γ′ⱼ| wⱼ  (i ≠ j)
Diagonal uses the same curvature limit as K′ (the limit is symmetric at coincident points).
The kernel has a removable singularity so the trapezoidal rule is spectrally accurate.
"""
function _laplace2d_double_layer(t, γ, γt, γtt, w)

    N = length(t)
    A = zeros(Float64, N, N)
    n̂ = (imag.(γt) .- im .* real.(γt)) ./ abs.(γt)

    for i in 1:N, j in 1:N
        if i == j
            A[i, j] = (1/(4π)) * real(conj(n̂[j]) * γtt[j]) / abs(γt[j]) * w[j]
        else
            δz = γ[i] - γ[j]
            A[i, j] = (1/(2π)) * real(conj(n̂[j]) / conj(δz)) * abs(γt[j]) * w[j]
        end
    end

    return A

end

"""
Cross adjoint double-layer matrix: sources on Γ_src, normal direction from Γ_tgt.

Normal convention: outward from the bounded domain of Γ_tgt.
A[i,j] = -1/(2π) (n_tgt[i]·(γ_tgt[i]-γ_src[j])) / |γ_tgt[i]-γ_src[j]|² |γ′_src[j]| w_src[j].
"""
function _laplace2d_double_layer_adjoint(Γ_src::Boundary, Γ_tgt::Boundary)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(Float64, N_tgt, N_src)

    for i in 1:N_tgt

        n̂ = (imag(Γ_tgt.γt[i]) - im * real(Γ_tgt.γt[i])) / abs(Γ_tgt.γt[i])

        for j in 1:N_src
            δz = Γ_tgt.γ[i] - Γ_src.γ[j]
            A[i, j] = -(1/(2π)) * (dot(n̂, δz) / abs2(δz)) * Γ_src.dΓ[j]
        end

    end

    return A

end

function _periodic_diff_matrix(M)

    m = [0:M÷2; -(M÷2-1):-1]
    ik = im .* Float64.(m)

    D = real(ifft(Diagonal(ik) * fft(I(M),1), 1))

    return D

end

function _laplace2d_hypersingular(t, γ, γt, γtt)

    N  = length(t)
    w  = fill(2π/N, N)
    dΓ = abs.(γt) .* w

    S  = _laplace2d_single_layer_kress(t, γ, γt)
    Dt = _periodic_diff_matrix(N)
    Ds = Diagonal(1.0 ./ abs.(γt)) * Dt
    return -Diagonal(1.0 ./ dΓ) * Ds' * Diagonal(dΓ) * S * Ds

end
