
using LinearAlgebra

"Gauss-Legendre nodes and weights on `dom`, scaled from [-1,1]."
function legpts(n; dom=[-1,1])
    x, w = gausslegendre(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

"Uniform trigonometric nodes and weights on `dom` (periodic trapezoid rule)."
function trigpts(n; dom=[0, 2π])
    a, b = dom
    x = range(a, stop=b, length=n+1)[1:end-1]
    w = fill((b - a) / n, n)
    return x, w
end

"""
Kress singular quadrature weights for the log-kernel on a closed curve.

Returns the N-vector R^(N) such that the (i,j) entry of the single-layer matrix is
`(R[mod(i-j,N)+1] + h * regularKernel[i,j]) * speed[j]`.
Reference: Kress, Linear Integral Equations, §12.
"""
function kressweights(N::Int)
    @assert iseven(N) "Kress quadrature requires even N."
    θ = 2π .* (0:N-1) ./ N
    w = zeros(Float64, N)
    for j in 0:N-1
        s = 0.0
        for m in 1:(N ÷ 2 - 1)
            s += cos(m * θ[j+1]) / m
        end
        s += cos((N ÷ 2) * θ[j+1]) / N
        w[j+1] = s / N
    end
    return w
end

"""
Single-layer matrix for the Laplace log-kernel on a smooth closed curve.

Uses Kress splitting: A[i,j] = (R[|i-j|] + h * regKernel[i,j]) * speed[j], where
the singular part is absorbed into the Kress weights and the remainder is smooth.
Arguments are the N uniform parameter nodes `s`, curve points `γ`, and tangents `γs`.
"""
function _laplace2d_single_layer_kress(s, γ, γs)
    N     = length(s)
    wsing = kressweights(N)
    h     = 2π / N
    A     = zeros(Float64, N, N)
    for i in 1:N
        for j in 1:N
            k     = mod(i - j, N) + 1
            speed = abs(γs[j])
            if i == j
                regkernel = -(1/(2π)) * log(speed)
            else
                denom     = 2 * abs(sin((s[i] - s[j]) / 2))
                regkernel = -(1/(2π)) * log(abs(γ[i] - γ[j]) / denom)
            end
            A[i, j] = (wsing[k] + h * regkernel) * speed
        end
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
    for i in 1:N
        for j in 1:N
            Rij = R[abs(i - j) + 1] + R[(i + j)]
            E2  = i == j ? log(abs(γt[j])) : log(abs(γ[i] - γ[j]) / abs(t[i] - t[j]))
            A[i, j] = -(1/(2π)) * (Rij + E2) * abs(γt[j]) * w[j]
        end
    end
    return A
end

"""
Adjoint double-layer (K′) matrix for the Laplace kernel on a single boundary.

K′[i,j] = -1/(2π) (n(γᵢ)·(γᵢ-γⱼ)) / |γᵢ-γⱼ|² |γ′ⱼ| wⱼ  (i ≠ j)
Diagonal entry uses the curvature limit.  Normal convention: outward from the bounded domain.
"""
function _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    N = length(t)
    A = zeros(Float64, N, N)
    n = (imag.(γt) - im * real.(γt)) ./ abs.(γt)
    for i in 1:N
        for j in 1:N
            if i == j
                A[i, j] = (1/(4π)) * dot(n[i], γtt[i]) / abs(γt[i])^3 * abs(γt[i]) * w[i]
            else
                δz       = γ[i] - γ[j]
                A[i, j]  = -(1/(2π)) * (dot(n[i], δz) / abs2(δz)) * abs(γt[j]) * w[j]
            end
        end
    end
    return A
end

"""
Cross single-layer matrix: sources on Γ_src evaluated at targets on Γ_tgt.

Off-diagonal block (no self-interaction); kernel is smooth so plain quadrature suffices.
A[i,j] = -1/(2π) log|γ_tgt[i] - γ_src[j]| |γ′_src[j]| w_src[j].
"""
function _laplace2d_single_layer(Γ_src::Boundary, Γ_tgt::Boundary)
    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(Float64, N_tgt, N_src)
    for i in 1:N_tgt, j in 1:N_src
        A[i,j] = -(1/(2π)) * log(abs(Γ_tgt.γ[i] - Γ_src.γ[j])) * abs(Γ_src.γt[j]) * Γ_src.w[j]
    end
    return A
end

"""
Cross adjoint double-layer matrix: sources on Γ_src, normal direction from Γ_tgt.

Off-diagonal block; kernel is smooth so plain quadrature suffices.
Normal convention: outward from the bounded domain of Γ_tgt.
A[i,j] = -1/(2π) (n_tgt[i]·(γ_tgt[i]-γ_src[j])) / |γ_tgt[i]-γ_src[j]|² |γ′_src[j]| w_src[j].
"""
function _laplace2d_double_layer_adjoint(Γ_src::Boundary, Γ_tgt::Boundary)
    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(Float64, N_tgt, N_src)
    for i in 1:N_tgt
        n = (imag(Γ_tgt.γt[i]) - im * real(Γ_tgt.γt[i])) / abs(Γ_tgt.γt[i])
        for j in 1:N_src
            δz      = Γ_tgt.γ[i] - Γ_src.γ[j]
            A[i,j]  = -(1/(2π)) * (dot(n, δz) / abs2(δz)) * abs(Γ_src.γt[j]) * Γ_src.w[j]
        end
    end
    return A
end
