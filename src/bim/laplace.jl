
using LinearAlgebra
using FFTW

include("closeeval.jl")

"Gauss-Legendre nodes and weights on `dom`, scaled from [-1,1]."
function _legpts(n; dom=[-1,1])
    x, w = gausslegendre(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

"Uniform trigonometric nodes and weights on `dom` (periodic trapezoid rule)."
function _trigpts(n; dom=[0, 2π])
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
function _kressweights(N::Int)

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
Cross single-layer matrix: sources on Γ_src evaluated at targets on Γ_tgt.

Off-diagonal block (no self-interaction); kernel is smooth so plain quadrature suffices.
A[i,j] = -1/(2π) log|γ_tgt[i] - γ_src[j]| |γ′_src[j]| w_src[j].
"""
function _laplace2d_single_layer(Γ_src::Boundary, Γ_tgt::Boundary; close_eval=false, close_eval_tol=0.01)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(N_tgt, N_src)

    for i in 1:N_tgt

        if minimum(abs.(Γ_tgt.γ[i] .- Γ_src.γ)) < close_eval_tol && close_eval

            t_tgt, isconverged = _find_preimage(Γ_tgt.γ[i], Γ_src)
            isclose = true

            if !isconverged
                @warn("Preimage did not converge for target $i; using direct quadrature.")
            end

        else
            isclose = false
        end

        for j in 1:N_src
            if isclose && isconverged
                ej = zeros(N_src); ej[j] = 1.0
                A[i, j] = 𝒱_close(Γ_tgt.γ[i], ej, Γ_src; t_tgt=t_tgt)
            else
                A[i,j] = -(1/(2π)) * log(abs(Γ_tgt.γ[i] - Γ_src.γ[j])) * Γ_src.dΓ[j]
            end
        end

    end

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
    w = _kressweights(N)
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
Cross adjoint double-layer matrix: sources on Γ_src, normal direction from Γ_tgt.

Uses close-evaluation (singularity swap) for targets within `tol` of the source.
Normal convention: outward from the bounded domain of Γ_tgt.
A[i,j] = -1/(2π) (n_tgt[i]·(γ_tgt[i]-γ_src[j])) / |γ_tgt[i]-γ_src[j]|² |γ′_src[j]| w_src[j].
"""
function _laplace2d_double_layer_adjoint(Γ_src::Boundary, Γ_tgt::Boundary; close_eval=false, close_eval_tol=0.01)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(Float64, N_tgt, N_src)

    for i in 1:N_tgt

        n̂ = (imag(Γ_tgt.γt[i]) - im * real(Γ_tgt.γt[i])) / abs(Γ_tgt.γt[i])

        if close_eval && minimum(abs.(Γ_tgt.γ[i] .- Γ_src.γ)) < close_eval_tol
            t_tgt, isconverged = _find_preimage(Γ_tgt.γ[i], Γ_src)
            isclose = isconverged

            if !isconverged
                @warn("Preimage did not converge for target $i; using direct quadrature.")
            end

        else
            isclose = false
        end

        for j in 1:N_src
            if isclose && isconverged
                ej = zeros(N_src); ej[j] = 1.0
                A[i, j] = K′_close(Γ_tgt.γ[i], n̂, ej, Γ_src; t_tgt=t_tgt)
            else
                δz = Γ_tgt.γ[i] - Γ_src.γ[j]
                A[i, j] = -(1/(2π)) * (dot(n̂, δz) / abs2(δz)) * Γ_src.dΓ[j]
            end
        end

    end

    return A

end
