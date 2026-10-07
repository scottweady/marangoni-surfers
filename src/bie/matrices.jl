
# -- Kernel matrix assembly ----------------------------------------------------

"""
    _laplace2d_single_layer_cross(Γ_src, Γ_tgt)

Cross single-layer matrix: sources on `Γ_src` evaluated at targets on `Γ_tgt`.

Off-diagonal block (no self-interaction); the kernel is smooth so plain quadrature
suffices.

# Arguments
- `Γ_src::Boundary`: source boundary.
- `Γ_tgt::Boundary`: target boundary.

# Returns
- `A::Matrix{Float64}` (`N_tgt × N_src`) with
  `A[i,j] = -1/(2π) log|γ_tgt[i] - γ_src[j]| |γ′_src[j]| w_src[j]`.
"""
function _laplace2d_single_layer_cross(Γ_src::Boundary, Γ_tgt::Boundary)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(N_tgt, N_src)

    for i in 1:N_tgt, j in 1:N_src
        A[i,j] = -(1/(2π)) * log(abs(Γ_tgt.γ[i] - Γ_src.γ[j])) * Γ_src.dΓ[j]
    end

    return A

end

"""
    _laplace2d_single_layer_smooth(t, γ, γt, w)

Single-layer matrix for the Laplace log-kernel on a smooth open arc with Fejér nodes, by
product integration.

The kernel is split as `log|γ(t) - γ(t′)| = log|t - t′| + log(|γ(t) - γ(t′)|/|t - t′|)`. The
smooth remainder is integrated with the Fejér weights. The `log|t - t′|` part is applied
exactly to the Chebyshev interpolant of the integrand, using the moments from
[`_log_moments`](@ref) and the transform from [`_cheb_matrix`](@ref).

# Arguments
- `t`: `N` Fejér nodes (the Gauss-Chebyshev points) in `(-1,1)`.
- `γ`, `γt`: curve points and tangents at `t`.
- `w`: Fejér quadrature weights.

# Returns
- `A::Matrix{Float64}` (`N × N`) self single-layer matrix.
"""
function _laplace2d_single_layer_smooth(t, γ, γt, w)

    N = length(t)

    arg = (γ .- transpose(γ)) ./ (t .- transpose(t))
    arg[diagind(arg)] .= γt
    A = -(1/(2π)) * log.(abs.(arg)) .* transpose(abs.(γt) .* w)

    _, T′ = _cheb_matrix(N)
    Mₖ = reduce(vcat, [transpose(_log_moments(t[i], N)) for i in 1:N])
    A .+= -(1/(2π)) * (Mₖ * T′) .* transpose(abs.(γt))

    return A

end

"""
    _laplace2d_single_layer_closed(t, γ, γt)

Single-layer matrix for the Laplace log-kernel on a smooth closed curve.

Uses Kress splitting: `A[i,j] = (R[|i-j|] + h * regKernel[i,j]) * speed[j]`, where the
singular part is absorbed into the Kress weights and the remainder is smooth.

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`: curve points and tangents at `t`.

# Returns
- `A::Matrix{Float64}` (`N × N`) self single-layer matrix.
"""
function _laplace2d_single_layer_closed(t, γ, γt)

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
    _laplace2d_single_layer_weighted(t, γ, γt, w)

Single-layer matrix for the Laplace log-kernel on an open arc with `1/√(1-t²)`-weighted
Gauss-Chebyshev nodes, via Bruno (2021) log-splitting.

The smooth remainder after subtracting the `log|t-t′|` singularity is applied exactly
in the Chebyshev basis using [`_cheb_matrix`](@ref) and [`_single_layer_eig_weighted`](@ref).

# Arguments
- `t`: `N` Gauss-Chebyshev nodes.
- `γ`, `γt`: curve points and tangents at `t`.
- `w`: Gauss-Chebyshev quadrature weights.

# Returns
- `A::Matrix{Float64}` (`N × N`) self single-layer matrix.
"""
function _laplace2d_single_layer_weighted(t, γ, γt, w)

    M = length(t)

    arg = (γ .- transpose(γ)) ./ (t .- transpose(t))
    arg[diagind(arg)] .= γt
    A = -(1/(2π)) * log.(abs.(arg)) .* transpose(abs.(γt) .* w)

    T, T′ = _cheb_matrix(M)
    Λ = _single_layer_eig_weighted(M)
    A .+= (T * Λ * T′) .* transpose(abs.(γt))

    return A

end

"""
    _helmholtz2d_single_layer_closed(t, γ, γt; k=1.0)

Single-layer matrix for the Helmholtz `K₀`-kernel on a smooth closed curve (Kress quadrature).

Splits `K₀(r)/(2π) = -log(r)/(2π) + [K₀(r)+log(r)]/(2π)`. The log part uses the same
Kress weights as the Laplace matrix. The smooth correction `[K₀(r)+log(r)]/(2π) →
(log 2 - γ_E)/(2π)` as `r → 0` is integrated by the trapezoidal rule (`h * smooth_value *
|γ′ⱼ|`).

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`: curve points and tangents at `t`.
- `k`: Helmholtz wavenumber (default `1.0`).

# Returns
- `H::Matrix{Float64}` (`N × N`) self single-layer matrix.
"""
function _helmholtz2d_single_layer_closed(t, γ, γt; k=1.0)

    N = length(t)
    H = _laplace2d_single_layer_closed(t, γ, γt)  # Laplace SL matrix (handles log singularity)
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
    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)

Adjoint double-layer (`K′`) matrix for the Laplace kernel on a single closed boundary.

Diagonal entry uses the curvature limit. Normal convention: outward from the bounded
domain.

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`, `γtt`: curve points and first/second derivatives at `t`.
- `w`: quadrature weights.

# Returns
- `A::Matrix{Float64}` (`N × N`) with, for `i ≠ j`,
  `A[i,j] = -1/(2π) (n(γᵢ)·(γᵢ-γⱼ)) / |γᵢ-γⱼ|² |γ′ⱼ| wⱼ`.
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
    _helmholtz2d_double_layer_adjoint(t, γ, γt, γtt, w; k=1.0)

Adjoint double-layer (`K_k′`) matrix for the Helmholtz kernel on a smooth closed curve.

Diagonal uses the curvature limit `κ(γᵢ)/(4π)`. The kernel has a removable singularity
(`n·r = O(|r|²)`, `K₁(kr)/r ~ 1/(kr²)`), so plain trapezoidal quadrature is spectrally
accurate.

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`, `γtt`: curve points and first/second derivatives at `t`.
- `w`: quadrature weights.
- `k`: Helmholtz wavenumber (default `1.0`).

# Returns
- `A::Matrix{Float64}` (`N × N`) with, for `i ≠ j`,
  `A[i,j] = -k/(2π) K₁(k|γᵢ-γⱼ|) (n(γᵢ)·(γᵢ-γⱼ))/|γᵢ-γⱼ| |γ′ⱼ| wⱼ`.
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
    _laplace2d_double_layer(t, γ, γt, γtt, w)

Self double-layer matrix `K` on a closed curve: normal taken at the source point.

Diagonal uses the same curvature limit as [`_laplace2d_double_layer_adjoint`](@ref) (the
limit is symmetric at coincident points). The kernel has a removable singularity so the
trapezoidal rule is spectrally accurate.

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`, `γtt`: curve points and first/second derivatives at `t`.
- `w`: quadrature weights.

# Returns
- `A::Matrix{Float64}` (`N × N`) with, for `i ≠ j`,
  `A[i,j] = +1/(2π) (n̂(γⱼ)·(γᵢ-γⱼ)) / |γᵢ-γⱼ|² |γ′ⱼ| wⱼ`.
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
    _laplace2d_double_layer_adjoint(Γ_src, Γ_tgt)

Cross adjoint double-layer matrix: sources on `Γ_src`, normal direction taken from `Γ_tgt`.

Normal convention: outward from the bounded domain of `Γ_tgt`.

# Arguments
- `Γ_src::Boundary`: source boundary.
- `Γ_tgt::Boundary`: target boundary (supplies the normal direction).

# Returns
- `A::Matrix{Float64}` (`N_tgt × N_src`) with
  `A[i,j] = -1/(2π) (n_tgt[i]·(γ_tgt[i]-γ_src[j])) / |γ_tgt[i]-γ_src[j]|² |γ′_src[j]| w_src[j]`.
"""
function _laplace2d_double_layer_adjoint(Γ_src::Boundary, Γ_tgt::Boundary)

    N_src, N_tgt = length(Γ_src), length(Γ_tgt)
    A = zeros(N_tgt, N_src)

    for i in 1:N_tgt

        n̂ = (imag(Γ_tgt.γt[i]) - im * real(Γ_tgt.γt[i])) / abs(Γ_tgt.γt[i])

        for j in 1:N_src
            δz = Γ_tgt.γ[i] - Γ_src.γ[j]
            A[i, j] = -(1/(2π)) * (dot(n̂, δz) / abs2(δz)) * Γ_src.dΓ[j]
        end

    end

    return A

end

"""
    _periodic_diff_matrix(M)

Spectral differentiation matrix for `M`-point periodic (trapezoidal) data, via the FFT.

# Arguments
- `M`: number of periodic grid points.

# Returns
- `D::Matrix{Float64}` (`M × M`) such that `D * f` approximates `df/dθ` at the grid points.
"""
function _periodic_diff_matrix(M)

    m = [0:M÷2; -(M÷2-1):-1]
    ik = im .* Float64.(m)

    D = real(ifft(Diagonal(ik) * fft(I(M),1), 1))

    return D

end

"""
    _laplace2d_hypersingular(t, γ, γt, γtt)

Hypersingular (`T`) operator matrix for the Laplace equation on a closed curve, built as
`-D_s' V D_s` acting on the single-layer matrix, where `D_s` is arc-length differentiation.

# Arguments
- `t`: `N` uniform parameter nodes.
- `γ`, `γt`, `γtt`: curve points and first/second derivatives at `t`.

# Returns
- `Matrix{Float64}` (`N × N`) hypersingular matrix.
"""
function _laplace2d_hypersingular(t, γ, γt, γtt)

    N  = length(t)
    w  = fill(2π/N, N)
    dΓ = abs.(γt) .* w

    S  = _laplace2d_single_layer_closed(t, γ, γt)
    Dt = _periodic_diff_matrix(N)
    Ds = Diagonal(1.0 ./ abs.(γt)) * Dt
    return -Diagonal(1.0 ./ dΓ) * Ds' * Diagonal(dΓ) * S * Ds

end

# -- Elliptical arcs parametrized by φ ------------------------------------------
#
# An arc s ∈ [sa, sb] of the ellipse γ(s) = a cos s + i b sin s is parametrized by
# s = sa + J(cos φ + 1), J = (sb - sa)/2, φ = π(1 - ψ)/2, with Fejér nodes in ψ. Densities
# with half-integer powers of the distance to the arc ends are smooth in φ.

_ellipse_speed(a, b, s) = abs(-a * sin(s) + im * b * cos(s))

# log|γ(s) - γ(s′)| from Δs = s - s′ (computed accurately by the caller) and m = (s + s′)/2,
# using γ(s) - γ(s′) = 2 sin(Δs/2) (-a sin m + i b cos m)
_ellipse_logdist(a, b, Δs, m) = log(2 * abs(sin(Δs / 2))) + log(abs(-a * sin(m) + im * b * cos(m)))

"""
    _laplace2d_single_layer_phi(a, b, sa, sb, ψ, w)

Self single-layer matrix on a φ-parametrized elliptical arc, acting on the physical density.

With `t = cos φ`, `log|cos φ - cos φ′|` has logarithmic singularities at `φ′ = φ`, `-φ` and
`2π - φ`, i.e. at `ψ′ = ψ`, `2 - ψ` and `-2 - ψ`. These are integrated exactly against the
Chebyshev interpolant of the integrand `σ |γ_ψ|`, using [`_log_moments`](@ref) for the
interior point and [`_log_moments_ext`](@ref) for the two images. The smooth remainder is
integrated with the Fejér weights.

# Arguments
- `a`, `b`: semi-axes of the ellipse.
- `sa`, `sb`: parameter interval of the arc.
- `ψ`, `w`: Fejér nodes and weights.

# Returns
- `A::Matrix{Float64}` (`N × N`) self single-layer matrix.
"""
function _laplace2d_single_layer_phi(a, b, sa, sb, ψ, w)

    N = length(ψ)
    J = (sb - sa) / 2
    φ = π .* (1 .- ψ) ./ 2
    s = sa .+ J .* (cos.(φ) .+ 1)
    γψ = _ellipse_speed.(a, b, s) .* J .* sin.(φ) .* (π / 2)    # |dγ/dψ|
    _, T′ = _cheb_matrix(N)

    A = zeros(N, N)
    Mom = zeros(N, N)
    for i in 1:N
        x1, x2 = -2 - ψ[i], 2 - ψ[i]
        Mom[i, :] = _log_moments(ψ[i], N) .+ _log_moments_ext(x1, N) .+ _log_moments_ext(x2, N)
        for j in 1:N
            images = log(abs(ψ[j] - x1)) + log(abs(ψ[j] - x2))
            if i == j
                R = log(γψ[i]) - images
            else
                Δs = -2J * sin((φ[i] + φ[j]) / 2) * sin((φ[i] - φ[j]) / 2)
                R = _ellipse_logdist(a, b, Δs, (s[i] + s[j]) / 2) - log(abs(ψ[j] - ψ[i])) - images
            end
            A[i, j] = R * w[j]
        end
    end
    A .+= Mom * T′

    return -(1/(2π)) .* A .* transpose(γψ)

end

"""
    _laplace2d_single_layer_phi_cross(a, b, sa, sb, ψ, w, s_tgt)

Single layer of a φ-parametrized elliptical arc, acting on its physical density, at
targets elsewhere on the same ellipse.

For a target with arc coordinate `c = cos φ*` outside `[-1, 1]`, `log|c - cos φ′|` is
singular at the complex conjugate points `ψ′ = 1 - 2φ*/π` and its conjugate, which lie close
to the arc end when the target is near it. These are integrated exactly against the
Chebyshev interpolant using [`_log_moments_ext`](@ref); the remainder uses Fejér weights.

# Arguments
- `a`, `b`: semi-axes of the ellipse.
- `sa`, `sb`: parameter interval of the source arc.
- `ψ`, `w`: Fejér nodes and weights of the source arc.
- `s_tgt`: ellipse parameters of the targets, off the source arc.

# Returns
- `A::Matrix{Float64}` (`N_tgt × N`).
"""
function _laplace2d_single_layer_phi_cross(a, b, sa, sb, ψ, w, s_tgt)

    N = length(ψ)
    J = (sb - sa) / 2
    φ = π .* (1 .- ψ) ./ 2
    s = sa .+ J .* (cos.(φ) .+ 1)
    γψ = _ellipse_speed.(a, b, s) .* J .* sin.(φ) .* (π / 2)
    _, T′ = _cheb_matrix(N)

    A = zeros(length(s_tgt), N)
    Mom = zeros(length(s_tgt), N)
    for (i, sx) in enumerate(s_tgt)
        # branch of the target parameter adjacent to the arc, split at the opposite midpoint
        while sx > sa + J + π;  sx -= 2π; end
        while sx <= sa + J - π; sx += 2π; end
        c = (sx - sa) / J - 1
        z = 1 - 2 * acos(complex(c)) / π
        Mom[i, :] = _log_moments_ext(z, N) .+ _log_moments_ext(conj(z), N)
        for j in 1:N
            # c - cos φ′, written without cancellation near the arc ends
            cm = c > 0 ? (c - 1) + 2sin(φ[j] / 2)^2 : (c + 1) - 2cos(φ[j] / 2)^2
            R = _ellipse_logdist(a, b, J * cm, (sx + s[j]) / 2) - log(abs(ψ[j] - z)) - log(abs(ψ[j] - conj(z)))
            A[i, j] = R * w[j]
        end
    end
    A .+= Mom * T′

    return -(1/(2π)) .* A .* transpose(γψ)

end

"""
    _laplace2d_double_layer_adjoint_phi(a, b, sa, sb, ψ, w, s_tgt)

Adjoint double-layer (`K′`) operator of a φ-parametrized elliptical arc, acting on its
physical density, at targets on the same ellipse (on or off the arc).

On the ellipse the kernel is `-ab / (4π |γ′(s_t)| (a² sin² m + b² cos² m))` with
`m = (s_t + s)/2`, which is smooth, includes the curvature limit at `s = s_t`, and avoids
cancellation for nearby points. It is integrated with the Fejér weights.

# Arguments
- `a`, `b`: semi-axes of the ellipse.
- `sa`, `sb`: parameter interval of the source arc.
- `ψ`, `w`: Fejér nodes and weights of the source arc.
- `s_tgt`: ellipse parameters of the targets.

# Returns
- `A::Matrix{Float64}` (`N_tgt × N`).
"""
function _laplace2d_double_layer_adjoint_phi(a, b, sa, sb, ψ, w, s_tgt)

    J = (sb - sa) / 2
    φ = π .* (1 .- ψ) ./ 2
    s = sa .+ J .* (cos.(φ) .+ 1)
    dΓ = _ellipse_speed.(a, b, s) .* J .* sin.(φ) .* (π / 2) .* w

    kernel(st, sj) = (m = (st + sj) / 2; -a * b / (4π * _ellipse_speed(a, b, st) * (a^2 * sin(m)^2 + b^2 * cos(m)^2)))
    return [kernel(st, s[j]) * dΓ[j] for st in s_tgt, j in eachindex(s)]

end
