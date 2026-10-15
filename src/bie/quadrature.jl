
"""
    dot(a, b)

Real inner product of two complex numbers, treating each as a 2-vector.

# Arguments
- `a`, `b`: complex numbers, interpreted as `(real, imag)` pairs.

# Returns
- `real(a)*real(b) + imag(a)*imag(b)`.
"""
dot(a::Complex, b::Complex)   = real(a)*real(b) + imag(a)*imag(b)

"""
    cross(a, b)

Real (scalar/z-component) cross product of two complex numbers, treating each as a 2-vector.

# Arguments
- `a`, `b`: complex numbers, interpreted as `(real, imag)` pairs.

# Returns
- `real(a)*imag(b) - imag(a)*real(b)`.
"""
cross(a::Complex, b::Complex) = real(a)*imag(b) - imag(a)*real(b)

"""
    integral(σ, Γ)

Boundary integral `∫_Γ σ ds = ∑ⱼ σⱼ |γ′ⱼ| wⱼ`.

# Arguments
- `σ`: density sampled at the quadrature nodes of `Γ`.
- `Γ::Boundary`: boundary to integrate over.

# Returns
- The scalar value of the integral.
"""
function integral(σ, Γ::Boundary)
    return sum(σ .* Γ.dΓ)
end

"""
    _fejerpts(n; dom=[-1,1])

Fejér first-rule nodes (Chebyshev nodes of the first kind) and weights, mapped onto `dom`.

# Arguments
- `n`: number of nodes.
- `dom`: target interval (default `[-1,1]`).

# Returns
- `(t, w)`: node and weight vectors.
"""
function _fejerpts(n; dom=[-1,1])
    θ = [π * (2j - 1) / (2n) for j in 1:n]
    t = cos.(θ[end:-1:1])
    w = [(2/n) * (1 - 2 * sum(cos(2k * θ[j]) / (4k^2 - 1) for k in 1:n÷2)) for j in n:-1:1]
    a, b = dom
    t = 0.5 * (b - a) * t .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return t, w
end

"""
    _legpts(n; dom=[-1,1])

Gauss-Legendre nodes and weights, scaled from `[-1,1]` onto `dom`.

# Arguments
- `n`: number of nodes.
- `dom`: target interval (default `[-1,1]`).

# Returns
- `(x, w)`: node and weight vectors.
"""
function _legpts(n; dom=[-1,1])
    x, w = gausslegendre(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

"""
    _chebpts(n; dom=[-1,1])

Gauss-Chebyshev nodes and weights (Chebyshev-T, weight `1/√(1-x²)`), scaled from `[-1,1]` onto `dom`.

# Arguments
- `n`: number of nodes.
- `dom`: target interval (default `[-1,1]`).

# Returns
- `(x, w)`: node and weight vectors.
"""
function _chebpts(n; dom=[-1,1])
    x, w = gausschebyshevt(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

"""
    _trigpts(n; dom=[0, 2π])

Uniform trigonometric nodes and weights on `dom` (periodic trapezoid rule).

# Arguments
- `n`: number of nodes.
- `dom`: target interval (default `[0, 2π]`).

# Returns
- `(x, w)`: node and weight vectors.
"""
function _trigpts(n; dom=[0, 2π])
    a, b = dom
    x = range(a, stop=b, length=n+1)[1:end-1]
    w = fill((b - a) / n, n)
    return x, w
end

"""
    _cheb_coeffs(σ; Nmax=length(σ))

Chebyshev coefficients of a function sampled at Chebyshev nodes, via the DCT.

# Arguments
- `σ`: (possibly complex) samples at the `_chebpts` nodes.
- `Nmax`: number of coefficients to keep (default: all of them).

# Returns
- The first `Nmax` Chebyshev coefficients.
"""
function _cheb_coeffs(σ; Nmax=length(σ))

    Nmax = min(Int(Nmax), length(σ))

    N = length(σ)
    c_real = FFTW.r2r(reverse(real.(σ)), FFTW.REDFT10) ./ N
    c_imag = FFTW.r2r(reverse(imag.(σ)), FFTW.REDFT10) ./ N
    c_real[1] /= 2
    c_imag[1] /= 2
    c = c_real .+ im .* c_imag
    return c[1:Nmax]

end

"""
    _cheb_eval(c, t)

Evaluate the Chebyshev series `Σ cₖ Tₖ(t)` at (possibly complex) `t` via Clenshaw recurrence.

# Arguments
- `c`: Chebyshev coefficients, as returned by [`_cheb_coeffs`](@ref).
- `t`: evaluation point(s), real or complex.

# Returns
- The value of the series at `t`.
"""
function _cheb_eval(c, t)
    n = length(c)
    b1, b2 = zero(t), zero(t)
    for k = n-1:-1:1
        b0 = c[k+1] + 2t * b1 - b2
        b2, b1 = b1, b0
    end
    return c[1] + t * b1 - b2
end

"""
    _kresspts(N)

Kress singular quadrature weights for the log-kernel on a closed curve.

# Arguments
- `N`: number of nodes (must be even).

# Returns
- `(θ, w)`: the N uniform angular nodes and the N-vector `R^(N)` such that the `(i,j)`
  entry of the single-layer matrix is `(R[mod(i-j,N)+1] + h * regularKernel[i,j]) * speed[j]`.
"""
function _kresspts(N::Int)

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

    return θ, w

end

"""
    _cheb_matrix(M)

Forward and inverse Chebyshev transform matrices on the `M`-point Gauss-Chebyshev grid.

# Arguments
- `M`: number of Chebyshev nodes.

# Returns
- `(T, T′)`: `T` evaluates a Chebyshev series at the nodes (`T[i,j] = T_{j-1}(t_i)`);
  `T′` maps node samples to Chebyshev coefficients (its quadrature-weighted pseudoinverse).
"""
function _cheb_matrix(M)
    t, w  = _chebpts(M)
    T     = [cos((j-1) * acos(t[i])) for i in 1:M, j in 1:M]
    norms = [k == 1 ? π : π/2 for k in 1:M]
    T′    = (T .* w)' ./ norms
    return T, T′
end

"""
    _single_layer_eig_weighted(M)

Eigenvalues of the Laplace single-layer operator restricted to the first `M` Chebyshev modes.

Diagonal in the Chebyshev basis: mode `1` (constant) maps to `log(2)/2`, mode `m` (`m>1`)
maps to `1/(2(m-1))`.

# Arguments
- `M`: number of Chebyshev modes.

# Returns
- `Diagonal` matrix of eigenvalues.
"""
function _single_layer_eig_weighted(M)
    return Diagonal([m == 1 ? log(2)/2 : 1/(2*(m-1)) for m in 1:M])
end

"""
    _log_moments(x, M)

Moments `mₖ(x) = ∫₋₁¹ log|x - t| Tₖ(t) dt` of the log kernel against the first `M`
Chebyshev polynomials, by a closed form.

Integrating `Re log(x - t)` by parts against an antiderivative `Aₖ` of `Tₖ` with `Aₖ(-1) = 0` gives

    mₖ(x) = Aₖ(1) log|1 - x| - Re[Aₖ(x) L(x)] - Re ∫₋₁¹ (Aₖ(t) - Aₖ(x))/(t - x) dt,

with `L(x) = log(x - 1) - log(x + 1)`, and the last integral follows from
`Dₙ(x) = ∫₋₁¹ (Tₙ(t) - Tₙ(x))/(t - x) dt`, which satisfies
`D_{n+1} = 2x Dₙ + 2∫₋₁¹ Tₙ dt - D_{n-1}` with `D₀ = 0`, `D₁ = 2`.

The formula holds for any `x`, but off `(-1, 1)` it cancels terms of size `|Tₖ(x)|`, losing a
factor of about `|r|⁻ᴹ` with `r = x - √(x² - 1)`, `|r| ≤ 1`. It is therefore used only on the
interval and very near it; see [`_log_moments_ext`](@ref).

# Arguments
- `x`: real or complex target.
- `M`: number of moments.

# Returns
- Vector `[m₀(x), …, m_{M-1}(x)]`.
"""
function _log_moments(x, M)

    xc = complex(float(x))
    L  = log(xc - 1) - log(xc + 1)
    l1 = log(abs(1 - x))

    # Tₙ(x) and Dₙ(x) for n = 0..M
    T = zeros(typeof(xc), M + 1)
    D = zeros(typeof(xc), M + 1)
    T[1] = 1.0
    T[2] = xc
    D[2] = 2.0
    intT(n) = n == 1 ? 0.0 : (1 + (-1)^n) / (1 - n^2)
    for n in 1:M-1
        T[n+2] = 2xc * T[n+1] - T[n]
        D[n+2] = 2xc * D[n+1] + 2 * intT(n) - D[n]
    end

    m = zeros(M)
    m[1] = 2l1 - real((1 + xc) * L) - 2
    M >= 2 && (m[2] = -real((xc^2 - 1) / 2 * L + xc))
    for k in 2:M-1
        A  = (T[k+2] / (k + 1) - T[k] / (k - 1)) / 2 - (-1)^k / (k^2 - 1)
        A1 = -(1 + (-1)^k) / (k^2 - 1)
        m[k+1] = A1 * l1 - real(A * L + (D[k+2] / (k + 1) - D[k] / (k - 1)) / 2)
    end

    return m

end

"""
    _log_moments_ext(z, M)

Moments `mₖ(z) = ∫₋₁¹ log|z - t| Tₖ(t) dt` for `z` off the open interval `(-1, 1)`, real or
complex.

The closed form of [`_log_moments`](@ref) cancels catastrophically here, since `Tₖ(z)` grows
geometrically. Instead, with `t = cos θ` and `r = z - √(z² - 1)` (`|r| < 1`),

    log|z - cos θ| = -Re log(2r) - 2 Σₙ Re(rⁿ) cos(nθ)/n,

and each term integrates against `cos(kθ) sin θ` in closed form via
`F(j) = ∫₀^π cos(jθ) sin θ dθ = (1 + (-1)ʲ)/(1 - j²)`. The series converges like `|r|ⁿ`,
so it is slow when `z` is close to an endpoint. That is exactly where the closed form of
[`_log_moments`](@ref) is accurate, since it loses only a factor of about `|r|⁻ᴹ`, so it is
used whenever that loss is small.

# Arguments
- `z`: real point with `|z| ≥ 1`, or complex point.
- `M`: number of moments.

# Returns
- Vector `[m₀(z), …, m_{M-1}(z)]`.
"""
function _log_moments_ext(z, M)

    zc = complex(float(z))
    r₀ = zc - sqrt(zc - 1) * sqrt(zc + 1)
    r = abs(r₀) > 1 ? 1 / r₀ : r₀
    -M * log(abs(r)) < 5 && return _log_moments(z, M)
    nmax = abs(r) < 1e-300 ? 1 : ceil(Int, log(1e-17) / log(abs(r)))

    F(j) = isodd(j) ? 0.0 : 2 / (1 - j^2)
    c = [real(r^n) / n for n in 1:nmax]
    l0 = -real(log(2r))
    m = zeros(M)
    for k in 0:M-1
        acc = l0 * F(k)
        for n in 1:nmax
            isodd(n + k) && continue
            acc -= c[n] * (F(n + k) + F(n - k))
        end
        m[k+1] = acc
    end

    return m

end
