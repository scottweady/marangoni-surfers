
# Real inner product and cross product treating complex numbers as 2-vectors.
dot(a::Complex, b::Complex)   = real(a)*real(b) + imag(a)*imag(b)
cross(a::Complex, b::Complex) = real(a)*imag(b) - imag(a)*real(b)

"Boundary integral: ∫_Γ σ ds = ∑ⱼ σⱼ |γ′ⱼ| wⱼ."
function integral(σ, Γ::Boundary)
    return sum(σ .* Γ.dΓ)
end

function integral(μ, Ω::Domain)
    return sum(μ .* Ω.dz)
end

"Fejér first-rule nodes (Chebyshev nodes of first kind) and weights on [-1,1]."
function _fejerpts(n; dom=[-1,1])
    θ = [π * (2j - 1) / (2n) for j in 1:n]
    t = cos.(θ[end:-1:1])
    w = [(2/n) * (1 - 2 * sum(cos(2k * θ[j]) / (4k^2 - 1) for k in 1:n÷2)) for j in n:-1:1]
    a, b = dom
    t = 0.5 * (b - a) * t .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return t, w
end

"Gauss-Legendre nodes and weights on `dom`, scaled from [-1,1]."
function _legpts(n; dom=[-1,1])
    x, w = gausslegendre(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

function _chebpts(n; dom=[-1,1])
    x, w = gausschebyshevt(n)
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

"Chebyshev coefficients of a function sampled at Chebyshev nodes, via DCT."
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

"Evaluate Chebyshev series Σ cₖ Tₖ(t) at (possibly complex) t via Clenshaw recurrence, skipping trailing small coefficients."
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
Kress singular quadrature weights for the log-kernel on a closed curve.

Returns the N-vector R^(N) such that the (i,j) entry of the single-layer matrix is
`(R[mod(i-j,N)+1] + h * regularKernel[i,j]) * speed[j]`.
Reference: Kress, Linear Integral Equations, §12.
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
