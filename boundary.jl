
using FastGaussQuadrature

struct Boundary
    s  :: Vector{Float64}
    x  :: Vector{ComplexF64}
    xs :: Vector{ComplexF64}
    V  :: Matrix{Float64}
end

function legpts(n; dom=[-1,1])
    x, w = gausslegendre(n)
    a, b = dom
    x = 0.5 * (b - a) * x .+ 0.5 * (a + b)
    w = 0.5 * (b - a) * w
    return x, w
end

function trigpts(n; dom=[0, 2π])
    a, b = dom
    x = range(a, stop=b, length=n+1)[1:end-1]
    w = fill((b - a) / n, n)
    return x, w
end

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


function _laplace2d_single_layer_matrix(s, x, xs)

    N = length(s)
    wsing = kressweights(N)

    h = 2π / N
    A = zeros(Float64, N, N)

    for i in 1:N
        for j in 1:N
            k = mod(i - j, N) + 1
            speed = abs(xs[j])
            if i == j
                regkernel = -(1 / (2π)) * log(speed)
            else
                denom = 2 * abs(sin((s[i] - s[j]) / 2))
                regkernel = -(1 / (2π)) * log(abs(x[i] - x[j]) / denom)
            end
            A[i, j] = (wsing[k] + h * regkernel) * speed
        end
    end

    return A
end

function circle(N::Int)
    s, _ = trigpts(N)
    x  = exp.(1im .* s)
    xs = 1im .* exp.(1im .* s)
    return Boundary(s, x, xs, _laplace2d_single_layer_matrix(s, x, xs))
end

Base.length(Γ::Boundary) = length(Γ.s)
