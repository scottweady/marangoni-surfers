
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

"""
Smooth quadrature, Kress
"""
function _laplace2d_single_layer_kress(s, x, xs)

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

"""
Quadrature for 1/sqrt(d) densities, Bruno 2021
"""
function _laplace2d_single_layer_bruno(s, x, xs)

    N = length(s)

    # R^(N)(l) = -∑_{m=0}^{N-1} (2-δ_{0m}) λ_m cos(lmπ/N), Bruno 2021 eq. 4.4
    R = zeros(2N)
    for l in 0:2N-1
        for m in 0:N-1
            coeff = m == 0 ? 1.0 : 2.0
            λm    = m == 0 ? log(2) / 2 : 1 / (2m)
            R[l+1] -= coeff * λm * cos(l * m * π / N)
        end
    end

    A = zeros(Float64, N, N)

    for i in 1:N
        for j in 1:N
            speed_j = abs(xs[j])
            Rij = R[abs(i - j) + 1] + R[i + j]
            if i == j
                E2 = -(1 / 2π) * log(speed_j)
            else
                E2 = -(1 / 2π) * log(abs(x[i] - x[j]) / abs(s[i] - s[j]))
            end
            A[i, j] = (π / N) * speed_j * (-(1 / 2π) * Rij + E2)
        end
    end

    return A
end

function _laplace2d_double_layer_adjoint(s, x, xs, xss, w)

    N = length(s)
    A = zeros(Float64, N, N)
    n = (imag.(xs) - im * real.(xs)) ./ abs.(xs)

    for i in 1:N

        for j in 1:N
            if i == j
                A[i, j] =  (1/(4π)) * dot(n[i], xss[i]) / abs(xs[i])^3 * abs(xs[i]) * w[i]
            else
                δz = x[i] - x[j]
                A[i, j] = -(1/(2π)) * (dot(n[i], δz) / abs2(δz)) * abs(xs[j]) * w[j]
            end
        end

    end

    return A

end

# Cross single-layer: sources on Γ_src, targets on Γ_tgt (smooth, no singularity)
function _laplace2d_single_layer(Γ_src::Boundary, Γ_tgt::Boundary)

  N_src, N_tgt = length(Γ_src), length(Γ_tgt)
  A = zeros(Float64, N_tgt, N_src)

  for i in 1:N_tgt
    for j in 1:N_src
      A[i,j] = -(1/(2π)) * log(abs(Γ_tgt.x[i] - Γ_src.x[j])) * abs(Γ_src.xs[j]) * Γ_src.w[j]
    end
  end

  return A

end

# Cross adjoint double-layer: sources on Γ_src, normal at Γ_tgt (smooth, no singularity)
# Γ.V′ stores the raw K̃', so V′_standard = -(1/2π) Γ.V′; cross matrix uses the same convention
function _laplace2d_double_layer_adjoint(Γ_src::Boundary, Γ_tgt::Boundary)

  N_src, N_tgt = length(Γ_src), length(Γ_tgt)
  A = zeros(Float64, N_tgt, N_src)

  for i in 1 : N_tgt
    n = (imag(Γ_tgt.xs[i]) - im * real(Γ_tgt.xs[i])) / abs.(Γ_tgt.xs[i])
    for j in 1 : N_src
      δz = Γ_tgt.x[i] - Γ_src.x[j]
      A[i,j] = -(1/(2π)) * (dot(n, δz) / abs2(δz)) * abs(Γ_src.xs[j]) * Γ_src.w[j]
    end
  end

  return A

end