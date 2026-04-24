
using FFTW
using LinearAlgebra
using QuadGK

function V(σ, Γ::Boundary)
    return Γ.V * σ
end

function V⁻¹(g, Γ::Boundary)

    N = length(Γ)
    w = (2π/N) * abs.(Γ.xs)

    Kaug = [Γ.V  w;
            w' 0.0]
    rhs  = [g; 0.0]

    sol = Kaug \ rhs

    return sol[1:N]

end

function V′(σ, Γ::Boundary)
    return Γ.V′ * σ
end

# Evaluate the single-layer potential at an off-curve target point r
function 𝒱(ζ, sigma, Γ::Boundary)
    N = length(Γ)
    h = 2π / N
    u = 0.0
    for j in 1:N
        u += -(1 / (2π)) * log(abs(ζ - Γ.x[j])) * sigma[j] * abs(Γ.xs[j]) * h
    end
    return u
end

function ℬ(ζ, sigma, Γ::Boundary)

    shp = size(ζ)
    ζ = vec(ζ)
    h = 2π / length(Γ)

    δz = ζ .- transpose(Γ.x)
    r = abs.(δz)
    A = r .* transpose(sigma .* abs.(Γ.xs) * h) / 4π

    return reshape(sum(A, dims=2), shp)

end

function ∇ℬ(ζ, sigma, Γ::Boundary; rtol=1e-8)

    N = length(Γ)
    freqs = vcat(0:N÷2-1, -(N÷2):-1)

    σ̂  = fft(sigma) ./ N
    x̂  = fft(Γ.x)   ./ N
    x̂s = fft(Γ.xs)  ./ N

    function trig_eval(ĉ, t)
        sum(ĉ[k+1] * cis(freqs[k+1] * t) for k in 0:N-1)
    end

    function integrand(t, ζi)
        xi  = trig_eval(x̂,  t)
        xsi = trig_eval(x̂s, t)
        σi  = real(trig_eval(σ̂, t))
        δz  = ζi - xi
        return (δz / abs(δz)) * σi * abs(xsi) / (4π)
    end

    shp = size(ζ)
    vals = [quadgk(t -> integrand(t, ζi), 0.0, 2π; rtol=rtol)[1] for ζi in vec(ζ)]
    return (real.(reshape(vals, shp)), imag.(reshape(vals, shp)))

end

function ∂∫(σ, Γ::Boundary)
  return sum(σ .* abs.(Γ.xs)) * (2π / length(Γ))
end

function solve_mixed(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)
    N_D, N_N = length(Γ_D), length(Γ_N)

    V_DD  = Γ_D.V
    V′_NN = -(1/2π) .* Γ_N.V′
    V_DN  = _single_layer_cross_matrix(Γ_N, Γ_D)
    V′_ND = _adjoint_double_layer_cross_matrix(Γ_D, Γ_N)

    jump_N = Diagonal(1 ./ (2 .* sqrt.(1 .- Γ_N.s .^ 2)))

    A = [V_DD                    V_DN;
         V′_ND  Matrix(jump_N) + V′_NN]

    α = A \ [g_D; g_N]
    return α[1:N_D], α[N_D+1:end]
end
