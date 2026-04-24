
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
