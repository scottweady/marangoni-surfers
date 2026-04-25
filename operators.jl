
using FFTW
using LinearAlgebra
using QuadGK

function dot(a::Complex, b::Complex)
    return real(a)*real(b) + imag(a)*imag(b)
end

function cross(a::Complex, b::Complex)
    return real(a)*imag(b) - imag(a)*real(b)
end

"""
Boundary-to-boundary single layer potential

# Arguments
- `Γ_src::Boundary`: source boundary
- `Γ_tgt::Boundary`: target boundary

# Returns
- Single layer potential evaluated on Γ
"""
function V(σ, Γ::Boundary)
  return Γ.V * σ
end

"""
Closed curves
"""
function V⁻¹(g, Γ::Boundary)

  N = length(Γ)
  e = ones(N)
  ds = abs.(Γ.xs) .* Γ.w

  A = [Γ.V   e;
       ds' 0.0]
  rhs  = [g; 0.0]

  sol = A \ rhs

  return sol[1:N]

end

function V′(σ, Γ::Boundary)
  return Γ.V′ * σ
end

function _smooth_evaluation(ζ, σ, kernel, Γ::Boundary)
  shp = size(ζ)
  ζ = vec(ζ)
  δz = ζ .- transpose(Γ.x)
  A = kernel.(δz) .* transpose(σ .* abs.(Γ.xs) .* Γ.w)
  return reshape(sum(A, dims=2), shp)
end

# Evaluate the single-layer potential at an off-curve target point r
function 𝒱(ζ, σ, Γ::Boundary)
  return _smooth_evaluation(ζ, σ, r -> -(1/(2π)) * log(abs(r)), Γ)
end

function ℬ(ζ, σ, Γ::Boundary)
  return _smooth_evaluation(ζ, σ, r -> abs(r)/(4π), Γ)
end

function ∇ℬ(ζ, σ, Γ::Boundary; rtol=1e-8)

    N = length(Γ)
    freqs = vcat(0:N÷2-1, -(N÷2):-1)

    σ̂  = fft(σ) ./ N
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

# Evaluate ∫_Γ (ζ-x')/|ζ-x'| α(x') ds'/(4π) for arc Boundary (Bruno amplitude α).
# Change variable s=cosφ removes the 1/√(1-s²) weight, leaving a smooth integrand for QuadGK.
function ∇ℬ_arc(ζ::Number, α, Γ::Boundary; rtol=1e-8)

    N = length(Γ)

    function dct_coeffs(v)
        raw = FFTW.r2r(collect(Float64, v), FFTW.REDFT10)
        ĉ = raw ./ N
        ĉ[1] /= 2
        return ĉ
    end

    function coseval(ĉ, φ)
        s = ĉ[1]
        for k in 1:N-1
            s += ĉ[k+1] * cos(k * φ)
        end
        return s
    end

    α̂  = dct_coeffs(α)
    x̂r = dct_coeffs(real.(Γ.x))
    x̂i = dct_coeffs(imag.(Γ.x))
    speed = abs(Γ.xs[1])  # |xs| is constant for unit-circle arc

    function integrand(φ)
        xi = coseval(x̂r, φ) + im * coseval(x̂i, φ)
        αi = coseval(α̂, φ)
        δz = ζ - xi
        return (δz / abs(δz)) * αi * speed / (4π)
    end

    val, _ = quadgk(integrand, 0.0, π; rtol=rtol)
    return real(val), imag(val)
end

function ∇ℬ_arc(ζ, α, Γ::Boundary; rtol=1e-8)
    shp = size(ζ)
    results = [∇ℬ_arc(ζi, α, Γ; rtol=rtol) for ζi in vec(ζ)]
    ux = reshape(first.(results), shp)
    uy = reshape(last.(results),  shp)
    return ux, uy
end

function ∫(σ, Γ::Boundary)
  return sum(σ .* abs.(Γ.xs)) * (2π / length(Γ))
end

function solve_mixed(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)

    N_D, N_N = length(Γ_D), length(Γ_N)

    V_DD  = Γ_D.V
    V′_NN = Γ_N.V′
    V_DN  = _laplace2d_single_layer(Γ_N, Γ_D)
    V′_ND = _laplace2d_double_layer_adjoint(Γ_D, Γ_N)

    I_N = Diagonal(1 ./ sqrt.(1 .- Γ_N.s .^ 2))

    # +1/2 jump for adjoint
    A = [V_DD            V_DN;
         V′_ND  I_N/2 + V′_NN]

    # Mean-zero condition on σ
    ds = vcat((π/N_D) .* abs.(Γ_D.xs),
              (π/N_N) .* abs.(Γ_N.xs))
    e = vcat(ones(N_D), zeros(N_N))
    A = [A      e; 
         ds'  0.0]
    σ = A \ [g_D; g_N; 0.0]

    display(V′_ND)

    return σ[1:N_D], σ[N_D+1:N_D+N_N]
end
