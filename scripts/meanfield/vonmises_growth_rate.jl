# Maximal growth rate of the linearized mean-field model about the von Mises state
# ψ₀ ∝ exp(ζ cos θ), as a function of wavenumber k and wavevector angle φ, for each
# alignment strength in `ζs`. Growth rates are roots of the dispersion relation, refined
# by Newton from the eigenvalues of the discretized operator.
#
# Output: output/vonmises_growth_rate.dat
#   Line 1: # Nk=... Nphi=... Nzeta=... a=... beta=... U=... eta=... N=...
#   Line 2: k values (Nk floats)
#   Line 3: φ values (Nphi floats)
#   Line 4: ζ values (Nzeta floats)
#   Then Nzeta blocks of Nphi rows × Nk columns: maximal growth rate at (φ, k)

using LinearAlgebra
using Printf
using QuadGK
using SpecialFunctions

# -- Parameters ----------------------------------------------------------------

aspect_ratio = 4.0
a = sqrt(aspect_ratio)
b = 1.0 / a

β = (a^2 - b^2) / (a^2 + b^2)
U = 1.0
η = 1.0

N  = 32   # angular resolution of the discretized operator
ζs = [5.0, 10.0, 15.0]
kh = range(0, 6, length=200)
ph = range(0, π/2, length=200)

# -- Dispersion relation -------------------------------------------------------

"Dispersion function `D(σ)` and its derivative `dD/dσ` at wavevector `k (cos φ, sin φ)`."
function dispersion(σ, k, φ, ζ, η, U, β)
    c(t)   = cos(t - φ)     # k̂ ⋅ p
    s(t)   = -sin(t - φ)    # k̂ ⋅ p⊥
    ψ₀(t)  = exp(ζ * (cos(t) - 1)) / (2π * besselix(0, ζ))
    h(t)   = -(1 + β) + 2β * c(t)^2 + ζ * β * sin(t) * s(t) * c(t)
    den(t) = σ + im * U * k * c(t)
    I1, = quadgk(t -> c(t) * h(t) * ψ₀(t) / den(t), 0, 2π; rtol=1e-6)
    I2, = quadgk(t -> c(t) * h(t) * ψ₀(t) / den(t)^2, 0, 2π; rtol=1e-6)
    return 1 - (im * η / 2) * I1, (im * η / 2) * I2
end

"Newton iteration for a root of `f`; gives up once the iterate leaves Re σ > 0, since only those roots are kept."
function newton(f, s)
    for _ in 1:25
        D, dD = f(s)
        ds = D / dD
        s -= ds
        real(s) <= 0 && return s, false
        abs(ds) < 1e-12 * max(1, abs(s)) && return s, true
    end
    return s, false
end

"Eigenvalues of the linearized operator discretized at `N` orientations."
function operator_eigs(k, φ, ζ, η, U, β, N)
    t  = range(-π, π, length=N+1)[1:end-1]
    w  = step(t)
    p  = [cos.(t) sin.(t)]
    p⊥ = [-sin.(t) cos.(t)]
    k̂  = [cos(φ), sin(φ)]

    ψ₀  = exp.(ζ .* (cos.(t) .- 1)) ./ (2π * besselix(0, ζ))
    dψ₀ = -ζ .* sin.(t) .* ψ₀

    m₁  = w * p'
    ku₁ = (η / 2) * k̂ * (k̂' * m₁)
    E₁(A, B) = (im / 2) * ((A * k̂) .* (B * ku₁) .+ (A * ku₁) .* (B * k̂))

    θ̇₁  = β * E₁(p⊥, p)
    dθ̇₁ = β * (E₁(p⊥, p⊥) - E₁(p, p))

    L = -im * k * Diagonal(U * p * k̂) - ψ₀ * (im * k̂' * ku₁) - (dθ̇₁ .* ψ₀ + θ̇₁ .* dψ₀)
    return eigvals(L)
end

"Root with the largest real part (NaN if none with Re σ > 0), refined from the operator eigenvalues."
function max_growth(k, φ, ζ, η, U, β, N)
    σ₀ = filter(σ -> real(σ) > 1e-12, operator_eigs(k, φ, ζ, η, U, β, N))
    f = σ -> dispersion(σ, k, φ, ζ, η, U, β)
    smax = complex(NaN)
    for s0 in σ₀
        s, ok = newton(f, s0)
        if ok && (isnan(smax) || real(s) > real(smax))
            smax = s
        end
    end
    return smax
end

"""
Maximal growth rate along a row of fixed φ. Marches up in k, continuing the root from the
previous k; falls back to the full search when that fails, and stops once `nstop`
consecutive points are stable after an unstable band.
"""
function growth_row(kh, φ, ζ, η, U, β, N)
    tol   = 1e-6   # growth rate below which a point doesn't mark the row unstable
    nstop = 10
    srow  = zeros(length(kh))
    s = complex(NaN)
    unstable = false
    nstable  = 0
    for (n, k) in enumerate(kh)
        ok = false
        if !isnan(s)
            s, ok = newton(σ -> dispersion(σ, k, φ, ζ, η, U, β), s)
        end
        ok || (s = max_growth(k, φ, ζ, η, U, β, N))
        if isnan(s)
            nstable += 1
            unstable && nstable >= nstop && break
            continue
        end
        nstable  = 0
        unstable = unstable || real(s) > tol
        srow[n]  = real(s)
    end
    return srow
end

# -- Sweep ---------------------------------------------------------------------

Smax = [zeros(length(ph), length(kh)) for _ in ζs]
t0 = time()
for (j, ζ) in enumerate(ζs), (i, φ) in enumerate(ph)
    Smax[j][i, :] = growth_row(kh, φ, ζ, η, U, β, N)
    @printf("ζ = %g   φ %3d/%d   max growth = %.3g   [%.0f s]\n",
            ζ, i, length(ph), maximum(Smax[j][i, :]), time() - t0)
end

# -- Write results -------------------------------------------------------------

write_row(io, v) = println(io, join([@sprintf("%.10e", x) for x in v], " "))

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
outfile = joinpath(outdir, "vonmises_growth_rate.dat")
open(outfile, "w") do io
    @printf(io, "# Nk=%d Nphi=%d Nzeta=%d a=%.10f beta=%.10f U=%.10f eta=%.10f N=%d\n",
            length(kh), length(ph), length(ζs), a, β, U, η, N)
    write_row(io, kh)
    write_row(io, ph)
    write_row(io, ζs)
    foreach(S -> foreach(i -> write_row(io, S[i, :]), axes(S, 1)), Smax)
end

println("Wrote ", relpath(outfile, @__DIR__))
