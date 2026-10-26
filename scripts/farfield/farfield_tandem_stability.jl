# Linear stability of an infinite tandem chain of surfers in the far-field model with the
# monopole closure: particles spaced `d` apart on the x axis, all pointing along e_x, with
# transverse displacements and orientations (y, θ) ∝ exp(i κ j + λ t).
#
# About the uniform chain the closure amplitudes αᵢ vanish, and by symmetry under y → -y
# they stay zero to first order in (y, θ), so each wavenumber κ = kd decouples into
#
#     d/dt [Y; Θ] = [-i a  B; C  i β a] [Y; Θ],
#
#     a = η Cl₂(κ)/(2π d²),   B = U - η log(2 sin(κ/2))/(2π d),   C = β η S₃(κ)/(π d³),
#
# with Cl₂(κ) = Σ sin(nκ)/n² and S₃(κ) = Σ (1 - cos nκ)/n³, whose eigenvalues are
#
#     λ± = [i a (β - 1) ± √(4BC - a²(1 + β)²)] / 2.
#
# The discriminant is real, so each branch is continuous in κ: growing and decaying for
# β > 0, and two travelling waves for β < 0. The matrix at -κ is the conjugate of that at κ,
# so a real perturbation (cos κj, sin κj), which combines the two, has the eigenvalues λ±
# and their conjugates, Re λ± ± i |Im λ±|. Wavenumbers are sampled as the modes κ = 2π m/N
# of a chain of `N` particles.
#
# Output: output/tandem_stability_<a/b>.dat, one row per mode
#   m  κ  Re λ₊  |Im λ₊|  Re λ₋  |Im λ₋|

using QuadGK

include("../../src/script_params.jl")
include("../../src/farfield.jl")

aspect_ratios = script_param(:aspect_ratios, [0.25, 1.0, 4.0])
N             = script_param(:N, 100)
d             = script_param(:d, 10.0)
outdir        = script_param(:outdir, joinpath(@__DIR__, "output"))

"""
    clausen2(κ)

Clausen function `Cl₂(κ) = Σ sin(nκ)/n² = -∫₀^κ log(2 sin(t/2)) dt`, for `0 < κ < 2π`.
"""
clausen2(κ) = -quadgk(t -> log(2sin(t / 2)), 0, κ)[1]

"""
    tandem_eigenvalues(κ, d, U, η, β)

Eigenvalues `(λ₊, λ₋)` of the linearized dynamics of the transverse displacement and
orientation amplitudes `(Y, Θ)` of mode `κ = kd` of an infinite tandem chain with spacing `d`.
"""
function tandem_eigenvalues(κ, d, U, η, β)
    a = η * clausen2(κ) / (2π * d^2)
    B = U - η * log(2sin(κ / 2)) / (2π * d)
    C = β * η * quadgk(clausen2, 0, κ)[1] / (π * d^3)
    s = sqrt(complex(4B * C - a^2 * (1 + β)^2))
    return (im * a * (β - 1) + s) / 2, (im * a * (β - 1) - s) / 2
end

mkpath(outdir)
for aspect_ratio in aspect_ratios
    U, η, β = farfield_constants(aspect_ratio)
    file = joinpath(outdir, "tandem_stability_$(aspect_ratio).dat")
    open(file, "w") do io
        for m in 1:N÷2
            κ = 2π * m / N
            λ₊, λ₋ = tandem_eigenvalues(κ, d, U, η, β)
            @printf(io, "%d  %.10e  %.10e  %.10e  %.10e  %.10e\n", m, κ,
                    real(λ₊), abs(imag(λ₊)), real(λ₋), abs(imag(λ₋)))
        end
    end
    println("Wrote ", relpath(file, @__DIR__))
end
