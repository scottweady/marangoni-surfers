# Single-particle marangoni_solve against exact swimming speeds.
#
# 1. Circular surfer with g(s) = cos(2s): convergence in the arc resolution `Ns`, with the
#    Stokes order fixed, against (π/16)(1 + cos s₀) sin² s₀.
# 2. Elliptical surfer with g(s) = -cos(s): convergence in the Stokes expansion order `M`,
#    with the arcs fixed at `Ns_ref`, against the exact speed (as in
#    scripts/convergence/convergence.jl). Prolate and oblate aspect ratios exercise both
#    branches of the mobility difference Δm.

using SpecialFunctions: ellipk, ellipe

@testset "Circular surfer" begin
    section("Circular surfer: speed against exact, refining boundary resolution N")
    s₀ = π/2
    g(s) = cos(2s)
    Uexact = (π/16) * (1 + cos(s₀)) * sin(s₀)^2
    M = 32
    Nspan = sweep(4:2:12, [10])
    Ω_ref = ellipse(1.0, 1.0, M)

    errs = map(Nspan) do Ns
        Γ_D, Γ_N = surfer_boundary(1.0, 1.0, s₀, Ns)
        U, _ = marangoni_solve(Ω_ref, Γ_D, Γ_N, g)
        abs(U[1] - Uexact)
    end
    check_errors(Nspan, errs)
end

@testset "Elliptical surfer" begin
    section("Elliptical surfer: speed against exact, refining body resolution M")
    s₀ = π/2
    g(s) = -cos(s)
    Ns_ref = setting(4096, 256)
    Mspan = sweep(4:4:20, [8])

    # U = Δm/(2πab) ∫ x σ ds, with the disk dipole moment scaled by (a+b)/2 under the
    # Joukowski map
    function exact_speed(a, b)
        Δm = if a > b
            ε² = 1 - (b/a)^2
            b/(2ε²) * (ellipk(ε²) - ellipe(ε²)) - b/2 * ellipk(ε²)
        else
            ε² = 1 - (a/b)^2
            a/(2ε²) * (ellipe(ε²) - (1 - ε²) * ellipk(ε²)) - a/2 * ellipk(ε²)
        end
        return -(a + b) * (1 + cos(s₀))^2 * Δm / (8 * a * b)
    end

    for aspect_ratio in (4.0, 0.25)
        a = sqrt(aspect_ratio)
        b = 1.0 / a
        Uexact = exact_speed(a, b)
        Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns_ref)

        subsection("a/b = $aspect_ratio")
        errs = map(Mspan) do M
            U, _ = marangoni_solve(ellipse(a, b, M), Γ_D, Γ_N, g)
            abs(U[1] - Uexact)
        end
        check_errors(Mspan, errs; label="M")
    end
end
