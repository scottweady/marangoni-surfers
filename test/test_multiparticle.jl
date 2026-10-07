# Multi-particle marangoni_solve on an elliptical surfer with data g that breaks the up-down
# symmetry, so all three components of U are nonzero.
#
# 1. One particle at the origin: U should match the single-particle marangoni_solve, since
#    the self flow is projected with the same reciprocal construction.
# 2. One particle placed at (xcm, θcm): U should be the origin result with its translational
#    part rotated by θcm, to within the FMM/GMRES tolerance.
# 3. Two particles: convergence in `M` against a reference at `M_ref`.
# 4. Two circular surfers: φ-arcs (the default) against the mixed arcs, which are
#    independent of the φ-arc operators and accurate for circles at high resolution. This
#    checks the particle-to-particle blocks, which use the arcs' tangents and normals.

@testset "Multi-particle solve" begin
    section("Multi-particle solve, refining body resolution M")
    s₀ = π/3
    g(s) = -cos(s) + 0.4 * sin(2s)
    a = sqrt(2.0)
    b = 1.0 / a
    Ns_ref = 256
    tol = 1e-10
    Mspan_single = sweep(4:4:20, [8])
    Mspan_pair = sweep(8:4:24, [12])
    M_ref = setting(32, 16)

    xcm, θcm = 1.3 - 0.7im, 0.9
    R = [cos(θcm) -sin(θcm); sin(θcm) cos(θcm)]

    xcm_pair, θcm_pair = [0.0im, 1.0 + 3.0im], [0.0, 2.0]

    Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns_ref)

    function pair_velocity(M)
        Ω = ellipse(a, b, M)
        self_op = marangoni_self_operator(Ω, Γ_D, Γ_N)
        U, = marangoni_solve(xcm_pair, θcm_pair, Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)
        return vcat(U...)
    end

    err_origin, err_placed = Float64[], Float64[]
    for M in Mspan_single
        Ω = ellipse(a, b, M)
        self_op = marangoni_self_operator(Ω, Γ_D, Γ_N)

        U_single, = marangoni_solve(Ω, Γ_D, Γ_N, g)
        U_origin, = marangoni_solve([0.0im], [0.0], Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)
        U_placed, = marangoni_solve([xcm], [θcm], Ω, Γ_D, Γ_N, g; tol=tol, self_op=self_op)
        U_rot = [R * U_origin[1][1:2]; U_origin[1][3]]

        push!(err_origin, norm(U_origin[1] .- U_single, Inf))
        push!(err_placed, norm(U_placed[1] .- U_rot, Inf))
    end

    U_pair_ref = pair_velocity(M_ref)
    err_pair = [norm(pair_velocity(M) .- U_pair_ref, Inf) for M in Mspan_pair]

    subsection("one particle at the origin, against single-particle solve")
    check_errors(Mspan_single, err_origin; label="M", every=true)

    subsection("one particle translated and rotated, against rotated origin result")
    check_errors(Mspan_single, err_placed; label="M", every=true)

    subsection("two particles, against high-M reference")
    check_errors(Mspan_pair, err_pair; label="M")
end

@testset "Multi-particle solve: φ-arcs against mixed arcs" begin
    section("Multi-particle solve: φ-arcs against mixed arcs, refining boundary resolution N")
    s₀ = π/2
    g(s) = -cos(s)
    M = 16
    tol = 1e-10
    Nspan = sweep(24:12:72, [48])
    Ns_mixed = 256

    xcm, θcm = [0.0im, 2.5 + 1.0im], [0.3, 2.0]
    Ω = ellipse(1.0, 1.0, M)
    pair_velocity(Γ_D, Γ_N) = vcat(first(marangoni_solve(xcm, θcm, Ω, Γ_D, Γ_N, g; tol=tol))...)

    U_mixed = pair_velocity(surfer_boundary(1.0, 1.0, s₀, Ns_mixed; quad=:mixed)...)
    errs = [norm(pair_velocity(surfer_boundary(1.0, 1.0, s₀, Ns)...) .- U_mixed, Inf) for Ns in Nspan]
    check_errors(Nspan, errs)
end
