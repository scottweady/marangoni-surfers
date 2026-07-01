
using ProjectedSphericalHarmonics

integral(f, Ω::Domain) = sum(f .* Ω.dz)

"Evaluate a domain field f at arc parameter values θ by tracing to the boundary and FFT-interpolating."
function _boundary_eval(f, Ω, s)
    vals  = trace(f, Ω)
    N     = length(vals)
    f̂     = fft(vals) ./ N
    freqs = fftfreq(N) * N
    return real.(sum(f̂[k] .* cis.(freqs[k] .* s) for k in 1:N))
end

"""
Solve the surfactant Laplace problem given reference arc boundaries.

`Γ_D` and `Γ_N` are reference (body-frame) Dirichlet and Neumann arc boundaries.
They are rigidly placed at each `(centers[k], orientations[k])` — the self-interaction
matrices `V`, `K′` are reused since the Laplace kernel is translation/rotation invariant.
`g` is the Dirichlet BC as a function of arc parameter s.
Returns `(densities, boundaries, ϕ∞)`.
"""
function surfactant_solve(xcm, θcm, Γ_D::Boundary, Γ_N::Boundary, g)

    Γ = Boundary[]
    bc = Symbol[]
    rhs = Vector{Float64}[]

    for (xcm_i, θcm_i) in zip(xcm, θcm)
        Γ_D_i = _place(Γ_D, xcm_i, θcm_i)
        Γ_N_i = _place(Γ_N, xcm_i, θcm_i)
        push!(Γ, Γ_D_i, Γ_N_i)
        push!(bc, :D, :N)
        push!(rhs, g.(Γ_D_i.s), zeros(length(Γ_N_i.s)))
    end

    σ, ϕinf = laplace_solve(Γ, bc, rhs)
    return σ, ϕinf, Γ

end

"""
Solve for the translational velocity U and angular velocity ω of a single surfer.

`Ω` is the Stokes domain boundary, `Γ_D` and `Γ_N` are the Dirichlet and Neumann
arc boundaries for the surfactant solve, and `g` is the Dirichlet BC as a function
of arc parameter s.
Returns `(U, σ, Γ, ϕ∞)`.
"""
function marangoni_solve(xcm, θcm, Ω, Γ_D::Boundary, Γ_N::Boundary, g; tol=1e-6)

    Nparticles = length(xcm)

    σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

    pts  = [exp(im * θcm[i]) .* Ω.z .+ xcm[i] for i in 1:Nparticles]
    Npts     = length(Ω.z)
    offs     = [0; cumsum(fill(Npts, Nparticles))]
    bx, by = ∇ℬ(vcat(vec.(pts)...), σ, Γ; eps=tol)

    uinf = Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, Nparticles)
    for i in 1:Nparticles
        ri = offs[i]+1:offs[i+1]
        shp = size(pts[i])
        uinf[i] = (reshape(copy(bx[ri]), shp), reshape(copy(by[ri]), shp))
    end

    F = [zeros(3) for _ in 1:Nparticles]
    _, U = suspension_mobility_solve(xcm, θcm, uinf, nothing, F, Ω; eps=tol, gmrestol=tol)

    return U, σ, Γ, ϕinf

end

"Single-particle Marangoni solve via mobility formulation."
function marangoni_solve(Ω, Γ_D, Γ_N, g)

    σ, ϕinf, Γ = surfactant_solve([0.0im], [0.0], Γ_D, Γ_N, g)
    σ_D, σ_N = σ[1], σ[2]
    x, y = real.(Ω.z), imag.(Ω.z)

    function _marangoni_force(e)
        μ = 𝒮_st⁻¹(e, Ω)
        h = 𝒮(μ[1] .* x + μ[2] .* y, Ω) - (x .* 𝒮(μ[1], Ω) + y .* 𝒮(μ[2], Ω))
        h_D = _boundary_eval(h, Ω, Γ_D.s)
        h_N = _boundary_eval(h, Ω, Γ_N.s)
        return integral(h_D .* σ_D, Γ_D) + integral(h_N .* σ_N, Γ_N)
    end

    F = [_marangoni_force((1.0, 0.0)), _marangoni_force((0.0, 1.0)), _marangoni_force((-y, x))]
    M = stokes_mobility_matrix(Ω)
    U = real.(M \ F)

    return U, σ, Γ, ϕinf

end

"""
Compute the Marangoni multipole moments for a single surfer.

Returns a NamedTuple with fields `velocity`, `force_monopole`, `force_dipole`, `surf_monopole`, `surf_dipole`.
"""
function marangoni_multipoles(Ω, Γ_D, Γ_N, g)

    U, σ, _, _ = marangoni_solve(Ω, Γ_D, Γ_N, g)

    ucm, ωcm = U[1:2], U[3]
    α_D, α_N = σ

    z    = Ω.z
    x, y = real.(z), imag.(z)

    bx = ∇ℬ(z, α_D, Γ_D)[1] .+ ∇ℬ(z, α_N, Γ_N)[1]
    by = ∇ℬ(z, α_D, Γ_D)[2] .+ ∇ℬ(z, α_N, Γ_N)[2]
    f  = 𝒮_st⁻¹((ucm[1] .- y*ωcm .- bx, ucm[2] .+ x*ωcm .- by), Ω)

    fm = [real.(sum(f[1] .* Ω.dz)), real.(sum(f[2] .* Ω.dz))]
    fd = [real.(sum(x .* f[1] .* Ω.dz))  real.(sum(x .* f[2] .* Ω.dz));
          real.(sum(y .* f[1] .* Ω.dz))  real.(sum(y .* f[2] .* Ω.dz))]

    xD, xN = cos.(Γ_D.s), cos.(Γ_N.s)
    yD, yN = sin.(Γ_D.s), sin.(Γ_N.s)

    sm = [integral(xD .* α_D, Γ_D) + integral(xN .* α_N, Γ_N),
          integral(yD .* α_D, Γ_D) + integral(yN .* α_N, Γ_N)]
    sd = [integral(xD.*xD .* α_D, Γ_D) + integral(xN.*xN .* α_N, Γ_N)  integral(xD.*yD .* α_D, Γ_D) + integral(xN.*yN .* α_N, Γ_N);
          integral(yD.*xD .* α_D, Γ_D) + integral(yN.*xN .* α_N, Γ_N)  integral(yD.*yD .* α_D, Γ_D) + integral(yN.*yN .* α_N, Γ_N)]

    return (velocity=U, force_monopole=fm, force_dipole=fd, surf_monopole=sm, surf_dipole=sd)

end
