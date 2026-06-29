
function _place(Γ::Boundary, xcm, θcm)
    eθ = exp(1im*θcm)
    return Boundary(Γ.quad, Γ.θ, Γ.t, Γ.s, xcm .+ eθ.*Γ.γ, eθ.*Γ.γt, eθ.*Γ.γtt, Γ.w, Γ.dΓ, Γ.V, Γ.K′)
end

function _jump_term(Γ::Boundary)
    if Γ.quad == :bruno
        return Diagonal(1 ./ (2 .* sqrt.(1 .- Γ.t.^2)))
    elseif Γ.quad == :theta
        return Diagonal(1 ./ (2 .* sin.(Γ.θ)))
    else
        return 0.5 * I(length(Γ))
    end
end

"Evaluate a domain field f at arc parameter values θ by tracing to the boundary and FFT-interpolating."
function _boundary_eval(f, Ω, s)
    vals  = trace(f, Ω)
    N     = length(vals)
    f̂     = fft(vals) ./ N
    freqs = fftfreq(N) * N
    return real.(sum(f̂[k] .* cis.(freqs[k] .* s) for k in 1:N))
end


"""
Assemble and solve the mixed Dirichlet/Neumann Laplace BIM system for a collection of boundaries.

`types[i]` is `:D` (Dirichlet, single-layer ansatz) or `:N` (Neumann, double-layer ansatz).
An extra unknown ϕ∞ (the far-field constant) is appended and fixed by enforcing ∫_all σ ds = 0.
Returns `(densities, ϕ∞)` where `densities[i]` is the layer density on boundary i.
"""
function laplace_solve(Γ::Vector{Boundary}, bc::Vector{Symbol}, g::Vector; decay=true, close_eval=false, close_eval_tol=0.01)

    Nparticles = length(Γ)
    offs  = [0; cumsum(length.(Γ))]
    Ntotal = offs[end]

    A   = zeros(Ntotal + 1, Ntotal + 1)
    rhs = zeros(Ntotal + 1)

    for i in 1:Nparticles
        ri = offs[i]+1 : offs[i+1]
        rhs[ri] = g[i]
        for j in 1:Nparticles
            cj = offs[j]+1 : offs[j+1]
            if i == j
                A[ri, cj] = bc[i] == :D ? Γ[i].V : Matrix(_jump_term(Γ[i])) + Γ[i].K′
            else
                A[ri, cj] = bc[i] == :D ? _laplace2d_single_layer(Γ[j], Γ[i]; close_eval=close_eval, close_eval_tol=close_eval_tol) : _laplace2d_double_layer_adjoint(Γ[j], Γ[i]; close_eval=close_eval, close_eval_tol=close_eval_tol)
            end
        end

        if bc[i] == :D
            A[ri, end] .= 1.0
        end

        A[end, ri] = Γ[i].dΓ

    end

    if decay
      σ_aug     = A \ rhs
      σ = [σ_aug[offs[i]+1:offs[i+1]] for i in 1:Nparticles]
      λ = σ_aug[end]
    else
      σ = A[1:Ntotal, 1:Ntotal] \ rhs[1:Ntotal]
      σ = [σ[offs[i]+1:offs[i+1]] for i in 1:Nparticles]
      λ = 0.0
    end

    return σ, λ

end

"""
Solve the two-region coupled pressure transmission problem:

    p1 = V_Γ[μ1] + f1          in Ω1 (interior of Γ)
    p2 = V_Γ[μ2] + V_Γext[τ] + f2   in Ω2 (between Γ and Γ_ext)

where f1 = -ξ1*(a12*δ²*H_Γ[σ] + a11/4*|x|²) and f2 analogously with ξ2, a21, a22.
H_Γ[σ] is the Helmholtz single-layer with k=1/δ, pre-solved so that H_Γ[σ]=1 on Γ.

Boundary conditions:
  - p1 = p2  on Γ
  - (1/ξ1) dp1/dn = (1/ξ2) dp2/dn  on Γ
  - dp2/dn = 0  on Γ_ext

Returns `(U1, U2)` where `Ui = -(1/ξi) dp/dn` evaluated on each boundary:
  - `U1` (length N):     normal velocity on Γ, limit from Ω1
  - `U2` (length N_ext): normal velocity on Γ_ext, ≈ 0 by construction
"""
function coupled_pressure_solve(Γ::Boundary, Γ_ext::Boundary;
                                 ξ1, ξ2, a11, a12, a21, a22, δ)

    N     = length(Γ)
    N_ext = length(Γ_ext)
    k     = 1 / δ

    # ── Helmholtz pre-solve: H[σ] = 1 on Γ ──────────────────────────────────
    H_self = _helmholtz2d_single_layer_kress(Γ.t, Γ.γ, Γ.γt; k=k)
    σ      = H_self \ ones(N)

    # ── Outward normals ───────────────────────────────────────────────────────
    n̂1    = (imag.(Γ.γt)     .- im.*real.(Γ.γt))     ./ abs.(Γ.γt)
    n̂_ext = (imag.(Γ_ext.γt) .- im.*real.(Γ_ext.γt)) ./ abs.(Γ_ext.γt)

    # ── Cross-boundary matrices ───────────────────────────────────────────────
    V_ext_on_Γ  = _laplace2d_single_layer(Γ_ext, Γ)
    K′_ext_on_Γ = _laplace2d_double_layer_adjoint(Γ_ext, Γ)
    K′_Γ_on_ext = _laplace2d_double_layer_adjoint(Γ, Γ_ext)

    K′_H_self = _helmholtz2d_double_layer_adjoint(Γ.t, Γ.γ, Γ.γt, Γ.γtt, Γ.w; k=k)

    K′_H_Γ_on_ext = zeros(N_ext, N)
    for i in 1:N_ext, j in 1:N
        r = Γ_ext.γ[i] - Γ.γ[j]
        ρ = abs(r)
        K′_H_Γ_on_ext[i, j] = -(k/(2π)) * besselk(1, k*ρ) * real(conj(n̂_ext[i]) * r) / ρ * Γ.dΓ[j]
    end

    # ── Forcing terms and their normal derivatives ────────────────────────────
    g1_Γ  = @. -ξ1*a12*δ^2 - ξ1*a11/4 * abs2(Γ.γ)      # f1 on Γ (H[σ]=1)
    g2_Γ  = @. -ξ2*a21*δ^2 - ξ2*a22/4 * abs2(Γ.γ)      # f2 on Γ

    # ∂f1/∂n1|_Γ from inside (- limit): ∂H/∂n|_- = (-½+K'_H)*σ
    dn1_H_in  = (-0.5 .* σ) .+ K′_H_self * σ
    dg1_dn1   = @. -ξ1*a12*δ^2 * dn1_H_in  - ξ1*a11/2 * real(conj(n̂1) * Γ.γ)

    # ∂f2/∂n1|_Γ from outside (+ limit): ∂H/∂n|_+ = (+½+K'_H)*σ
    dn1_H_out = (0.5 .* σ) .+ K′_H_self * σ
    dg2_dn1   = @. -ξ2*a21*δ^2 * dn1_H_out - ξ2*a22/2 * real(conj(n̂1) * Γ.γ)

    # ∂f2/∂n2|_Γext (smooth cross-boundary Helmholtz ADL)
    dn2_H_ext  = K′_H_Γ_on_ext * σ
    dg2_dn_ext = @. -ξ2*a21*δ^2 * dn2_H_ext - ξ2*a22/2 * real(conj(n̂_ext) * Γ_ext.γ)

    # ── Block system [μ1 (N), μ2 (N), τ (N_ext)] ─────────────────────────────
    K′_Γ   = Γ.K′
    K′_ext = Γ_ext.K′

    A = [
        Γ.V                                 -Γ.V                                   -V_ext_on_Γ;
        (1/ξ1).*(-0.5.*I(N) .+ K′_Γ)       -(1/ξ2).*(0.5.*I(N) .+ K′_Γ)         -(1/ξ2).*K′_ext_on_Γ;
        zeros(N_ext, N)                       K′_Γ_on_ext                           -0.5.*I(N_ext) .+ K′_ext
    ]
    b = [
        g2_Γ .- g1_Γ;
        (1/ξ2).*dg2_dn1 .- (1/ξ1).*dg1_dn1;
        -dg2_dn_ext
    ]

    sol = A \ b
    μ1  = sol[1:N]
    μ2  = sol[N+1:2N]
    τ   = sol[2N+1:end]

    # ── U = -(1/ξ) dp/dn on each boundary ────────────────────────────────────
    U1 = -(1/ξ1) .* ((-0.5 .* μ1 .+ K′_Γ * μ1) .+ dg1_dn1)
    U2 = -(1/ξ2) .* (K′_Γ_on_ext * μ2 .+ (-0.5 .* τ .+ K′_ext * τ) .+ dg2_dn_ext)

    return U1, U2

end

"Two-boundary convenience wrapper: Dirichlet arc Γ_D and Neumann arc Γ_N."
function laplace_solve(Γ_D::Boundary, Γ_N::Boundary, g_D, g_N)
    σ, ϕinf = laplace_solve([Γ_D, Γ_N], [:D, :N], [g_D, g_N])
    return σ[1], σ[2], ϕinf
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

"""
Solve the coupled Laplace–Helmholtz system on a closed boundary Γ:

    p = V[σ] + V[τ] - H[τ],   p|_Γ = 0
    φ = H[τ],                  φ|_Γ = 1

where V is the Laplace single-layer (kernel -log(r)/2π) and H is the Helmholtz single-layer
(kernel K₀(r)/2π).  Both boundary conditions give the decoupled system:

    H_mat τ = 1   (Helmholtz first-kind IE, enforces φ = 1 on Γ)
    V_mat σ = 1 - V_mat τ   (Laplace first-kind IE, enforces p = 0 on Γ)

Returns `(σ, τ, p_fun, φ_fun, velocity)` where the last three are callable at interior points:
  - `p_fun(ζ)`      : pressure field
  - `φ_fun(ζ)`      : Helmholtz field
  - `velocity(ζ)`   : (ux, uy) = -∇p at ζ (not valid on or near Γ)
"""
function helmholtz_laplace_solve(Γ::Boundary)

    N = length(Γ)
    H_mat = _helmholtz2d_single_layer_kress(Γ.t, Γ.γ, Γ.γt)

    # Enforces φ = 1 on Γ via Helmholtz first-kind IE: H_mat τ = 1.
    τ = H_mat \ ones(N)

    # Enforces p = 0 on Γ via Laplace first-kind IE: V_mat μ = H_mat τ = 1.
    # Augmented with a free constant c and gauge ∫μ dΓ = 0 to handle boundaries
    # of logarithmic capacity 1 (where V_mat has a zero eigenvalue).
    # Then p = V[μ] - H[τ] + c satisfies p = (1-c) - 1 + c = 0 on Γ.
    A_aug   = [Γ.V ones(N); Γ.dΓ' 0.0]
    sol     = A_aug \ [ones(N); 0.0]
    μ, c    = sol[1:N], sol[N+1]
    σ       = μ .- τ

    p_fun(ζ)  = 𝒱(ζ, μ, Γ) .- ℋ(ζ, τ, Γ) .+ c
    φ_fun(ζ)  = ℋ(ζ, τ, Γ)

    function velocity(ζ)
        # u = -∇p;  ∇c = 0 so c drops out
        vx, vy = ∇𝒱(ζ, μ, Γ)
        hx, hy = ∇ℋ(ζ, τ, Γ)
        return -(vx .- hx), -(vy .- hy)
    end

    return σ, τ, p_fun, φ_fun, velocity

end
