#=
Ewald + NUFFT summation of planar kernels from sources `yⱼ` with real weights `qⱼ` to targets
`xᵢ` (adapted from the self-mode Stokeslet Ewald sum in ProjectedSphericalHarmonics,
`suspensions/ewald.jl`),

    u(xᵢ) = Σⱼ G(xᵢ − yⱼ) qⱼ,   r = xᵢ − yⱼ ∈ R², ρ = |r|,

for a kernel G with K components, split as G = G_L + G_S with G_L smooth and G_S decaying
like exp(−ρ²/σ²) or erfc(ρ/σ):

 * Long range: G_L is smooth but need not decay. Multiplied by a box window W(r) = w(r₁)w(r₂)
   that equals 1 (to tolerance) for every separation in the bounding box of sources and
   targets, it has compact support, so its free-space convolution equals a periodic one on
   twice the support. Its Fourier multiplier comes from an FFT of samples; the sum is a
   type-1 NUFFT from the sources, the multiplier, and a type-2 NUFFT at the targets.
 * Short range: source-target pairs closer than r_c = cσ (erfc(c) ≈ tol) are summed directly
   on a cell list. Coincident pairs are skipped, and in self mode (targets = sources) their
   long-range part W(0) G_L(0) qᵢ is subtracted, so they contribute zero.

Kernels:
 * `_gradB_ewald`: G = r / (4πρ), the biharmonic gradient, with 1/ρ split as
   erf(ρ/σ)/ρ + erfc(ρ/σ)/ρ.
 * `_laplace_ewald_plan`: G = (log ρ, r/ρ²), the log potential and its gradient, with
   log ρ = [log ρ + E₁(ρ²/σ²)/2] − E₁(ρ²/σ²)/2.
=#

# ── Kernels ───────────────────────────────────────────────────────────────────

# erf(ρ/σ)/ρ, with its series near 0
@inline function _ewald_phiL(ρ, σ)
    x = ρ / σ
    x < 1e-3 && return (2 / (σ * sqrt(π))) * (1 - x^2 / 3 + x^4 / 10)
    return erf(x) / ρ
end

@inline _gradB_long(rx, ry, σ) = (f = _ewald_phiL(hypot(rx, ry), σ) / 4π; (f * rx, f * ry))
@inline _gradB_short(rx, ry, ρ, σ) = (f = erfc(ρ / σ) / (4π * ρ); (f * rx, f * ry))

# Ein(u) = Σₖ (−1)^{k+1} uᵏ / (k k!) = E₁(u) + log u + γ, for 0 ≤ u ≤ 1
@inline function _ein(u)
    s, t = 0.0, 1.0
    for k in 1:20
        t *= -u / k
        s -= t / k
    end
    return s
end

# log ρ + E₁(ρ²/σ²)/2 and its gradient r (1 − e^{−ρ²/σ²}) / ρ²
@inline function _laplace_long(rx, ry, σ)
    ρ2 = rx^2 + ry^2
    u  = ρ2 / σ^2
    if u < 1
        v = log(σ) + (_ein(u) - MathConstants.eulergamma) / 2
        f = u < 1e-8 ? (1 - u / 2) / σ^2 : -expm1(-u) / ρ2
    else
        v = log(ρ2) / 2 + expint(u) / 2
        f = -expm1(-u) / ρ2
    end
    return (v, f * rx, f * ry)
end

# −E₁(ρ²/σ²)/2 and its gradient r e^{−ρ²/σ²} / ρ²
@inline function _laplace_short(rx, ry, ρ, σ)
    u = ρ^2 / σ^2
    f = exp(-u) / ρ^2
    return (-expint(u) / 2, f * rx, f * ry)
end

# ── Neighbor search on a cell list ────────────────────────────────────────────

"Points sorted into square cells of side `h`: CSR `(start, idx)` over cells, and the cell map."
function _ewald_cells(X, lo, h)
    n  = size(X, 2)
    ci(x, d) = floor(Int, (x - lo[d]) / h) + 1
    nc = (ci(maximum(@view X[1, :]), 1), ci(maximum(@view X[2, :]), 2))
    cid = [ci(X[1, j], 1) + nc[1] * (ci(X[2, j], 2) - 1) for j in 1:n]
    counts = zeros(Int, nc[1] * nc[2])
    for c in cid; counts[c] += 1; end
    start = cumsum([1; counts])       # cell c holds idx[start[c]:start[c+1]-1]
    next  = start[1:end-1]
    idx   = Vector{Int}(undef, n)
    for j in 1:n
        idx[next[cid[j]]] = j
        next[cid[j]] += 1
    end
    return start, idx, nc, ci
end

# Number of sources within `rc` of each of up to 2000 strided targets, on average
function _ewald_mean_neighbors(Xs, Xt, lo, rc)
    start, idx, nc, ci = _ewald_cells(Xs, lo, rc)
    nt = size(Xt, 2)
    sample = nt <= 2000 ? (1:nt) : round.(Int, range(1, nt, length=2000))
    cnt = 0
    for i in sample
        x1, x2 = Xt[1, i], Xt[2, i]
        a, b = ci(x1, 1), ci(x2, 2)
        for cb in max(b - 1, 1):min(b + 1, nc[2]), ca in max(a - 1, 1):min(a + 1, nc[1])
            c = ca + nc[1] * (cb - 1)
            for p in start[c]:start[c + 1] - 1
                j = idx[p]
                cnt += ((x1 - Xs[1, j])^2 + (x2 - Xs[2, j])^2 < rc^2)
            end
        end
    end
    return cnt / length(sample)
end

# ── Plan and apply ────────────────────────────────────────────────────────────

"""
    _EwaldPlan

Ewald + NUFFT plan for a `K`-component planar kernel from fixed sources to fixed targets,
built by [`_ewald_plan`](@ref) and applied to weights by [`_ewald_apply`](@ref).
"""
mutable struct _EwaldPlan{K}
    ns::Int
    nt::Int
    M::NTuple{2,Int}
    Ghat::Array{ComplexF64,3}               # (M₁, M₂, K) multiplier, FFT mode order, scaled
    plan1::FINUFFT.finufft_plan{Float64}    # type 1 from the sources
    plan2::FINUFFT.finufft_plan{Float64}    # type 2 (K transforms) at the targets
    ptr::Vector{Int}                        # short range, CSR over targets:
    col::Vector{Int}                        #   sources idx col[ptr[i]:ptr[i+1]-1]
    val::Matrix{Float64}                    #   K × nnz kernel values
    self::Bool                              # targets = sources
    selfval::NTuple{K,Float64}              # own term removed in self mode
end

function _ewald_destroy!(E::_EwaldPlan)
    FINUFFT.finufft_destroy!(E.plan1)
    FINUFFT.finufft_destroy!(E.plan2)
    return nothing
end

"""
    _ewald_plan(long, short, Val(K), Xs, Xt; tol=1e-6, nneigh=100, self=false)

Ewald + NUFFT plan for a `K`-component planar kernel from the sources `Xs` (2 × ns) to the
targets `Xt` (2 × nt). `long(rx, ry, σ)` and `short(rx, ry, ρ, σ)` return the long- and
short-range parts of the kernel as `NTuple{K}`; `short` is called only for `0 < ρ < r_c`.
With `self=true`, `Xt` must be `Xs`, and the long-range part of each point's own term is
removed.

The splitting parameter is set so each target has about `nneigh` sources in its short-range
sum, trading the direct sum against the size of the Fourier grid.

The truncation errors scale with `Σ|qⱼ|`, which exceeds `max |u|` when the weights nearly
cancel, as for a surfactant layer density (by a factor of about 75 for a chain of 50 surfers).
The sum is therefore carried out at `tol / 10`, which keeps the error below `tol · max |u|`
in that case, at little extra cost.
"""
function _ewald_plan(long::FL, short::FS, ::Val{K}, Xs, Xt;
                     tol=1e-6, nneigh=100, self=false) where {FL, FS, K}

    ns, nt = size(Xs, 2), size(Xt, 2)
    lo = (min(minimum(@view Xs[1, :]), minimum(@view Xt[1, :])),
          min(minimum(@view Xs[2, :]), minimum(@view Xt[2, :])))
    hi = (max(maximum(@view Xs[1, :]), maximum(@view Xt[1, :])),
          max(maximum(@view Xs[2, :]), maximum(@view Xt[2, :])))
    c0 = (lo .+ hi) ./ 2
    L  = max.(hi .- lo, 1e-12)

    # Splitting parameter σ, from a few rescalings toward about `nneigh` short-range sources
    # per target (the count grows like σ² for spread sources, like σ along curves)
    tol = tol / 10
    c = max(erfcinv(tol), 3.0)
    σ = sqrt(nneigh * L[1] * L[2] / (π * ns)) / c
    for _ in 1:3
        m = _ewald_mean_neighbors(Xs, Xt, lo, c * σ)
        σ *= clamp(sqrt(nneigh / max(m, 1.0)), 0.25, 4.0)
    end
    rc = c * σ

    # Samples of the windowed long-range kernel on the periodic grid, and its multiplier.
    # w(t) = ½ erfc((|t| − L − rc)/σ) is 1 − O(tol) for |t| ≤ L and O(tol) beyond L + 2rc
    P = 2 .* (L .+ 2rc)
    h = π * σ / (2c)
    M = ntuple(d -> 2 * ceil(Int, P[d] / (2h)), 2)
    hx, hy = P[1] / M[1], P[2] / M[2]
    G = zeros(ComplexF64, M[1], M[2], K)
    Threads.@threads for b in 0:M[2]-1
        ry = (b < M[2] ÷ 2 ? b : b - M[2]) * hy
        wy = erfc((abs(ry) - L[2] - rc) / σ) / 2
        for a in 0:M[1]-1
            rx = (a < M[1] ÷ 2 ? a : a - M[1]) * hx
            wx = erfc((abs(rx) - L[1] - rc) / σ) / 2
            g  = long(rx, ry, σ)
            for k in 1:K
                G[a+1, b+1, k] = wx * wy * g[k]
            end
        end
    end
    FFTW.fft!(G, (1, 2))
    G .*= hx * hy / (P[1] * P[2])     # u = (1/P₁P₂) Σₘ Ĝ(kₘ) Q̂ₘ e^{i kₘ·x}, Ĝ ≈ h₁h₂ FFT

    # NUFFT plans at the sources (type 1) and targets (type 2)
    ntol = max(tol / 10, 100eps())
    nthr = Threads.nthreads()
    plan1 = FINUFFT.finufft_makeplan(1, [M[1], M[2]], -1, 1, ntol; modeord=1, nthreads=nthr)
    FINUFFT.finufft_setpts!(plan1, 2π .* (Xs[1, :] .- c0[1]) ./ P[1], 2π .* (Xs[2, :] .- c0[2]) ./ P[2])
    plan2 = FINUFFT.finufft_makeplan(2, [M[1], M[2]], +1, K, ntol; modeord=1, nthreads=nthr)
    FINUFFT.finufft_setpts!(plan2, 2π .* (Xt[1, :] .- c0[1]) ./ P[1], 2π .* (Xt[2, :] .- c0[2]) ./ P[2])

    # Own term in self mode: the windowed long-range kernel at r = 0
    W0 = (erfc((-L[1] - rc) / σ) / 2) * (erfc((-L[2] - rc) / σ) / 2)
    selfval = W0 .* long(0.0, 0.0, σ)

    # Short-range table: neighbors of each target (threaded), then CSR
    start, idx, nc, ci = _ewald_cells(Xs, lo, rc)
    rc2  = rc^2
    nbrs = Vector{Vector{Int}}(undef, nt)
    Threads.@threads for i in 1:nt
        x1, x2 = Xt[1, i], Xt[2, i]
        a, b = ci(x1, 1), ci(x2, 2)
        lst = Int[]
        for cb in max(b - 1, 1):min(b + 1, nc[2]), ca in max(a - 1, 1):min(a + 1, nc[1])
            cc = ca + nc[1] * (cb - 1)
            for p in start[cc]:start[cc + 1] - 1
                j = idx[p]
                ρ2 = (x1 - Xs[1, j])^2 + (x2 - Xs[2, j])^2
                (ρ2 < rc2 && ρ2 > 0) && push!(lst, j)
            end
        end
        nbrs[i] = lst
    end
    ptr = cumsum([1; length.(nbrs)])
    col = Vector{Int}(undef, ptr[end] - 1)
    val = Matrix{Float64}(undef, K, ptr[end] - 1)
    Threads.@threads for i in 1:nt
        x1, x2 = Xt[1, i], Xt[2, i]
        for (o, j) in enumerate(nbrs[i])
            p = ptr[i] + o - 1
            rx, ry = x1 - Xs[1, j], x2 - Xs[2, j]
            col[p] = j
            v = short(rx, ry, sqrt(rx^2 + ry^2), σ)
            for k in 1:K
                val[k, p] = v[k]
            end
        end
    end

    E = _EwaldPlan{K}(ns, nt, M, G, plan1, plan2, ptr, col, val, self, selfval)
    finalizer(_ewald_destroy!, E)
    return E

end

"""
    _ewald_apply(E, q)

Kernel sums of the plan `E` (see [`_ewald_plan`](@ref)) for the real source weights `q`.

# Returns
- `nt × K` matrix of the kernel sums at the targets.
"""
function _ewald_apply(E::_EwaldPlan{K}, q) where {K}

    Q = Array{ComplexF64}(undef, E.M[1], E.M[2])
    FINUFFT.finufft_exec!(E.plan1, ComplexF64.(q), Q)
    F = E.Ghat .* Q
    U = Matrix{ComplexF64}(undef, E.nt, K)
    FINUFFT.finufft_exec!(E.plan2, F, U)

    u = Matrix{Float64}(undef, E.nt, K)
    ptr, col, val, sv = E.ptr, E.col, E.val, E.selfval
    Threads.@threads for i in 1:E.nt
        s = ntuple(_ -> 0.0, Val(K))
        @inbounds for p in ptr[i]:ptr[i + 1] - 1
            qj = q[col[p]]
            s = s .+ ntuple(k -> val[k, p] * qj, Val(K))
        end
        @inbounds for k in 1:K
            u[i, k] = real(U[i, k]) + s[k] - (E.self ? sv[k] * q[i] : 0.0)
        end
    end
    return u

end

# ── Kernel sums ───────────────────────────────────────────────────────────────

"""
    _gradB_ewald(Xs, q, Xt; tol=1e-6, nneigh=100)

Ewald + NUFFT sum `u(xᵢ) = Σⱼ (xᵢ − yⱼ) qⱼ / (4π |xᵢ − yⱼ|)` from the sources `Xs` (2 × ns)
with real weights `q` to the targets `Xt` (2 × nt); see [`_ewald_plan`](@ref). Coincident
source-target pairs contribute zero.

# Returns
- `(ux, uy)`: vectors of length `nt`.
"""
function _gradB_ewald(Xs, q, Xt; tol=1e-6, nneigh=100)
    E = _ewald_plan(_gradB_long, _gradB_short, Val(2), Xs, Xt; tol=tol, nneigh=nneigh)
    u = _ewald_apply(E, q)
    _ewald_destroy!(E)
    return u[:, 1], u[:, 2]
end

"""
    _laplace_ewald_plan(X; tol=1e-6, nneigh=100)

Ewald + NUFFT plan (see [`_ewald_plan`](@ref)) for the sums `Σ_{j≠i} log|xᵢ − xⱼ| qⱼ` and
`Σ_{j≠i} (xᵢ − xⱼ) qⱼ / |xᵢ − xⱼ|²` over the distinct points `X` (2 × n). Applied to weights
`q` by [`_ewald_apply`](@ref), it gives the `n × 3` matrix with columns `(ϕ, gx, gy)`.
"""
_laplace_ewald_plan(X; tol=1e-6, nneigh=100) =
    _ewald_plan(_laplace_long, _laplace_short, Val(3), X, X; tol=tol, nneigh=nneigh, self=true)
