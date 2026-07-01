
"""
Discretized boundary curve for 2D Laplace BIM.

Fields
- `quad` : quadrature type — `:kress` (closed), `:bruno` (open, 1/√(1-t²) density), `:smooth` (open, smooth density)
- `t`   : reference arc-length parameter (Chebyshev in (-1,1) for arcs, uniform in [0,2π) for circles)
- `s`   : physical arc-length parameter
- `γ`   : curve points γ(s) ∈ ℂ
- `γt`  : first derivative
- `γtt` : second derivative
- `w`   : reference quadrature weights (in `t`)
- `dΓ`  : physical quadrature weights (in `s`, i.e. |γ′(t)| dt)
- `V`   : self single-layer matrix (N × N)
- `K′`  : self adjoint double-layer matrix (N × N)
"""
struct Boundary
    quad :: Symbol
    θ   :: Vector{Float64}
    t   :: Vector{Float64}
    s   :: Vector{Float64}
    γ   :: Vector{ComplexF64}
    γt  :: Vector{ComplexF64}
    γtt :: Vector{ComplexF64}
    w   :: Vector{Float64}
    dΓ  :: Vector{Float64}
    V   :: Matrix{Float64}
    K′  :: Matrix{Float64}
end

"""
Closed curve on physical parameter interval `dom` with N Gauss-Chebyshev nodes; uses Kress quadrature.

`γ`, `γs`, `γss` are functions of the physical parameter s.
The Chebyshev nodes t ∈ (-1,1) are mapped to s ∈ dom via s = a + (b-a)(t+1)/2.
`γt` and `γtt` are chain-ruled derivatives with respect to t.
"""
function curve(γ, γs, γss, N::Int; dom=[0.0, 2π])
    t, _ = _trigpts(N; dom=dom)
    s = t
    γ   = γ.(t)
    γt  = γs.(t)
    γtt = γss.(t)
    w   = fill(2π / N, N)
    dΓ = abs.(γt) .* w
    V = _laplace2d_single_layer_kress(t, γ, γt)
    K′ = _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    return Boundary(:kress, t, t, s, γ, γt, γtt, w, dΓ, V, K′)
end

circle(N::Int) = curve(s -> exp(1im * s), s -> 1im * exp(1im * s), s -> -exp(1im * s), N)

"""
Open arc on physical parameter interval `dom` with N quadrature nodes.

`γ`, `γs`, `γss` are functions of the physical parameter s.
The nodes t ∈ (-1,1) are mapped to s ∈ dom via s = a + (b-a)(t+1)/2.
`γt` and `γtt` are chain-ruled derivatives with respect to t.

`quad` selects the quadrature scheme:
- `:bruno`  (default) — Gauss-Chebyshev nodes, Bruno (2021) log-splitting for 1/√(1-t²) densities
- `:smooth` — Fejér first-rule nodes, singularity subtraction for smooth densities
"""
function arc(γ_fun, γs_fun, γss_fun, N::Int; dom=[-1.0, 1.0], quad=:bruno)

    if quad == :smooth
        t, w = _fejerpts(N)
        θ = t
        _single_layer = _laplace2d_single_layer_smooth
    elseif quad == :bruno
        t, w = _chebpts(N)
        θ = t
        _single_layer = _laplace2d_single_layer_bruno
    elseif quad == :theta
        θ, w = _fejerpts(N; dom=[0.0, π])
        t = cos.(θ)
    else
        error("Unknown quadrature type: $quad")
    end

    a, b = dom
    s    = a .+ (b - a) .* (t .+ 1) ./ 2
    γ   = γ_fun.(s)
    γt  = γs_fun.(s)  .* ((b - a) / 2)
    γtt = γss_fun.(s) .* ((b - a) / 2)^2
    dΓ  = abs.(γt) .* w

    if quad == :theta
        t̃, w̃ = _chebpts(N)
        s̃ = a .+ (b - a) .* (t̃ .+ 1) ./ 2
        γ̃ = γ_fun.(s̃)
        γ̃t = γs_fun.(s̃) .* ((b - a) / 2)
        V = _laplace2d_single_layer_theta(θ, γ, γt, w, t̃, γ̃, γ̃t, w̃)
    else
        V = _single_layer(t, γ, γt, w)
    end

    K′  = _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    return Boundary(quad, θ, t, s, γ, γt, γtt, w, dΓ, V, K′)
end

segment(N::Int; dom=[-1.0, 1.0], quad=:bruno) = arc(s -> ComplexF64(s), s -> 1.0 + 0.0im, s -> 0.0im, N; dom=dom, quad=quad)

Base.length(Γ::Boundary) = length(Γ.s)
