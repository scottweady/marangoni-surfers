
"""
    Boundary

Discretized boundary curve for 2D Laplace BIM.

# Fields
- `quad::Symbol`: quadrature type — `:closed` (closed), `:open_weight` (open, `1/√(1-t²)` density),
  `:open` (open, smooth density).
- `t::Vector{Float64}`: reference parameter (Chebyshev nodes in `(-1,1)` for arcs, uniform in
  `[0,2π)` for circles).
- `s::Vector{Float64}`: physical arc-length parameter.
- `γ::Vector{ComplexF64}`: curve points `γ(s) ∈ ℂ`.
- `γt::Vector{ComplexF64}`: first derivative with respect to `t`.
- `γtt::Vector{ComplexF64}`: second derivative with respect to `t`.
- `w::Vector{Float64}`: reference quadrature weights (in `t`).
- `dΓ::Vector{Float64}`: physical quadrature weights (in `s`, i.e. `|γ′(t)| dt`).
- `V::Matrix{Float64}`: self single-layer matrix (`N × N`).
- `K′::Matrix{Float64}`: self adjoint double-layer matrix (`N × N`).
"""
struct Boundary
    quad :: Symbol
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
    closedcurve(γ, γs, γss, N; dom=[0.0, 2π])

Closed curve on physical parameter interval `dom` with `N` uniform trigonometric nodes;
uses Kress quadrature.

The Chebyshev nodes `t ∈ (-1,1)` are mapped to `s ∈ dom` via `s = a + (b-a)(t+1)/2`, and
`γt`, `γtt` are the corresponding chain-ruled derivatives with respect to `t`.

# Arguments
- `γ`, `γs`, `γss`: functions of the physical parameter `s` giving the curve, its first
  derivative, and its second derivative.
- `N::Int`: number of quadrature nodes.
- `dom`: physical parameter interval (default `[0.0, 2π]`).

# Returns
- `Boundary` with `quad = :closed`.
"""
function closedcurve(γ, γs, γss, N::Int; dom=[0.0, 2π])
    t, _ = _trigpts(N; dom=dom)
    s = t
    γ   = γ.(t)
    γt  = γs.(t)
    γtt = γss.(t)
    w   = fill(2π / N, N)
    dΓ = abs.(γt) .* w
    V = _laplace2d_single_layer_closed(t, γ, γt)
    K′ = _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    return Boundary(:closed, t, s, γ, γt, γtt, w, dΓ, V, K′)
end

"""
    circle(N)

Unit circle with `N` quadrature nodes, via [`closedcurve`](@ref).

# Arguments
- `N::Int`: number of quadrature nodes.

# Returns
- `Boundary` for the unit circle.
"""
circle(N::Int) = closedcurve(s -> exp(1im * s), s -> 1im * exp(1im * s), s -> -exp(1im * s), N)

"""
    openarc(γ_fun, γs_fun, γss_fun, N; dom=[-1.0, 1.0], quad=:open_weight)

Open arc on physical parameter interval `dom` with `N` quadrature nodes.

The nodes `t ∈ (-1,1)` are mapped to `s ∈ dom` via `s = a + (b-a)(t+1)/2`, and `γt`, `γtt`
are the corresponding chain-ruled derivatives with respect to `t`.

# Arguments
- `γ_fun`, `γs_fun`, `γss_fun`: functions of the physical parameter `s` giving the curve,
  its first derivative, and its second derivative.
- `N::Int`: number of quadrature nodes.
- `dom`: physical parameter interval (default `[-1.0, 1.0]`).
- `quad`: quadrature scheme — `:open_weight` (default) uses Gauss-Chebyshev nodes with 
  log-splitting for `1/√(1-t²)` densities; `:open` uses Fejér first-rule nodes
  with singularity subtraction for smooth densities.

# Returns
- `Boundary` with the requested `quad` type.
"""
function openarc(γ_fun, γs_fun, γss_fun, N::Int; dom=[-1.0, 1.0], quad=:open_weight)

    if quad == :open
        t, w = _fejerpts(N)
        _single_layer = _laplace2d_single_layer_open
    elseif quad == :open_weight
        t, w = _chebpts(N)
        _single_layer = _laplace2d_single_layer_open_weight
    else
        error("Unknown quadrature type: $quad")
    end

    a, b = dom
    s    = a .+ (b - a) .* (t .+ 1) ./ 2
    γ   = γ_fun.(s)
    γt  = γs_fun.(s)  .* ((b - a) / 2)
    γtt = γss_fun.(s) .* ((b - a) / 2)^2
    dΓ  = abs.(γt) .* w
    V = _single_layer(t, γ, γt, w)
    K′  = _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w)
    return Boundary(quad, t, s, γ, γt, γtt, w, dΓ, V, K′)
end

"""
    segment(N; dom=[-1.0, 1.0], quad=:open_weight)

Straight line segment on `dom` with `N` quadrature nodes, via [`openarc`](@ref).

# Arguments
- `N::Int`: number of quadrature nodes.
- `dom`: physical parameter interval (default `[-1.0, 1.0]`).
- `quad`: quadrature scheme, passed through to [`openarc`](@ref) (default `:open_weight`).

# Returns
- `Boundary` for the line segment.
"""
segment(N::Int; dom=[-1.0, 1.0], quad=:open_weight) = openarc(s -> ComplexF64(s), s -> 1.0 + 0.0im, s -> 0.0im, N; dom=dom, quad=quad)

"""
    length(Γ)

Number of quadrature nodes on `Γ`.
"""
Base.length(Γ::Boundary) = length(Γ.s)
