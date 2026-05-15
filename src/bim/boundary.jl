
using FastGaussQuadrature

"""
Discretized boundary curve for 2D Laplace BIM.

Fields
- `t`   : quadrature parameter nodes (Chebyshev in (-1,1) for arcs, uniform in [0,2π) for circles)
- `s`   : physical arc-length parameter (equals `t` for circles and segments)
- `γ`   : curve points γ(s) ∈ ℂ
- `γt`  : tangent dγ/dt (chain-ruled through the parameter map)
- `γtt` : second derivative d²γ/dt²
- `w`   : quadrature weights (in the `t` variable)
- `V`   : pre-built self single-layer matrix (N×N)
- `V′`  : pre-built self adjoint double-layer matrix (N×N)
"""
struct Boundary
    t   :: Vector{Float64}
    s   :: Vector{Float64}
    γ   :: Vector{ComplexF64}
    γt  :: Vector{ComplexF64}
    γtt :: Vector{ComplexF64}
    w   :: Vector{Float64}
    V   :: Matrix{Float64}
    V′  :: Matrix{Float64}
end

include("matrices.jl")

"Unit circle with N uniform nodes; uses Kress quadrature for the self single-layer."
function circle(N::Int)
    t, _ = trigpts(N)
    γ   = exp.(1im .* t)
    γt  = 1im .* exp.(1im .* t)
    γtt = -exp.(1im .* t)
    w   = fill(2π / N, N)
    return Boundary(t, t, γ, γt, γtt, w,
                    _laplace2d_single_layer_kress(t, γ, γt),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

"""
Closed curve on physical parameter interval `dom` with N Gauss-Chebyshev nodes; uses Kress quadrature.

`γ`, `γs`, `γss` are functions of the physical parameter s.
The Chebyshev nodes t ∈ (-1,1) are mapped to s ∈ dom via s = a + (b-a)(t+1)/2.
`γt` and `γtt` are chain-ruled derivatives with respect to t.
"""
function curve(γ, γs, γss, N::Int)
    t, _ = trigpts(N)
    s = t
    γ   = γ.(t)
    γt  = γs.(t)
    γtt = γss.(t)
    w   = fill(2π / N, N)
    return Boundary(t, s, γ, γt, γtt, w,
                    _laplace2d_single_layer_kress(t, γ, γt),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

"""
Straight segment on `dom` with N Gauss-Chebyshev nodes; uses Bruno quadrature for the self single-layer.

The density is assumed to have 1/√(1-t²) endpoint behavior.
"""
function segment(N::Int; dom=[-1.0, 1.0])
    t, w = gausschebyshevt(N)
    a, b = dom
    γ    = ComplexF64.(0.5 * (b - a) * t .+ 0.5 * (a + b))
    γt   = fill(ComplexF64(0.5 * (b - a)), N)
    γtt  = zeros(ComplexF64, N)
    return Boundary(t, t, γ, γt, γtt, w,
                    _laplace2d_single_layer_bruno(t, γ, γt, w),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

"""
Open arc on physical parameter interval `dom` with N Gauss-Chebyshev nodes; uses Bruno quadrature.

`γ`, `γs`, `γss` are functions of the physical parameter s.
The Chebyshev nodes t ∈ (-1,1) are mapped to s ∈ dom via s = a + (b-a)(t+1)/2.
`γt` and `γtt` are chain-ruled derivatives with respect to t.
"""
function arc(γ, γs, γss, N::Int; dom=[-1.0, 1.0])
    t, w = gausschebyshevt(N)
    a, b = dom
    s    = a .+ (b - a) .* (t .+ 1) ./ 2
    γ   = γ.(s)
    γt  = γs.(s)  .* ((b - a) / 2)
    γtt = γss.(s) .* ((b - a) / 2)^2
    return Boundary(t, s, γ, γt, γtt, w,
                    _laplace2d_single_layer_bruno(t, γ, γt, w),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

"Unit-circle arc on `dom`; convenience wrapper with γ(s) = exp(is)."
function arc(N::Int; dom=[-1.0, 1.0])
    return arc(s -> exp(1im * s), s -> 1im * exp(1im * s), s -> -exp(1im * s), N; dom=dom)
end

Base.length(Γ::Boundary) = length(Γ.s)
