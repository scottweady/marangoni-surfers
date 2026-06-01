
using FastGaussQuadrature
using FFTW

"""
Discretized boundary curve for 2D Laplace BIM.

Fields
- `t`   : reference arc-length parameter (Chebyshev in (-1,1) for arcs, uniform in [0,2π) for circles)
- `s`   : physical arc-length parameter
- `γ`   : curve points γ(s) ∈ ℂ
- `γt`  : first derivative
- `γtt` : second derivative
- `w`   : reference quadrature weights (in `t`)
- `dΓ`  : physical quadrature weights (in `s`, i.e. |γ′(t)| dt)
- `V`   : self single-layer matrix (N × N)
- `V′`  : self adjoint double-layer matrix (N × N)
"""
struct Boundary
    open :: Bool
    t   :: Vector{Float64}
    s   :: Vector{Float64}
    γ   :: Vector{ComplexF64}
    γt  :: Vector{ComplexF64}
    γtt :: Vector{ComplexF64}
    w   :: Vector{Float64}
    dΓ  :: Vector{Float64}
    V   :: Matrix{Float64}
    V′  :: Matrix{Float64}
end

include("laplace.jl")

"""
Closed curve on physical parameter interval `dom` with N Gauss-Chebyshev nodes; uses Kress quadrature.

`γ`, `γs`, `γss` are functions of the physical parameter s.
The Chebyshev nodes t ∈ (-1,1) are mapped to s ∈ dom via s = a + (b-a)(t+1)/2.
`γt` and `γtt` are chain-ruled derivatives with respect to t.
"""
function curve(γ, γs, γss, N::Int)
    t, _ = _trigpts(N)
    s = t
    γ   = γ.(t)
    γt  = γs.(t)
    γtt = γss.(t)
    w   = fill(2π / N, N)
    dΓ = abs.(γt) .* w
    return Boundary(false, t, s, γ, γt, γtt, w, dΓ,
                    _laplace2d_single_layer_kress(t, γ, γt),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

circle(N::Int) = curve(s -> exp(1im * s), s -> 1im * exp(1im * s), s -> -exp(1im * s), N)

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
    dΓ  = abs.(γt) .* w
    return Boundary(true, t, s, γ, γt, γtt, w, dΓ,
                    _laplace2d_single_layer_bruno(t, γ, γt, w),
                    _laplace2d_double_layer_adjoint(t, γ, γt, γtt, w))
end

segment(N::Int; dom=[-1.0, 1.0]) = arc(s -> ComplexF64(s), s -> 1.0 + 0.0im, s -> 0.0im, N; dom=dom)

Base.length(Γ::Boundary) = length(Γ.s)
