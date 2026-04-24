
using FastGaussQuadrature

struct Boundary
    s   :: Vector{Float64}
    x   :: Vector{ComplexF64}
    xs  :: Vector{ComplexF64}
    xss :: Vector{ComplexF64}
    V   :: Matrix{Float64}
    V′  :: Matrix{Float64}
end

include("matrices.jl")

function circle(N::Int)
    s, _ = trigpts(N)
    x   = exp.(1im .* s)
    xs  = 1im .* exp.(1im .* s)
    xss = -exp.(1im .* s)
    return Boundary(s, x, xs, xss,
                    _laplace2d_single_layer_matrix(s, x, xs),
                    _laplace2d_adjoint_double_layer_matrix(s, x, xs, xss))
end

function segment(N::Int)
    j   = 1:N
    θ   = π .* (2j .- 1) ./ (2N)
    s   = cos.(θ)
    x   = ComplexF64.(s)
    xs  = ones(ComplexF64, N)
    xss = zeros(ComplexF64, N)
    return Boundary(s, x, xs, xss,
                    _laplace2d_single_layer_matrix_arc(s, x, xs),
                    _laplace2d_adjoint_double_layer_matrix_arc(s, x, xs, xss))
end

function arc(N::Int; dom =[θ_a, θ_b])
    θ_a, θ_b = dom
    j    = 1:N
    s    = cos.(π .* (2j .- 1) ./ (2N))   # Chebyshev-Gauss nodes on (-1,1)
    θ    = θ_a .+ (θ_b - θ_a) .* (s .+ 1) ./ 2
    dθds = (θ_b - θ_a) / 2
    x    = exp.(1im .* θ)
    xs   = 1im .* exp.(1im .* θ) .* dθds
    xss  = -exp.(1im .* θ) .* dθds^2
    return Boundary(s, x, xs, xss,
                    _laplace2d_single_layer_matrix_arc(s, x, xs),
                    _laplace2d_adjoint_double_layer_matrix_arc(s, x, xs, xss))
end

Base.length(Γ::Boundary) = length(Γ.s)