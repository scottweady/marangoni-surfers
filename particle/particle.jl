
using LinearAlgebra
using Plots

function uinf(xᵢ, xⱼ, pᵢ, pⱼ)
    u = zeros(Float64, 2)
    rvec = xᵢ - xⱼ
    r = norm(rvec)
    r̂ = rvec ./ r
    M = (I(2) .- r̂ * r̂') ./ (4π * r)
    return M * pⱼ
end

function ωinf(xᵢ, xⱼ, pᵢ, pⱼ)
    rvec = xᵢ - xⱼ
    r = norm(rvec)
    r̂ = rvec ./ r
    M = ((I(2) .- 3 * r̂ * r̂') .* (r̂' * pⱼ) .+ r̂ * pⱼ' .+ pⱼ * r̂') ./ (4π * r.^2)
    qᵢ = [-pᵢ[2], pᵢ[1]]
    return -(qᵢ' * M * pᵢ)
end

function rhs(x, θ, β, λ, uS)
    N = size(x, 1)
    p = [cos.(θ) sin.(θ)]
    u = uS * p
    ω = zeros(N)
    for nᵢ = 1:N
        for nⱼ = 1:N
            if nᵢ != nⱼ
                u[nᵢ,:] .+= β * uinf(x[nᵢ,:], x[nⱼ,:], p[nᵢ,:], p[nⱼ,:])
                ω[nᵢ]   += β * λ * ωinf(x[nᵢ,:], x[nⱼ,:], p[nᵢ,:], p[nⱼ,:])
            end
        end
    end
    return u, ω
end

function simulate(x, θ, β, λ, uS; dt=0.01, tf=100.0)

    Nt = tf ÷ dt
    nt = 0

    while nt < Nt

        k1u, k1ω = rhs(x,           θ,           β, λ, uS)
        k2u, k2ω = rhs(x .+ dt/2 * k1u, θ .+ dt/2 * k1ω, β, λ, uS)
        k3u, k3ω = rhs(x .+ dt/2 * k2u, θ .+ dt/2 * k2ω, β, λ, uS)
        k4u, k4ω = rhs(x .+ dt   * k3u, θ .+ dt   * k3ω, β, λ, uS)

        x .+= dt/6 * (k1u .+ 2k2u .+ 2k3u .+ k4u)
        θ .+= dt/6 * (k1ω .+ 2k2ω .+ 2k3ω .+ k4ω)
        display((k1u .+ 2k2u .+ 2k3u .+ k4u))
        nt += 1

        if nt % 100 == 0
            plt = scatter(vec(x[:,1]),vec(x[:,2]); aspect_ratio=:equal)
            xlims!(-4.0, 4.0)
            ylims!(-4.0, 4.0)
            display(plt)
        end

    end

end

dt = 0.01
tf = 10000.0

# a = 0.25
# b = 1.0/a
# β = 0.048896
# uS = 0.003081

a = 1.0
b = 1.0
β =  1.570796
uS = 0.098175

# a = 2.0
# b = 0.5
# β =  1.689178
# uS = 0.128646

λ = (a^2 - b^2) / (a^2 + b^2)

R = β/(4π * uS) * (1 - λ/2)

display(R)
x = [-R/2 0.0;
      R/2 0.0]
θ = [-π/2, π/2]

simulate(x, θ, β, λ, uS, dt=dt, tf=tf)
