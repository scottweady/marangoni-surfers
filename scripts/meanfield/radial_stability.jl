# Linear stability of an expanding disk of radially polarized surfers in the mean-field model
# (meanfield.jl). The base state is a uniform disk (an annulus Rin < r < R, avoiding the
# orientation singularity at the center) with orientations von Mises distributed about the
# outward radial direction; it expands and is not steady, so stability is measured by the
# amplification of perturbations along the base trajectory. The perturbations are tilts
# δθ = cos(mφ) of the local orientation, with φ the polar angle about the disk center, one per
# azimuthal mode m; since the base state is axisymmetric, the modes evolve independently.
#
# Units: lengths are scaled by R and times by 1/(ηq), with q the number density, so the
# single remaining parameter besides the shape (β) and the spread (ζ) is Λ = ηqR/U, the ratio
# of the collective expansion speed to the self-propulsion speed.
#
# Output: <outdir>/radial_stability.dat, one row per saved time
#   t  R/R0  [log A_m  |α̂_m|  |ρ̂_m|] for each m
# where A_m is the amplification of the norm of the perturbation, α̂_m the amplitude of mode
# m of the tilt of the perturbed polarization, q₁·θ̂, relative to the base radial
# polarization q₀·r̂, and ρ̂_m that of the perturbed density relative to the base density,
# both multiplied by A_m and normalized by their initial tilt amplitude.

include(joinpath(@__DIR__, "meanfield.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "script_params.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "farfield.jl"))

aspect_ratio = script_param(:aspect_ratio, 4.0)
Λ            = script_param(:Λ, 10.0)
ζ            = script_param(:ζ, 25.0)
Rin          = script_param(:Rin, 0.1)
L            = script_param(:L, 10.0)
N            = script_param(:N, 256)
Nθ           = script_param(:Nθ, 32)
dT           = script_param(:dT, 1e-5)
dR           = script_param(:dR, 1e-5)
ms           = script_param(:ms, [0, 1, 2, 4, 8, 16])
t_end        = script_param(:t_end, 20.0)
tplt         = script_param(:tplt, 0.5)
outdir       = script_param(:outdir, joinpath(@__DIR__, "output", "radial_stability_$(aspect_ratio)"))

FFTW.set_num_threads(Threads.nthreads())

_, _, β = farfield_constants(aspect_ratio)
U, η = 1 / Λ, 1.0
G = SpectralGrid(N, N, Nθ; Lx=L, Ly=L)
xc, yc = G.Lx / 2, G.Ly / 2
polar(x, y) = atan(y - yc, x - xc)

base(δθ) = ic_patch(G; ρ0=1.0, R=1.0, Rin=Rin, width=2G.dx, ζ=ζ, δθ=δθ)
ψ0 = base((x, y) -> 0.0)

# tilt perturbation per unit amplitude, -cos(mφ) ∂ψ₀/∂θ, by a centered difference in the tilt
h = 1e-4
φ0s = [(base((x, y) -> h * cos(m * polar(x, y))) .- base((x, y) -> -h * cos(m * polar(x, y)))) ./ 2h
       for m in ms]

@printf("a/b = %g   β = %.4f   Λ = %g   U = %.4f   ζ = %g   dx = %.4f   modes = %s\n",
        aspect_ratio, β, Λ, U, ζ, G.dx, join(ms, ", "))

rr = [hypot(G.x[i] - xc, G.y[j] - yc) for i in 1:G.Nx, j in 1:G.Ny]
er = [cis(polar(G.x[i], G.y[j])) for i in 1:G.Nx, j in 1:G.Ny]

"Radius of the disk about the center containing a fraction `q` of the mass."
function mass_radius(ρ; q=0.99)
    k = sortperm(vec(rr))
    c = cumsum(vec(ρ)[k])
    return vec(rr)[k][findfirst(c .>= q * c[end])]
end

"Amplitudes of mode m of the tilt q₁·θ̂ relative to ∫q₀·r̂, and of ρ₁ relative to ∫ρ₀."
function mode_amplitudes(ψ, φ, m)
    ρ0, a1, a2 = moments(ψ, G)
    ρ1, b1, b2 = moments(φ, G)
    q0r = sum(real.(conj.(er) .* (a1 .+ im .* a2)))
    q1θ = imag.(conj.(er) .* (b1 .+ im .* b2))
    em = conj.(er) .^ m
    return abs(sum(q1θ .* em)) / q0r, abs(sum(ρ1 .* em)) / sum(ρ0)
end

mkpath(outdir)
io = open(joinpath(outdir, "radial_stability.dat"), "w")
println(io, "# t  R/R0  " * join(("logA_$m  alpha_$m  rho_$m" for m in ms), "  "))
R_init, α_init = Ref(NaN), zeros(length(ms))

function diagnostics(t, ψ, φs, logA)
    ρ0, _, _ = moments(ψ, G)
    R = mass_radius(ρ0)
    isnan(R_init[]) && (R_init[] = R)
    @printf(io, "%.4f  %.6f", t, R / R_init[])
    for (k, m) in enumerate(ms)
        α̂, ρ̂ = mode_amplitudes(ψ, φs[k], m)
        t == 0 && (α_init[k] = α̂)
        A = exp(logA[k]) / α_init[k]
        @printf(io, "  %+.6f  %.6e  %.6e", logA[k], A * α̂, A * ρ̂)
    end
    println(io)
    flush(io)
end

simulate_tangent(ψ0, φ0s, G; tf=t_end, U=U, β=β, η=η, dT=dT, dR=dR, tplt=tplt,
                 callback=diagnostics)
close(io)
