# Surfactant concentration for a random particle configuration, for the schematic figure
# (fig_schematic.m).
#
# Output: output/schematic_snapshot.dat
#   Line 1: # Nx=... Ny=... N=... a=... b=... k=... Nz=... Nb=...
#   Line 2: x values (Nx floats)
#   Line 3: y values (Ny floats)
#   Line 4: particle x centers (N floats)
#   Line 5: particle y centers (N floats)
#   Line 6: particle orientations θ (N floats)
#   Next Ny lines: ϕ matrix, Ny rows × Nx columns (NaN inside particles)
#   Next line: zoom x values (Nz floats), centered on particle k
#   Next line: zoom y values (Nz floats)
#   Next Nz lines: zoom ϕ matrix, Nz rows × Nz columns (NaN inside particles)
#   Next line: boundary parameter s of particle k (Nb floats)
#   Next line: ϕ on the boundary of particle k (Nb floats)

using Random, Statistics

using MarangoniSurfers
using MarangoniSurfers: _find_preimage
using Printf

# -- Geometry ------------------------------------------------------------------

aspect_ratio = 2.0
a = sqrt(aspect_ratio)
b = 1.0 / a

s₀   = π/2
g(s) = -cos(s)

Ns = 64

γ(s) = a*cos(s) + im*b*sin(s)
Γ_D, Γ_N = surfer_boundary(a, b, s₀, Ns; quad=:mixed)   # 𝒱_close needs Gauss-Chebyshev arcs

# -- Random initial configuration (rejection sampling) -------------------------

Random.seed!(3)
Nparticles = 14
L     = 12.0          # half-width of the box
d_min = 2a + 3.0      # minimum center-to-center distance

xcm = ComplexF64[]
while length(xcm) < Nparticles
    z = L * (2rand() - 1) + im * L * (2rand() - 1)
    all(abs(z - w) > d_min for w in xcm) && push!(xcm, z)
end
θcm = 2π .* rand(Nparticles)

# -- Surfactant solve ----------------------------------------------------------

σ, ϕinf, Γ = surfactant_solve(xcm, θcm, Γ_D, Γ_N, g)

inside(z) = any(abs2(real(w) / a) + abs2(imag(w) / b) < 1.0
                for w in (exp(-im * θcm[k]) * (z - xcm[k]) for k in 1:Nparticles))

"Concentration at `z`, using singularity-swap close evaluation for arcs within `rclose`."
function ϕ_eval(z; rclose=0.75)
    val = ϕinf
    for (σk, Γk) in zip(σ, Γ)
        d, j = findmin(abs.(z .- Γk.γ))
        if d < rclose
            t0, ok = _find_preimage(z, Γk; tol=1e-12, t_init=Γk.t[j])
            val += ok ? 𝒱_close(z, σk, Γk; t_tgt=t0) : 𝒱(z, σk, Γk)
        else
            val += 𝒱(z, σk, Γk)
        end
    end
    return real(val)
end

field(xs, ys) = [inside(x + im*y) ? NaN : ϕ_eval(x + im*y) for y in ys, x in xs]

# Square plotting window around the configuration
c    = mean(xcm)
hw   = maximum(max.(abs.(real.(xcm .- c)), abs.(imag.(xcm .- c)))) + 4.0
Nplt = 400
xplt = range(real(c) - hw, real(c) + hw, length=Nplt)
yplt = range(imag(c) - hw, imag(c) + hw, length=Nplt)
ϕ = field(xplt, yplt)

# Zoom on the particle closest to the center of the configuration
k     = argmin(abs.(xcm .- c))
hz    = 2.5
Nz    = 300
xz    = range(real(xcm[k]) - hz, real(xcm[k]) + hz, length=Nz)
yz    = range(imag(xcm[k]) - hz, imag(xcm[k]) + hz, length=Nz)
ϕz = field(xz, yz)

# Boundary concentration of particle k (midpoints avoid the arc endpoints)
Nb = 512
sb = range(-π/2, 3π/2, length=Nb+1)
sb = (sb[1:end-1] .+ sb[2:end]) ./ 2
ϕb = [ϕ_eval(xcm[k] + exp(im * θcm[k]) * γ(s)) for s in sb]


# -- Write results -------------------------------------------------------------

write_row(io, v) = println(io, join([@sprintf("%.10e", x) for x in v], " "))
write_matrix(io, A) = foreach(i -> write_row(io, A[i, :]), axes(A, 1))

outdir = script_param(:outdir, joinpath(@__DIR__, "output"))
mkpath(outdir)
outfile = joinpath(outdir, "schematic_snapshot.dat")
open(outfile, "w") do io
    @printf(io, "# Nx=%d Ny=%d N=%d a=%.10f b=%.10f k=%d Nz=%d Nb=%d\n",
            Nplt, Nplt, Nparticles, a, b, k, Nz, Nb)
    write_row(io, xplt)
    write_row(io, yplt)
    write_row(io, real.(xcm))
    write_row(io, imag.(xcm))
    write_row(io, θcm)
    write_matrix(io, ϕ)
    write_row(io, xz)
    write_row(io, yz)
    write_matrix(io, ϕz)
    write_row(io, sb)
    write_row(io, ϕb)
end

println("Wrote ", relpath(outfile, @__DIR__))
