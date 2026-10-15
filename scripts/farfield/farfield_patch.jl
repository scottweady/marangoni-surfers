# Far-field model of a patch of surfers, started from the initial configuration of a
# boundary integral suspension run (scripts/suspension/suspension.jl) so the two can be
# compared directly, with the dipole-only model or, if `closure` is set, the monopole closure.
# Stops at `t_end` or on contact.
#
# Output: <outdir>/frame_<n>.dat, one row per particle
#   t  x  y  θ  ux  uy  ω

include("../../src/script_params.jl")
include("../../src/farfield.jl")

aspect_ratio = script_param(:aspect_ratio, 1.0)
t_end        = script_param(:t_end, 500.0)
closure      = script_param(:closure, false)
framedir     = script_param(:framedir, joinpath(@__DIR__, "..", "suspension", "output",
                                                "frames_$(aspect_ratio)"))
outdir       = script_param(:outdir, joinpath(@__DIR__, "output", "patch_$(aspect_ratio)"))

frame0 = [parse.(Float64, split(line)) for line in eachline(joinpath(framedir, "frame_000000.dat"))]
xcm = [row[2] + im * row[3] for row in frame0]
θcm = [row[4] for row in frame0]

U, η, β = farfield_constants(aspect_ratio)
@printf("a/b = %g   N = %d   U = %.6f   η = %.6f   β = %.6f\n",
        aspect_ratio, length(xcm), U, η, β)

isdir(outdir) && rm(outdir; recursive=true)
farfield_simulate(xcm, θcm, aspect_ratio; t_end=t_end, closure=closure, outdir=outdir)
