# Convergence tests for the boundary integral solvers and the Marangoni solve. Each test
# computes its error against an exact or high-resolution reference and asserts on the
# finest resolution. By default only the asserted resolutions are run; with `--verbose`
# each test runs its full span of resolutions and prints the convergence table.
#
#     julia --project=.. runtests.jl              # quick pass/fail summary
#     julia --project=.. runtests.jl --verbose    # full convergence study

using Test

include("../src/bie/bie.jl")

const VERBOSE = "--verbose" in ARGS

section(title)    = VERBOSE && println("\n", title)
subsection(title) = VERBOSE && println("  ", title)

# Resolutions to run: the full `span` for the convergence study, otherwise only `quick`
sweep(span, quick=[last(span)]) = VERBOSE ? span : quick

# Setting for the convergence study, or a cheaper one that still passes for the quick run
setting(full, quick) = VERBOSE ? full : quick

# `label` names the resolution: "N" for the boundary resolution, "M" for the body (Stokes)
# expansion order
function printerror(n, err; label="N")
    color = if err < 1e-12
        :green
    elseif err > 1e-4
        :red
    else
        :yellow
    end
    printstyled("$label = $n, error = $err\n", color=color)
end

# Error bound for every assertion, well above the converged errors and well below any real
# failure
const TOL = 1e-6

# Print the convergence table when verbose, and test the error at the finest resolution
# (or at every resolution, for checks that hold exactly at any resolution)
function check_errors(span, errs; label="N", every=false)
    VERBOSE && foreach((n, err) -> printerror(n, err; label=label), span, errs)
    if every
        @test maximum(errs) < TOL
    else
        @test last(errs) < TOL
    end
end

@testset "Marangoni" begin

    # Closed-curve solvers
    include("test_closedcurve.jl")

    # Open-arc solvers
    include("test_segment.jl")
    include("test_openarc.jl")
    include("test_close_eval.jl")

    # Mixed boundary condition solvers
    include("test_mixed_exterior.jl")

    # Marangoni solvers
    include("../src/marangoni.jl")
    include("test_surfer.jl")
    include("test_multiparticle.jl")

end
