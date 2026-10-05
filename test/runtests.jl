# Convergence tests for the boundary integral solvers and the Marangoni solve. Each test
# prints its error against an exact or high-resolution reference over `Nspan`.

include("../src/bie/bie.jl")

function printerror(N, err)
    color = if err < 1e-12
        :green
    elseif err > 1e-4
        :red
    else
        :yellow
    end
    printstyled("N = $N, error = $err\n", color=color)
end

Nspan = 4:4:24

# Closed-curve solvers
include("test_circle.jl")
include("test_closedcurve.jl")

# Open-arc solvers
include("test_segment.jl")
include("test_openarc.jl")
include("test_close_eval.jl")

# Mixed boundary condition solvers
include("test_mixed_openarc.jl")
include("test_mixed_exterior.jl")

# Marangoni solvers
include("../src/marangoni.jl")
include("test_surfer.jl")
include("test_ellipse_surfer.jl")
