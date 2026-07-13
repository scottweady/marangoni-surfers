
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

Nspan = 20:20:100

# Closed-curve solvers
include("test_circle.jl")
include("test_closedcurve.jl")

# Open-arc solvers
include("test_segment.jl")
include("test_openarc.jl")

# Mixed boundary condition solvers
include("test_mixed_openarc.jl")

# Marangoni solvers
include("../src/marangoni.jl")
include("test_surfer.jl")
include("test_ellipse_surfer.jl")