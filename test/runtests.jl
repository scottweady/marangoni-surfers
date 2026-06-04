

include("../src/marangoni.jl")

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

Nspan = 10:10:100

# Smooth solvers
include("test_circle.jl")
include("test_curve.jl")

# Open-arc solvers
include("test_segment.jl")
include("test_arc.jl")

# Mixed boundary condition solvers
include("test_mixed_arc.jl")

# Marangoni solvers
include("test_surfer.jl")
include("test_ellipse_surfer.jl")
