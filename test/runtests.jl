

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

# Smooth solvers
include("test_circle.jl")
include("test_curve.jl")
include("test_hypersingular.jl")

# Open-arc solvers
include("test_segment.jl")
include("test_arc.jl")

# Mixed boundary condition solvers
include("test_mixed_arc.jl")
