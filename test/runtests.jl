
using Plots
using Printf

include("../src/bim/boundary.jl")
include("../src/bim/operators.jl")
include("../src/bim/solvers.jl")

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
include("test_segment.jl")
include("test_circle.jl")
include("test_arc.jl")
include("test_curve.jl")

include("test_close_eval.jl")
include("test_mixed_arc.jl")

include("test_surfer.jl")
include("test_ellipse_surfer.jl")
