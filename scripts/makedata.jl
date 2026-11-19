# Generate the data for the paper figures (figs/makefigs.m), writing into data/, or into
# data-DD-MM-YYYY (today's date) if data/ already exists.
#
#     julia --project=scripts scripts/makedata.jl [single] [pairs] [suspension] [tandem] [resume] [data=<dir>]
#
# With no groups, all groups are run. Each generator runs in its own headless Julia process
# with the settings used in the figures, passed as `name=value` overrides (see
# `script_param`); run on their own, the generators keep their exploratory defaults.
#
#   single      single-particle and pair-gap convergence, monopole sweep, schematic
#   pairs       pair fields and pair trajectories
#   suspension  patch expansion runs, their far-field counterparts, and the expansion snapshot
#   tandem      tandem chain runs, their linear stability, and the chain snapshots
#
# Options:
#   resume      skip runs whose output is already complete; a run that was interrupted
#               (a trajectory or suspension that did not reach its end time) is redone
#   data=<dir>  write into <dir> (relative to the repository root, or absolute), e.g. to
#               resume into an existing data-DD-MM-YYYY directory

using Dates
using Printf

const ROOT = dirname(@__DIR__)

const OPTIONS = filter(arg -> arg == "resume" || startswith(arg, "data="), ARGS)
const GROUPS  = filter(arg -> !(arg in OPTIONS), ARGS)
const RESUME  = "resume" in OPTIONS

const DATA = let given = filter(arg -> startswith(arg, "data="), OPTIONS)
    if !isempty(given)
        dir = last(given)[length("data=")+1:end]
        isabspath(dir) ? dir : joinpath(ROOT, dir)
    elseif isdir(joinpath(ROOT, "data"))
        joinpath(ROOT, "data-" * Dates.format(today(), "dd-mm-yyyy"))
    else
        joinpath(ROOT, "data")
    end
end

"""
Run `script` (relative to the repository root) with `name=value` overrides, unless resuming
and its output is already complete (`done`).
"""
function run_script(script, params...; done=false)
    args = [string(k, "=", v isa AbstractString ? v : repr(v)) for (k, v) in params]
    if RESUME && done
        println("\n== skip (complete) ", script, " ", join(args, " "))
        return
    end
    println("\n== ", script, " ", join(args, " "))
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "scripts")) $(joinpath(ROOT, script)) $args`
    t = @elapsed run(pipeline(addenv(cmd, "GKSwstype" => "100"); stdin=devnull))
    println("== done in ", round(t / 60, digits=1), " min")
end

# -- Completeness checks for resume --------------------------------------------

"Trajectory file written by pair.jl, for aspect ratio `aspect_ratio` and separation `Δx`."
function trajectory_file(dir, aspect_ratio, Δx)
    a = sqrt(aspect_ratio)
    return joinpath(dir, @sprintf("a%.4f_b%.4f_d%.1f.dat", a, 1/a, Δx))
end

"Whether a pair trajectory reached `t_max` or stopped on contact (gap below `stop_dist`)."
function trajectory_done(file, t_max; stop_dist=1e-2)
    isfile(file) || return false
    lines = filter(!isempty, readlines(file))
    isempty(lines) && return false
    row = parse.(Float64, split(last(lines)))
    return row[1] >= t_max - 1e-6 || row[end] < stop_dist
end

"Whether the frames written by suspension.jl in `dir` reach `t_end`."
function frames_done(dir, t_end)
    isdir(dir) || return false
    files = filter(f -> startswith(f, "frame_") && endswith(f, ".dat"), readdir(dir))
    isempty(files) && return false
    t_last = maximum(parse(Float64, split(readline(joinpath(dir, f)))[1]) for f in files)
    return t_last >= t_end - 1e-6
end

# -- Groups --------------------------------------------------------------------

function single()
    convdir = joinpath(DATA, "convergence")
    run_script("scripts/convergence/convergence.jl", :plot => false, :outdir => convdir;
               done=all(isfile, joinpath.(convdir, ["convergence_expansion.dat", "convergence_boundary.dat"])))
    run_script("scripts/convergence/pair_gap_convergence.jl", :plot => false, :outdir => convdir;
               done=isfile(joinpath(convdir, "pair_gap_expansion.dat")))
    run_script("scripts/farfield/monopole_vs_separation.jl",
               :outdir => joinpath(DATA, "multipole");
               done=isfile(joinpath(DATA, "multipole", "monopole_vs_separation.dat")))
    run_script("figs/schematic_snapshot.jl",
               :outdir => joinpath(DATA, "schematic");
               done=isfile(joinpath(DATA, "schematic", "schematic_snapshot.dat")))
end

function pairs()
    field_ratios = [0.25, 1.0, 4.0]
    for orient in (:parallel, :antiparallel)
        dir = joinpath(DATA, string(orient))
        files = [joinpath(dir, @sprintf("fields_a%.4f_b%.4f.dat", sqrt(r), 1/sqrt(r))) for r in field_ratios]
        run_script("scripts/pair/pair_fields.jl", :orient => orient,
                   :aspect_ratios => field_ratios, :outdir => dir; done=all(isfile, files))
    end

    # fig_pair_parallel plots every trajectory file in data/parallel, so clear stale ones
    # (unless resuming, when completed trajectories are kept)
    pardir = mkpath(joinpath(DATA, "parallel"))
    RESUME || foreach(f -> rm(joinpath(pardir, f)),
                      filter(f -> startswith(f, "a") && occursin("_d", f), readdir(pardir)))
    t_max = 200.0
    for aspect_ratio in (0.25, 0.5, 1.0, 2.0, 4.0)
        run_script("scripts/pair/pair.jl", :orient => :parallel, :plot => false,
                   :aspect_ratios => [aspect_ratio], :Δx_range => [8.0],
                   :t_max => t_max, :outdir => pardir;
                   done=trajectory_done(trajectory_file(pardir, aspect_ratio, 8.0), t_max))
    end

    # antiparallel trajectories are run per (aspect ratio, separation) pair
    antidir = joinpath(DATA, "antiparallel")
    t_max = 200.0
    for (aspect_ratio, Δx) in ((0.25, 8.0), (1.0, 8.0), (4.0, 8.0), (16.0, 2.0))
        run_script("scripts/pair/pair.jl", :orient => :antiparallel, :plot => false,
                   :aspect_ratios => [aspect_ratio], :Δx_range => [Δx],
                   :t_max => t_max, :outdir => antidir;
                   done=trajectory_done(trajectory_file(antidir, aspect_ratio, Δx), t_max))
    end
end

function suspension()
    patchdir = joinpath(DATA, "patch")
    for aspect_ratio in (0.25, 1.0, 4.0)
        t_end = aspect_ratio == 1.0 ? 500.0 : 180.0
        framedir = joinpath(patchdir, "frames_$(aspect_ratio)")
        run_script("scripts/suspension/suspension.jl", :aspect_ratio => aspect_ratio,
                   :t_end => t_end, :plot => false, :outdir => framedir;
                   done=frames_done(framedir, t_end))

        # far-field model with the monopole closure, from the same initial configuration
        # and up to the final time of the direct run
        ffdir = joinpath(patchdir, "farfield_$(aspect_ratio)")
        run_script("scripts/farfield/farfield_patch.jl", :aspect_ratio => aspect_ratio,
                   :closure => true, :framedir => framedir, :outdir => ffdir;
                   done=frames_done(ffdir, t_end))
    end

    # the expansion snapshot is taken from the disk run
    run_script("figs/expansion_snapshot.jl", :framedir => joinpath(patchdir, "frames_1.0"),
               :outdir => patchdir; done=isfile(joinpath(patchdir, "expansion_snapshot.dat")))
end

function tandem()
    tandemdir = joinpath(DATA, "tandem")
    N, d, m, ε, t_end = 50, 10.0, 5, 0.01, 1000.0
    for aspect_ratio in (0.25, 1.0, 4.0)
        framedir = joinpath(tandemdir, "tandem_$(aspect_ratio)_N$(N)_d$(d)_m$(m)_ε$(ε)")
        run_script("scripts/suspension/tandem.jl", :aspect_ratio => aspect_ratio,
                   :N => N, :d => d, :m => m, :ε => ε, :δθ => 0.0, :t_end => t_end,
                   :save_every => 10.0, :plot => false, :outdir => framedir;
                   done=frames_done(framedir, t_end))
    end

    run_script("scripts/farfield/farfield_tandem_stability.jl", :N => N, :d => d, :m => m,
               :ε => ε, :t_end => t_end, :outdir => tandemdir;
               done=all(isfile, [joinpath(tandemdir, "tandem_$(f)_$(r).dat")
                                 for f in ("stability", "prediction"), r in (0.25, 1.0, 4.0)]))

    run_script("figs/tandem_snapshot.jl", :tandemdir => tandemdir, :outdir => tandemdir;
               done=all(isfile, [joinpath(tandemdir, "tandem_snapshot_$(r)_t$(t).dat")
                                 for r in (0.25, 1.0, 4.0), t in (0, 1000)]))

    # taller field around the unstable chain, for the graphical abstract
    gadir = joinpath(DATA, "graphical_abstract")
    run_script("figs/tandem_snapshot.jl", :aspect_ratios => [4.0], :t_snaps => [1000.0],
               :yh_min => 75.0, :tandemdir => tandemdir, :outdir => gadir;
               done=isfile(joinpath(gadir, "tandem_snapshot_4.0_t1000.dat")))
end

# -- Run -----------------------------------------------------------------------

groups = Dict("single" => single, "pairs" => pairs, "suspension" => suspension, "tandem" => tandem)
selected = isempty(GROUPS) ? ["single", "pairs", "suspension", "tandem"] : GROUPS
for name in selected
    haskey(groups, name) || error("Unknown group $name; choose from $(join(keys(groups), ", "))")
end
println("Writing data to ", DATA, RESUME ? " (resuming)" : "")
for name in selected
    println("\n######## ", name)
    groups[name]()
end
