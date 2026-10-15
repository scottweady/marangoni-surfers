
"""
    script_param(name, default)

Value of a script parameter, overridden from the command line by an argument `name=value`.

The value is parsed as a Julia expression (e.g. `orient=:parallel`,
`aspect_ratios=[0.25, 1.0]`), except when `default` is a string, in which case it is taken
verbatim (e.g. `outdir=/path/to/data`).

# Arguments
- `name`: parameter name, as a `Symbol` or string.
- `default`: value used when `name` is not given on the command line.

# Returns
- The command line value if given, otherwise `default`.
"""
function script_param(name, default)
    for arg in ARGS
        key, val = occursin("=", arg) ? split(arg, "="; limit=2) : (arg, "")
        key == string(name) || continue
        return default isa AbstractString ? String(val) : Core.eval(Main, Meta.parse(val))
    end
    return default
end
