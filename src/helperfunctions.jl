"""
Obtain time information from the netCDF file `ncfile`.
"""
function get_time(ncfile::AbstractString; timevar::String="ocean_time")
    NCDatasets.Dataset(ncfile) do nc
        return nc[timevar][:]
    end
end

"""
Set the output file name variables in the ROMS input files specified in `config`.
"""
function set_output_names!(config::Dict{String, Any}, suffix::String)
    # set local filenames
    rif_ocean = ROMSInputFile(config["ocean_in"])
    for ft in ("gst", "rst", "his", "tlm", "tlf", "adj", "avg", "dia", "sta", "flt")
        k = uppercase(ft)*"NAME"
        #v = "$(config["file_prefix"])_$(ft)_$(suffix).nc"
        if ft == "rst"
            v = "$(config["file_prefix"])_$(ft)_$(suffix).nc"
        else
            # only interested in mod and rst files, leave out suffix, so that files are automatically overwritten
            v = "$(config["file_prefix"])_$(ft).nc"
        end
        #@info "for template: setting $k=$v"
        rif_ocean[k] = v
    end
    if length(config["s4dvar_in"]) > 0
        ft = "mod"
        k = uppercase(ft) * "name"
        v = "$(config["file_prefix"])_$(ft)_$(suffix).nc"
        #@info "for template: setting $k=$v"
        rif_4dvar = ROMSInputFile(config["s4dvar_in"])
        rif_4dvar[k] = v
    end
    nothing
end

"""
Get parameter values from the ROMS output netCDF file `ncfile`.
"""
function get_parameter_values(ncfile::String, varnames::Array{String}; checkagainst=nothing, thresh::Float64=1e-6) :: Array{Float64}
    NCDatasets.Dataset(ncfile) do nc
        values = [ncread[v][1] for v in varnames]
    end
    if checkagainst !== nothing
        if any(abs.(values-checkagainst).>thresh)
            println(@sprintf("%20s %13s %13s", "", "current value", "expected"))
            for (ip, p) in enumerate(varnames)
                if abs(values[ip]-checkagainst[ip]) <= thresh
                    s = "✓"
                else
                    s = "x"
                end
                println(@sprintf("%20s %13.9f %13.9f %s", p, values[ip], checkagainst[ip], s))
            end
            error("Failed parameter value check for file \"$ncfile\".")
        end
    end
    return values
end

"""
Turn the path `p` into an absolute one.
"""
function make_abs(p::String) :: String
    if p[1] ≠ '/'
        return joinpath(pwd(), p)
    else
        return p
    end
end


function _replace_patterns(d::Dict, patterns::Tuple{Tuple{String, String}})
    for (k,v) in d
        if v isa AbstractString
            for pattern in patterns
                # replace both $HOME and ${HOME}
                d[k] = replace(d[k], "\$" * pattern[1] => pattern[2])
                d[k] = replace(d[k], "\${" * pattern[1] * "}" => pattern[2])
            end
        elseif v isa Dict
            _replace_patterns(d[k], patterns)
        end
    end
end


"""
Parse a configuration file to obtain the configuration in `Dict{String, Any}` format.
"""
function parse_config(input::AbstractString) :: Dict{String, Any}
    config = YAML.load(open(input); dicttype=Dict{String, Any}) :: Dict{String, Any}
    # search and replace $HOME or ${HOME}
    patterns = (("HOME", homedir()),)
    _replace_patterns(config, patterns)

    # initial_conditions
    if haskey(config, "initial_conditions") && config["initial_conditions"] isa AbstractString
        config["initial_conditions"] = glob(basename(config["initial_conditions"]), dirname(config["initial_conditions"]))
    end
    # assimilation_function
    if haskey(config, "assimilation_function") && config["assimilation_function"] isa AbstractString
        config["assimilation_function"] = eval(Meta.parse(config["assimilation_function"])) :: Function
    end

    return config
end
