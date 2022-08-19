function get_time(ncfile::AbstractString; timevar::String="ocean_time")
    NCDatasets.Dataset(ncfile) do nc
        return nc[timevar][:]
    end
end

function set_output_names!(config::Dict, suffix::String)
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
        k = uppercase(ft)*"name"
        v = "$(config["file_prefix"])_$(ft)_$(suffix).nc"
        #@info "for template: setting $k=$v"
        rif_4dvar = ROMSInputFile(config["s4dvar_in"])
        rif_4dvar[k] = v
    end
    nothing
end

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

function make_abs(s::String) :: String
    if s[1] != '/'
        return joinpath(pwd(), s)
    else
        return s
    end
end

function parse_config(input::AbstractString) :: Dict{String, Any}
    config = YAML.load(open(input); dicttype=Dict{String, Any}) :: Dict{String, Any}
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
