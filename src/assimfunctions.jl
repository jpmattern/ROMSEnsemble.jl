function af_stateestimation_stub(cycle::Int, iter::Int, mod_values::Array{Float64, 2}, obs_values::Array{Float64, 1}, rstfiles::Array{String, 1}, obs_info::Dict{String, Any})
    @info """input:
               cycle: $(cycle)
           iteration: $(iter)
          mod_values: $(size(mod_values)) matrix
          obs_values: $(size(obs_values)) vector
            rstfiles: $(size(rstfiles)) vector of file names
            obs_type: $(size(obs_info["obs_type"]))
    """
end
