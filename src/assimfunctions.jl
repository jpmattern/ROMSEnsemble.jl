function assimilate_stub(mod_values::Array{Float64, 2}, obs_values::Array{Float64, 1}, rstfiles::Array{String, 1}, obs_info::Dict{String, Any})
    @info """input:
          mod_values: $(size(mod_values))
          obs_values: $(size(obs_values))
            rstfiles: $(size(rstfiles))
            obs_type: $(size(obs_info["obs_type"]))
    """
end
