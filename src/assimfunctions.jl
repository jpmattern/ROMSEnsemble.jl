function assimilate_stub(mod_values::Array{Float64, 2}, obs_values::Array{Float64, 1}, rstfiles::Array{String, 1}, obs_info::Dict{String, Any}; i_bg::Int=1, localization_mode::Int=0, num_inner::Int=3, logtransform::Bool=true, minval::Float64=1e-6)
    @info """input:
          mod_values: $(size(mod_values))
          obs_values: $(size(obs_values))
            rstfiles: $(size(rstfiles))
            obs_type: $(size(obs_info["obs_type"]))
         obs_type==9: $(sum(obs_info["obs_type"] .== 9))
    """
end
