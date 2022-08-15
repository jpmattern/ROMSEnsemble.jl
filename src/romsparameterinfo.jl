#
# general parameter info
#

abstract type ROMSParameterInfo end

function get_names(rpi::ROMSParameterInfo) end

function get_values(rpi::ROMSParameterInfo, num::Int) end

function get_normalizing_values(rpi::ROMSParameterInfo)
    return nothing
end

#
# latin hypercube sampling for normalized parameters
#


mutable struct LHSROMSParameterInfo  <: ROMSParameterInfo
    parameter_names :: Array{String, 1}
    normalizing_values :: Array{Float64, 1}
    normrange :: Array{Float64}
    function LHSROMSParameterInfo(parameter_names, normalizing_values, normrange)
        @assert(length(parameter_names)==length(normalizing_values), "Arrays containing parameter names and normalizing values must have the same number of entries.")
        @assert(length(normrange)==2, "Input normrange requires 2 entries.")
        rpi = new()
        rpi.parameter_names = parameter_names
        rpi.normalizing_values = normalizing_values
        rpi.normrange = normrange
        return rpi
    end
end

function get_names(rpi::LHSROMSParameterInfo)
    return rpi.parameter_names
end

function get_normalizing_values(rpi::LHSROMSParameterInfo)
    return rpi.normalizing_values
end

function get_values(rpi::LHSROMSParameterInfo, num::Int)
    # latin hypercube sampling
    num_p = length(rpi.parameter_names)
    v = zeros(num, num_p)
    v[1, :] .= 1.0
    for icol in 1:num_p
        v[2:end, icol] = Random.shuffle((Random.rand(num-1)+collect(0:num-2))./((num-1)/(rpi.normrange[2]-rpi.normrange[1])).+rpi.normrange[1])
    end
    # de-normalize
    return v.*transpose(rpi.normalizing_values)
end

function get_values(rpi::LHSROMSParameterInfo, oldvalues::Array{Float64, 2})
    # latin hypercube sampling
    (num, num_p) = size(oldvalues)
    @assert(num_p == length(rpi.parameter_names), "Invalid input size, number of parameters does not match.")
    v = ones(num, num_p)
    for icol in 1:num_p
        v[2:end, icol] = Random.shuffle((Random.rand(num-1)+collect(0:num-2))./((num-1)/(rpi.normrange[2]-rpi.normrange[1])).+rpi.normrange[1])
        v[:, icol] .*= oldvalues[1, icol]
    end
    return v
end

