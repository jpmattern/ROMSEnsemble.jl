"""
    ConfigurationError(message::String, key::String="")

An error or problem related to the ROMS configuration.
"""
struct ConfigurationError <: Exception
    message::String
    key::String
    function ConfigurationError(message::String, key::String="")
        return new(message, key)
    end
end

function Base.showerror(io::IO, e::ConfigurationError)
    if length(e.key) > 0
        print(io, "ConfigurationError: $(e.message)\nProblem related to entry: \"$(e.key)\".")
    else
        print(io, "ConfigurationError: $(e.message)")
    end
end

"""
    ROMSError(message::String, filename::String="")

A ROMS error related to a failure of ROMS to run, or failure to produce desired output.
"""
struct ROMSError <: Exception
    message::String
    filename::String
    function ROMSError(message::String, filename::String="")
        return new(message, filename)
    end
end

function Base.showerror(io::IO, e::ROMSError)
    if length(e.filename) > 0
        print(io, "ROMSError: $(e.message)\nProblem in file: \"$(e.filename)\".")
    else
        print(io, "ROMSError: $(e.message)")
    end
end
