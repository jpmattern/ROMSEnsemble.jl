struct ConfigurationError <: Exception
    message::String
    key::String
    function ConfigurationError(message::String, key::String)
        return new(message, key)
    end
    function ConfigurationError(message::String)
        return new(message, "")
    end
end

function Base.showerror(io::IO, e::ConfigurationError)
    if length(e.key) > 0
        print(io, "ConfigurationError: $(e.message)\nProblem related to entry: \"$(e.key)\".")
    else
        print(io, "ConfigurationError: $(e.message)")
    end
end

struct ROMSError <: Exception
    message::String
    filename::String
    function ROMSError(message::String, filename::String)
        return new(message, filename)
    end
    function ROMSError(message::String)
        return new(message, "")
    end
end

function Base.showerror(io::IO, e::ROMSError)
    if length(e.filename) > 0
        print(io, "ROMSError: $(e.message)\nProblem in file: \"$(e.filename)\".")
    else
        print(io, "ROMSError: $(e.message)")
    end
end
