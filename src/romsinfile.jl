struct ROMSInputFile
    filename :: String
    function ROMSInputFile(filename::String)
        if ! isfile(filename)
            error("ROMS input file \"$(filename)\" does not exist.")
        end
        new(filename)
    end
end

function get(rif::ROMSInputFile, varname::String)
    re = Regex("^ *$(varname) *==? *([^! ]+)")
    open(rif.filename) do file
        for l in eachline(file)
            m = match(re, l)
            if m !== nothing
                return m.captures[1]
            end
        end
    end
end

function set!(rif::ROMSInputFile, varname::String, newvalue::Int)
    set!(rif, varname, "$(newvalue)")
end

function set!(rif::ROMSInputFile, varname::String, newvalue::Float64)
    set!(rif, varname, @sprintf("%fd0", newvalue))
end

function set!(rif::ROMSInputFile, varname::String, newvalue::String)
    txt = open(rif.filename) do file
        read(file, String)
    end
    # note the use of "m" flag here
    txt = replace(txt, Regex("(^ *$(varname) *==? *)[^! \n]+([^\n]*)", "m") => SubstitutionString("\\g<1>$(newvalue)\\g<2>") )
    open(rif.filename, "w") do file
        write(file, txt)
    end
    nothing
end

function Base.getindex(rif::ROMSInputFile, varname::String)
    return get(rif, varname)
end

function Base.setindex!(rif::ROMSInputFile, newvalue, varname::String)
    set!(rif, varname, newvalue)
    nothing
end

function Base.show(io::IO, rif::ROMSInputFile)
    print(io, "ROMS input file \"$(rif.filename)\"")
end
