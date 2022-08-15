function get_variable(fname::String, varname::String)
    re = Regex("^ *$(varname) *==? *([^! ]+)")
    open(fname) do file
        for l in eachline(file)
            m = match(re, l)
            if m !== nothing
                return m.captures[1]
            end
        end
    end
end

function set_variable(fname::String, varname::String, newvalue::Int)
    set_variable(fname, varname, "$(newvalue)")
end

function set_variable(fname::String, varname::String, newvalue::Float64)
    set_variable(fname, varname, @sprintf("%fd0", newvalue))
end

function set_variable(fname::String, varname::String, newvalue::String)
    txt = open(fname) do file
        read(file, String)
    end
    # note the use of "m" flag here
    txt = replace(txt, Regex("(^ *$(varname) *==? *)[^! \n]+([^\n]*)", "m") => SubstitutionString("\\g<1>$(newvalue)\\g<2>") )
    open(fname, "w") do file
        write(file, txt)
    end
    nothing
end
