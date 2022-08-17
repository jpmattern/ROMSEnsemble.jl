mutable struct ROMSFileManager
    rundirs :: Array{String, 1}
    files :: Dict{String, String}
    inkeys :: Dict{String, Tuple{String, String}}
    shortnames :: Dict{String, String}
    localnames :: Dict{String, String}
    logfilename :: String
    function ROMSFileManager(rundirs,
                     oceanin::String;
                     bioin::String="",
                     s4dvarin::String="",
                     obs::String="",
                     logfilename::String="roms.log",
                     useshortnames::Bool=true
                     )
        rfm = new()
        rfm.inkeys = Dict("bioin"=>("oceanin", "BPARNAM"), "s4dvarin"=>("oceanin", "APARNAM"), "obs"=>("s4dvarin", "OBSname"))
        rfm.shortnames = Dict("oceanin"=>"ocean.in", "bioin"=>"bio.in", "s4dvarin"=>"s4dvar.in", "obs"=>"obs.nc")
        rfm.localnames = Dict("oceanin"=>"ocean.in", "bioin"=>"bio.in", "s4dvarin"=>"s4dvar.in", "obs"=>"obs.nc")
        rfm.logfilename = logfilename

        if rundirs isa String
            rfm.rundirs = glob(basename(rundirs), dirname(rundirs))
            for d in rfm.rundirs
                @assert(isdir(d), "\"$d\" must be a directory.")
            end
        elseif rundirs isa Array{String, 1}
            for d in rundirs
                bd = dirname(d)
                @assert(isdir(bd), "Base directory \"$bd\" must be a directory.")
            end
            rfm.rundirs = rundirs
        else
            error("Invalid input for rundirs.")
        end

        @assert(isfile(oceanin), "File \"$oceanin\" does not exist.")
        rfm.files = Dict("oceanin" => oceanin)

        if length(bioin) > 0
            @assert(isfile(bioin), "File \"$bioin\" does not exist.")
            rfm.files["bioin"] = bioin
        end

        if length(s4dvarin) > 0
            @assert(isfile(s4dvarin), "File \"$s4dvarin\" does not exist.")
            rfm.files["s4dvarin"] = s4dvarin
        end

        if length(obs) > 0
            @assert(isfile(obs), "File \"$obs\" does not exist.")
            @assert("s4dvarin" in keys(rfm.files), "Obervation file specification requires s4dvarin specification.")
            rfm.files["obs"] = obs
        end

        if ! useshortnames
            for (k, f) in rfm.files
                rfm.localnames[k] = basename(f)
            end
        end

        rfm
    end
end

function Base.show(io::IO, rfm::ROMSFileManager)
    println(io, "ROMSFileManager")
    println(io, "directories")
    for (i, d) in enumerate(rfm.rundirs)
        println(io, "$(i): \"$(d)\"")
    end
end

function create_directories(rfm::ROMSFileManager)
    for d in rfm.rundirs
        if ! isdir(d)
            mkdir(d)
        end
    end
    nothing
end

function create_files(rfm::ROMSFileManager; verbose::Bool=false)
    create_directories(rfm)

    currentfiles = Dict()

    for d in rfm.rundirs
        if verbose
            println("directory \"$d\"")
        else
            @debug("directory \"$d\"")
        end

        # copy files
        for (k, f) in rfm.files
            currentfiles[k] = joinpath(d, rfm.localnames[k])
            if ! isfile(f)
                @error("Cannot find \"$(f)\".")
            end
            cp(f, currentfiles[k], force=true)
        end

        # set ROMS variables
        for (linktarget, linkinfo) in rfm.inkeys
            if haskey(currentfiles, linktarget)
                if verbose
                    println("   setting $(basename(currentfiles[linkinfo[1]])): $(linkinfo[2]) = $(currentfiles[linktarget])")
                else
                    @debug("   setting $(basename(currentfiles[linkinfo[1]])): $(linkinfo[2]) = $(currentfiles[linktarget])")
                end
                rif = ROMSInputFile(currentfiles[linkinfo[1]])
                rif[linkinfo[2]] = currentfiles[linktarget]
            end
        end
    end
end

function set!(rfm::ROMSFileManager, dir::String, infilekey::String, key::String, value::String)
    # set_variable(joinpath(dir, rfm.localnames[infilekey]), key, value)
    rif = ROMSInputFile(joinpath(dir, rfm.localnames[infilekey]))
    rif[key] = value
    nothing
end

function set!(rfm::ROMSFileManager, infilekey::String, key::String, value::String)
    for d in rfm.rundirs
        set!(rfm, d, infilekey, key, value)
    end
    nothing
end

function set!(rfm::ROMSFileManager, infilekey::String, key::String, vals::Array{String, 1})
    if length(vals) == 1
        for (i, d) in enumerate(rfm.rundirs)
            set!(rfm, d, infilekey, key, vals[1])
        end
        return
    end
    @assert(length(rfm.rundirs) == length(vals), "Length of values ($(length(vals))) must match number of runs ($(length(rfm.rundirs))).")
    for (i, d) in enumerate(rfm.rundirs)
        set!(rfm, d, infilekey, key, vals[i])
    end
    nothing
end

function move_output(rfm::ROMSFileManager, todir::String; filenamemask::String="*.nc", verbose::Bool=false)
    @assert(isdir(todir), "Directory \"$(todir)\" must exist.")
    for (i, d) in enumerate(rfm.rundirs)
        suffix = @sprintf("%03d", i)
        move_if_exists(joinpath(d, rfm.logfilename), suffix, todir, verbose=verbose)
        fnames = glob(filenamemask, d)
        for fname in fnames
            move_if_exists(fname, suffix, todir, verbose=verbose)
        end
    end
    nothing
end

function move_if_exists(fname::String, addsuffix::String, todir::String; verbose::Bool=false)
    if isfile(fname)
        parts = splitext(basename(fname))
        fname_new = joinpath(todir, "$(parts[1])_$(addsuffix)$(parts[2])")
        if verbose
            @info("moving \"$(fname)\" => \"$(fname_new)\"")
        end
        mv(fname, fname_new, force=true)
    end
    nothing
end
