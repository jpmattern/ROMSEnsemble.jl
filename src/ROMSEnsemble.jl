module ROMSEnsemble

using Statistics
using LinearAlgebra

import ArgParse
import ConfParser
import YAML
import NCDatasets
import Dates

import Printf: @sprintf
import Glob: glob

export run, parse_config

include("assimfunctions.jl")
include("helperfunctions.jl")
include("romsfilemanagers.jl")
include("romsinputfile.jl")
include("romsparameterinfo.jl")
include("romsstarters.jl")
include("errors.jl")

"""
    run(configfile::AbstractString; kwargs...)

Start a ROMS ensemble simulation, based on the configuration specified in the file `configfile`.
"""
function run(configfile::AbstractString; kwargs...)
    config = parse_config(configfile)
    return run(config; kwargs...)
end

"""
    run(config::Dict{String, Any}; allow_skip_da::Bool=true, display_function::Union{Nothing, Function}=nothing)

Start a ROMS ensemble simulation, based on the configuration specified in `config`.
"""
function run(config::Dict{String, Any}; allow_skip_da::Bool=true, display_function::Union{Nothing, Function}=nothing)

    #
    # get information / perform checks
    #

    cdate = nothing
    minval = 1e-3

    # get current date from initial file
    if config["initial_conditions"] isa AbstractString
        cdates = get_time(config["initial_conditions"])
        if length(cdates) ≠ 1
            throw(ConfigurationError("Initial conditions file must have one time step.", "initial_conditions"))
        end
        cdate = cdates[1]
    elseif config["initial_conditions"] isa Array{String, 1}
        if length(config["initial_conditions"]) > 1 && length(config["initial_conditions"]) ≠ config["n_ens"]
            throw(ConfigurationError("Number of initial conditions files must be 1 or the number of ensemble members ($(config["n_ens"])); received $(length(config["initial_conditions"])).", "initial_conditions"))
        end
        for (i, f) in enumerate(config["initial_conditions"])
            cdates = get_time(f)
            if length(cdates) ≠ 1
                throw(ConfigurationError("Initial conditions file \"$(f)\" must have one time step.", "initial_conditions"))
            end
            if i == 1
                cdate = cdates[1]
            elseif cdate ≠ cdates[1]
                throw(ConfigurationError("Initial conditions files must all have the same time (problem with \"$(f)\").", "initial_conditions"))
            end
        end
    else
        throw(ConfigurationError("Invalid format for \"initial_conditions\"."))
    end

    if ! haskey(config, "assimilation_function")
        throw(ConfigurationError("No assimilation function specified.", "assimilation_function"))
    end

    if ! haskey(config, "costtype")
        config["costtype"] = "normcost"
    end

    if ! haskey(config, "estimationtype")
        config["estimationtype"] = "bioparameter"
    end

    if ! haskey(config, "file_prefix")
        config["file_prefix"] = "romsensemble"
    end

    if ! haskey(config, "spinup_days")
        config["spinup_days"] = 0
    end

    if ! haskey(config, "assimilation_function_usefileinput")
        config["assimilation_function_usefileinput"] = false
    end

    if ! haskey(config, "cleanup_on_success")
        config["cleanup_on_success"] = false
    end

    if config["estimationtype"] == "bioparameter"
        if ! haskey(config, "negative_parameter_response")
            config["negative_parameter_response"] = "adjust"
        else
            if config["negative_parameter_response"] ∉ ("warn", "error", "adjust")
                throw(ConfigurationError("Value of negative_parameter_response must be one of \"warn\", \"error\", \"adjust\"."))
            end
        end
        if config["negative_parameter_response"] == "adjust"
            if ! haskey(config, "negative_parameter_minval")
                config["negative_parameter_minval"] = 0.0
            end
        end
    end

    # make some paths abolute
    for key in ["rundir", "storagedir", "obsfile", "modfile", "executable", "ocean_in", "bio_in", "s4dvar_in", "initial_conditions"]
        if haskey(config, key)
            if config[key] isa AbstractString
                #@assert(isfile(config[key]), "File \"$(config[key])\" does not exist.")
                config[key] = make_abs(config[key])
            else
                config[key] = [make_abs(f) for f in config[key]]
            end
        end
    end


    #
    # create rundir and templatedir
    #

    if ! isdir(config["rundir"])
        mkdir(config["rundir"])
    end
    if ! isdir(config["storagedir"])
        mkdir(config["storagedir"])
    end
    if ! haskey(config, "templatedir")
        config["templatedir"] = joinpath(config["rundir"], ".templates")
    end
    if ! isdir(config["templatedir"])
        mkdir(config["templatedir"])
    end

    #
    # extract data
    #

    if ! haskey(config, "obsfile_variables")
        config["obsfile_variables"] = ["obs_time", "obs_value", "obs_error", "obs_type", "obs_provenance"]
    end

    obs_info = Dict{String, Any}()
    NCDatasets.Dataset(config["obsfile"]) do nc
        for v in config["obsfile_variables"]
            obs_info[v] = nc[v][:]
        end
        for v in ("obs_value", "obs_time")
            if ! haskey(obs_info, v)
                obs_info[v] = nc[v][:]
            end
        end
    end

    # TODO-low generalize this
    #obsindex_general = obs_info["obs_type"] .== 9
    obsindex_general = trues(size(obs_info["obs_type"]))

    #
    # obtain obs
    #

    if haskey(config, "use_modforobs") && config["use_modforobs"]
        if ! haskey(config, "modfile")
            throw(ConfigurationError("A mod-file must be specified when \"use_modforobs\" is active.", "modfile"))
        end
        @info "using NLmodel_value from mod-file as observations (source: \"$(config["modfile"])\")."
        NCDatasets.Dataset(config["modfile"]) do nc
            obs_values = nc["NLmodel_value"][:]
        end
        if length(obs_info["obs_time"]) ≠ length(obs_values)
            throw(ConfigurationError("The mod-file is not compatible with the specified observation file."))
        end
        #=
        # TODO-low implement obs error here
        obs_values = mod_values_goal
        if obsnoise > 0.0
            @warn "adding $(obsnoise*100)% error to observations"
            scenario_desc *= ", $(obsnoise*100)% observation noise"
            obs_values .+= obsnoise .* obs_values .* Random.randn(length(obs_values))
        end
        =#
    else
        @info "using observations from observation file \"$(config["obsfile"])\"."
        obs_values = obs_info["obs_value"]
    end
    #obs_values = max.(minval, obs_values)

    #
    # create and manipulate template files
    #

    # copy files to templatedir
    for k in ("ocean_in", "bio_in", "s4dvar_in", "obsfile")
        if (k in keys(config)) & (length(config[k]) > 0)
            f = joinpath(config["templatedir"], basename(config[k]))
            cp(config[k], f, force=true)
            config[k] = f
        end
    end

    rif_ocean = ROMSInputFile(config["ocean_in"])
    for (k, v) in config["paramchanges_ocean"]
        @debug "for template: setting $k=$v"
        rif_ocean[k] = v
    end

    dt = parse(Float64, replace(rif_ocean["DT"], "d"=>"e"))
    refdate = Dates.DateTime(replace(rif_ocean["TIME_REF"], r"\.0*d0*" => ""), "yyyymmdd")
    ntasks = parse(Int, rif_ocean["NtileI"])
    ntasks *= parse(Int, rif_ocean["NtileJ"])
    @debug """dt=$dt
    refdate=$refdate
    ntasks=$ntasks"""

    #
    # compute stop dates
    #

    if haskey(config, "cycle_length")
        if haskey(config, "stopdates")
            @warn("Using \"cycle_length\" to determine stop dates, ignoring value of \"stopdates\" in configuration.")
        end
        if ! haskey(config, "num_cycles")
            throw(ConfigurationError("Variable \"num_cycles\" expected in configuration.", "num_cycles"))
        end

        stopdates = [cdate+Dates.Day(config["spinup_days"] + i*config["cycle_length"]) for i in 1:config["num_cycles"]]

        @debug("stop dates: $(stopdates)")
    else
        stopdates = [Dates.DateTime(d) for d in config["stopdates"]]

        keep = trues(length(stopdates))
        laststop = cdate
        for (i, stopdate) in enumerate(stopdates)
            if stopdate <= cdate
                keep[i] = false
                println("stop date $i ($(stopdate)) is in the past and will be eliminated")
                continue
            elseif config["spinup_days"] > 0 && stopdate <= cdate + Dates.Day(config["spinup_days"])
                keep[i] = false
                println("stop date $i ($(stopdate)) is in the spinup period and will be eliminated")
                continue
            end
            nsteps_float = (stopdate-laststop).value/(dt*1000)
            nsteps = Int(round(nsteps_float))
            # TODO-low deal with this
            #if nsteps_float - round(nsteps_float) > 1e-10
            #
            #end
            if nsteps < config["timestep_threshold"]
                keep[i] = false
                println("stop date $i ($(stopdate)) is too close to the previous one ($(laststop)) and will be eliminated")
                continue
            end
            laststop = stopdate
        end
        stopdates = stopdates[keep]
        @info "stop dates after elimination: $(stopdates)"
    end
    num_cycles = length(stopdates)
    num_ens = config["n_ens"]

    #
    # determine parameter/state estimation
    #

    if config["estimationtype"] == "bioparameter"
        rpi = config["ROMSParameterInfo"]
        parameter_names = get_names(rpi)
        parameter_values = get_values(rpi, num_ens)
    elseif config["estimationtype"] == "state"
        # TODO preparation needed?
    elseif config["estimationtype"] != "spinup-only"
        throw(ConfigurationError("Invalid estimation type \"$(config["estimationtype"])\".", "estimationtype"))
    end

    refillensemble = false
    startfrom1 = false
    iterconstini = false
    if config["cycle_setup"] isa Array
        cycle_setup = config["cycle_setup"]
        startfrom1 = config["startfrom1"]::Bool
        if haskey(config, "noassimlastiter")
            noassimlastiter = config["noassimlastiter"]::Bool
        else
            noassimlastiter = true
        end
        iterconstini = config["iterconstini"]::Bool
        if haskey(config, "refillensemble")
            refillensemble = config["refillensemble"]::Bool
            @debug "refillensemble: $(config["refillensemble"]::Bool)"
        end
    elseif config["cycle_setup"] == "4DEnOI"
        cycle_setup = [num_ens, 1]
        if config["estimationtype"] == "bioparameter"
            startfrom1 = true
            iterconstini = true
        else
            refillensemble = true
        end
        noassimlastiter = true
    elseif config["cycle_setup"] == "EnKF"
        cycle_setup = [num_ens, num_ens] # NOTE changed 2020-05-04
        if config["estimationtype"] == "bioparameter"
            iterconstini = true
        end
        noassimlastiter = true
    else
        throw(ConfigurationError("Invalid cycle setup \"$(config["cycle_setup"])\".", "cycle_setup"))
    end
    num_iter = length(cycle_setup)

    #
    # create directory structure and set up ROMSFileManager and ROMSStarter
    #

    rundirs = [joinpath(config["rundir"], @sprintf("%03d", i)) for i in 1:num_ens]

    if haskey(config, "bio_in")
        rfm = ROMSFileManager(rundirs, config["ocean_in"], bioin=config["bio_in"], s4dvarin=config["s4dvar_in"], obs=config["obsfile"])
    else
        rfm = ROMSFileManager(rundirs, config["ocean_in"], s4dvarin=config["s4dvar_in"], obs=config["obsfile"])
    end

    if ! haskey(config, "ROMSStarter")
        throw(ConfigurationError("The ROMSStarter variable is missing from the configuration.", "ROMSStarter"))
    end

    if config["ROMSStarter"] isa ROMSStarter
        rs = config["ROMSStarter"]::ROMSStarter
    elseif config["ROMSStarter"] == "sbatch"
        if ! haskey(config, "batchscript")
            throw(ConfigurationError("Using the sbatch ROMS starter requires specification of \"batchscript\" variable.", "batchscript"))
        end
        if ! isfile(config["batchscript"])
            throw(ConfigurationError("The batch script file \"$(config["batchscript"])\" does not exist.", "batchscript"))
        end
        rs = SbatchROMSStarter(num_ens, config["executable"], ntasks, config["batchscript"], sleeptime=1.0, jobnameprefix="ROMSEnsemble.jl")
    elseif config["ROMSStarter"] == "srun"
        rs = SrunROMSStarter(num_ens, config["executable"], ntasks; sleeptime=1.0)
    end

    #
    # perform spinup
    #

    if config["spinup_days"] > 0
        suffix = "000-spinup"
        num_ens_curr = num_ens
        if config["initial_conditions"] isa  AbstractString || length(config["initial_conditions"]) == 1
            num_ens_curr = 1
        else
            if length(config["initial_conditions"]) ≠ num_ens
                throw(ConfigurationError("Number of initial conditions files must be 1 or the number of ensemble members ($num_ens)."))
            end
        end

        if isfile(joinpath(config["storagedir"], "$(config["file_prefix"])_rst_$(suffix)_001.nc"))
            println("Found output for spinup -- not starting ROMS.")
        else

            #
            # set parameters
            #

            nsteps = Int(config["spinup_days"]*86400/dt)
            dstart = (cdate-refdate).value/(86400*1000)

            rif_ocean["DSTART"] = dstart
            rif_ocean["NTIMES"] = nsteps
            rif_ocean["NRST"] = nsteps

            set_output_names!(config, suffix)
            rfm.logfilename = "roms_$(suffix).log"
            # move from template to individual files
            create_files(rfm)
            set!(rfm, "oceanin", "ININAME", config["initial_conditions"])

            if config["estimationtype"] == "bioparameter"
                for (i, name) in enumerate(parameter_names)
                    tmp = [@sprintf("%f", x) for x in parameter_values[:, i]]
                    #@info "$(i), $(name): $(tmp)"
                    set!(rfm, "bioin", name, tmp)
                end
            end

            #
            # start ROMS
            #

            start_jobs(rs, rfm, numjobs=num_ens_curr)
            move_output(rfm, config["storagedir"], filenamemask="*_rst_*.nc")

        end

        #
        # check rst files
        #

        rstfiles = glob("$(config["file_prefix"])_rst_$(suffix)_*.nc", config["storagedir"])
        if length(rstfiles) ≠ num_ens_curr
            throw(ROMSError("Expected $(num_ens_curr) restart files but found $(length(rstfiles))."))
        end

        cdate += Dates.Day(config["spinup_days"])
        for (irst, rstfile) in enumerate(rstfiles)
            cdates = get_time(rstfile)
            @info cdates
            if length(cdates) ≠ 1
                throw(ROMSError("Restart file must have one time step, but found $(length(cdates)).", rstfile))
            end
            if cdates[1] ≠ cdate
                throw(ROMSError("Date in restart file ($(cdates[1])) was expected to be $(cdate).", rstfile))
            end
            if config["estimationtype"] == "bioparameter"
                # perform check to see that parameter values are correct
                get_parameter_values(rstfile, varnames=parameter_names, checkagainst=parameter_values)
            end
        end
    else
        # no spinup, convert from initial_conditions to rstfiles
        if config["initial_conditions"] isa AbstractString
            rstfiles = [config["initial_conditions"] for i in 1:num_ens]
        else
            rstfiles = config["initial_conditions"]
        end
    end
    #error(" --- THE END FOR NOW --------------------------------------------------------------------------------")
    if config["estimationtype"] == "spinup-only"
        @info("Spinup-only mode -- quitting.")
        return nothing
    end

    #
    # start loop
    #

    obs_info_sub = Dict{String, Any}()
    suffix_last = nothing
    stats_old = Dict{Int, Dict{String, Float64}}()
    stopdate = stopdates[1]

    for i_cycle in 1:num_cycles
        if i_cycle > 1
            cdate = stopdate
        end
        stopdate = stopdates[i_cycle]
        nsteps = Int(round((stopdate-cdate).value/(dt*1000)))
        dstart = (cdate-refdate).value/(86400*1000)
        #dstart_next = (stopdate-refdate).value/(86400*1000)

        @info("""

        cycle $(i_cycle)/$(num_cycles)
            start: $(cdate)
              end: $(stopdate)
        """)

        #obs_time_min = dstart
        #obs_time_max = dstart_next
        obs_time_min = cdate
        obs_time_max = stopdate

        @debug("""
        dstart=$(dstart)
        nsteps=$(nsteps)
        """)

        #
        # set up new cycle
        #

        rif_ocean["DSTART"] = dstart
        rif_ocean["NTIMES"] = nsteps
        rif_ocean["NRST"] = nsteps

        # keep copy of rstfiles to be used a initial conditions throughout
        # all inner loops if iterconstini is active (for parameter estimation)
        #inifiles = [joinpath(rundirs[i], "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d_%03d", i_cycle, 1, i)).nc") for i in 1:num_ens]
        inifiles = [joinpath(rundirs[i], "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d", i_cycle, 0)).nc") for i in 1:num_ens]
        if length(rstfiles) ≠ num_ens
            throw(ROMSError("Expected to find $num_ens restart files but found $(length(rstfiles))."))
        end
        create_directories(rfm)
        for (irst, rstfile) in enumerate(rstfiles)
            @info "copying $rstfile => $(inifiles[irst])"
            cp(rstfile, inifiles[irst], force=true, follow_symlinks=true)
            if startfrom1
                break
            end
        end

        obsindex = copy(obsindex_general)
        obsindex .&= obs_info["obs_time"] .>= obs_time_min
        obsindex .&= obs_info["obs_time"] .<  obs_time_max
        mod_values = fill(NaN, (num_ens, sum(obsindex)))

        for i_iter in 1:num_iter
            @info("""

            cycle $i_cycle iteration/outer loop $i_iter""")
            suffix = @sprintf("%03d_%03d", i_cycle, i_iter-1) # note that first run starts at 000

            num_ens_curr = cycle_setup[i_iter]
            rfm.rundirs = rundirs[1:num_ens_curr]

            skip_roms = isfile(joinpath(config["storagedir"], "$(config["file_prefix"])_mod_$(suffix)_001.nc"))
            if skip_roms
                # make sure that all mod-files are present
                modfiles = glob("$(config["file_prefix"])_mod_$(suffix)_*.nc", config["storagedir"])
                skip_roms = length(modfiles) == num_ens_curr
            end

            if skip_roms
                println("Found output for this iteration -- not starting ROMS.")
            else
                rfm.logfilename = "roms_$(suffix).log"

                set_output_names!(config, suffix)
                create_files(rfm)
                #=
                @info begin
                    s = ""
                    for (i, f) in enumerate(inifiles[1:num_ens_curr])
                        s *= @sprintf("%2d", i) * ") $(f) [$(isfile(f))]\n"
                    end
                    s
                end
                =#


                # move from template to individual files
                if startfrom1
                    set!(rfm, "oceanin", "ININAME", inifiles[1])
                else
                    for (i, f) in enumerate(inifiles[1:num_ens_curr])
                        if ! isfile(f)
                            error("Initial file \"$(f)\" does not exist.")
                        end
                    end
                    set!(rfm, "oceanin", "ININAME", inifiles[1:num_ens_curr])
                end
                if config["estimationtype"] == "bioparameter"

                    #
                    # set parameters
                    #

                    for (i, name) in enumerate(parameter_names)
                        tmp = [@sprintf("%f", x) for x in parameter_values[1:num_ens_curr, i]]
                        #@info "$(i), $(name): $(tmp)"
                        set!(rfm, "bioin", name, tmp)
                    end
                end

                #
                # start ROMS
                #

                start_jobs(rs, rfm)
                @info "moving model output to \"$(config["storagedir"])\""
                move_output(rfm, config["storagedir"], filenamemask="*_[mri]*_*.nc", verbose=false)
            end

            #
            # extract rst and mod values
            #

            if refillensemble
                #inifiles_save = copy(inifiles)
                rstfiles_save = copy(rstfiles)
            end
            modfiles = glob("$(config["file_prefix"])_mod_$(suffix)_*.nc", config["storagedir"])
            rstfiles = glob("$(config["file_prefix"])_rst_$(suffix)_*.nc", config["storagedir"])
            if length(modfiles) ≠ num_ens_curr
                throw(ROMSError("Expected $num_ens_curr mod files but found $(length(modfiles))."))
            end
            if length(rstfiles) ≠ num_ens_curr
                throw(ROMSError("Expected $num_ens_curr restart files but found $(length(rstfiles))."))
            end
            @assert(issorted(modfiles), "The mod files are not in correct order.")
            @assert(issorted(rstfiles), "The restart files are not in correct order.")
            for (irst, rstfile) in enumerate(rstfiles)
                cdates = get_time(rstfile)
                if length(cdates) ≠ 1
                    throw(ROMSError("Restart file must have one time step, but found $(length(cdates)).", rstfile))
                end
                if cdates[1] ≠ stopdate
                    throw(ROMSError("Date in restart file ($(cdates[1])) was expected to be $(stopdate).", rstfile))
                end
                if config["estimationtype"] == "bioparameter"
                    # perform check to see that parameter values are correct
                    get_parameter_values(rstfile, varnames=parameter_names, checkagainst=parameter_values[irst, :])
                #elseif config["estimationtype"] == "state"
                #    inifile = joinpath(dirname(rstfile), replace(basename(rstfile), "$(config["file_prefix"])_rst_"=>"$(config["file_prefix"])_ini_"))
                #    @assert(inifile != rstfile, "Failed to generate initial file name (\"$inifile\")")
                #    @info "copying $rstfile => $inifile"
                #    cp(rstfile, inifile, force=true)
                end
            end
            #if config["estimationtype"] == "state"
            #    # switch to ini files for state estimation
            #    rstfiles = glob("$(config["file_prefix"])_ini_$(suffix)_*.nc", config["storagedir"])
            #end

            obs_scale = nothing
            obsindex_sub = nothing
            @info "extracting $(sum(obsindex)) mod file values"
            for (imod, modfile) in enumerate(modfiles)
                @debug("opening \"$(modfile)\"")
                NCDatasets.Dataset(modfile) do nc
                    ## apply minval
                    #mod_values[imod, :] = max.(minval, nc["NLmodel_value"][obsindex])
                    if length(obsindex) != length(nc["NLmodel_value"])
                        error("Mod file \"$(modfile)\" is no compatible with the observation file \"$(config["obsfile"])\".")
                    end
                    mod_values[imod, :] = nc["NLmodel_value"][obsindex]
                    if imod == 1
                        obs_scale = nc["obs_scale"][obsindex]
                        # obs_scale == 0 values were already removed from obsindex in step above but may pop up due to bioObsThresh
                        #@assert(all(obs_scale.>0), "Found obs_scale = 0.")
                        if any(obs_scale .=== missing)
                            @info "Found $(sum(obs_scale .=== missing)) values with obs_scale = 0 which will be excluded."
                        end
                        if length(unique(obs_scale)) > 2
                            throw(ROMSError("Expected 1.0 and  missing in obs_scale but found $(join(unique(obs_scale), ", ", " and ")).", modfile))
                        end
                        obsindex_sub = obs_scale .!== missing
                    else
                        tmp = nc["obs_scale"][obsindex]
                        if !all(obs_scale .=== tmp)
                            @warn """Found differing number of values obs_scale = 0
                            in reference (\"$(basename(modfiles[1]))\"): $(sum(obs_scale .=== missing))"
                            in file \"$(basename(modfile))\":        $(sum(tmp .=== missing))"""
                            obsindex_sub .&= tmp .!== missing
                        end
                    end
                end
            end

            for k in keys(obs_info)
                obs_info_sub[k] = obs_info[k][obsindex][obsindex_sub]
            end
            obs_values_sub = obs_values[obsindex][obsindex_sub]

            #
            # compute cost functions
            #

            #=
            @warn "writing J[$(i_j)]"
            js_select["cycle"][i_j] = i_cycle
            js_select["iter"][i_j] = i_iter

            #js_select["paramrmse1"][i_j] = sqrt(mean((parameter_values_curr[1, :].-(parameter_values_goal./parameter_values_norm)).^2))
            js_select["meanrmse"][i_j] = mean(sqrt.(mean((mod_values[:, obsindex_sub].-transpose(obs_values[obsindex][obsindex_sub])).^2, dims=2)))
            js_select["rmsemean"][i_j] = sqrt(mean((dropdims(mean(mod_values[:, obsindex_sub], dims=1), dims=1).-obs_values[obsindex][obsindex_sub]).^2))
            js_select["rmse1"][i_j] =    sqrt(mean((mod_values[1, obsindex_sub].-obs_values[obsindex][obsindex_sub]).^2))
            js_select["logrmse1"][i_j] = sqrt(mean((log.(mod_values[1, obsindex_sub]./obs_values[obsindex][obsindex_sub])).^2))
            i_j += 1
            =#

            stats = Dict{Int, Dict{String, Float64}}()
            obs_type_unique = sort(unique(obs_info_sub["obs_type"]))
            for ot in obs_type_unique
                cind = obs_info_sub["obs_type"] .== ot
                stats[ot] = Dict{String, Float64}("numobs" => sum(cind),
                                                 "rmse1" => sqrt(mean((mod_values[1, obsindex_sub][cind].-obs_values_sub[cind]).^2)),
                                                 "bias1" => mean(mod_values[1, obsindex_sub][cind])-mean(obs_values_sub[cind]),
                                                )
                if num_ens_curr > 1
                    stats[ot]["meanrmse"] = mean(sqrt.(mean((mod_values[:, obsindex_sub][:, cind].-transpose(obs_values_sub[cind])).^2, dims=2)))
                    stats[ot]["meanbias"] = mean(mean(mod_values[:, obsindex_sub][:, cind], dims=2).-mean(obs_values_sub[cind]))
                    stats[ot]["rmsemean"] = sqrt(mean((dropdims(mean(mod_values[:, obsindex_sub][:, cind], dims=1), dims=1).-obs_values_sub[cind]).^2))
                    stats[ot]["biasmean"] = mean(mod_values[:, obsindex_sub][:, cind])-mean(obs_values_sub[cind])
                end
            end

            @info begin
                s = "summary:\n"
                s *= "number of observations: $(length(obs_info_sub["obs_type"]))\n"
                rows = ["                 ", "number of obs.: ", "RMSE (index 1): ", "bias (index 1): "]
                for ot in obs_type_unique
                    cind = obs_info_sub["obs_type"] .== ot
                    rows[1] *= @sprintf("obs_type%3d  ", ot)
                    rows[2] *= @sprintf("%12d ",   stats[ot]["numobs"])
                    rows[3] *= @sprintf("%12.6f ", stats[ot]["rmse1"])
                    rows[4] *= @sprintf("%12.6f ", stats[ot]["bias1"])
                end
                irow = 5
                if num_ens_curr > 1
                    append!(rows, ["ens. mean RMSE: ", "ens. mean bias: ", "RMSE ens. mean: ", "bias ens. mean: "])
                    for ot in obs_type_unique
                        cind = obs_info_sub["obs_type"] .== ot
                        rows[5] *= @sprintf("%12.6f ", stats[ot]["meanrmse"])
                        rows[6] *= @sprintf("%12.6f ", stats[ot]["meanbias"])
                        rows[7] *= @sprintf("%12.6f ", stats[ot]["rmsemean"])
                        rows[8] *= @sprintf("%12.6f ", stats[ot]["biasmean"])
                    end
                    irow = 9
                end
                if i_iter > 1
                    append!(rows, ["change RMSE(i1):", "change bias(i1):"])
                    for ot in obs_type_unique
                        rows[irow]   *= @sprintf("%+12.6f ", stats[ot]["rmse1"]-stats_old[ot]["rmse1"])
                        rows[irow+1] *= @sprintf("%+12.6f ", stats[ot]["bias1"]-stats_old[ot]["bias1"])
                    end
                end
                for row in rows
                    s *= row * "\n"
                end
                s
            end

            stats_old = copy(stats)

            #
            #
            #

            skip_da = false
            if allow_skip_da
                @warn("allow_skip_da is still experimental")
                if i_iter < num_iter
                    # note that first run starts at 000
                    suffix_next = @sprintf("%03d_%03d", i_cycle, i_iter)
                else
                    suffix_next = @sprintf("%03d_%03d", i_cycle+1, 0)
                end
                if isfile(joinpath(config["storagedir"], "$(config["file_prefix"])_mod_$(suffix_next)_001.nc"))
                    println("Found output for next iteration (and \"allow_skip_da\" is active) -- skipping ahead.")
                    skip_da = true
                end
            end

            # create new initial files
            if !iterconstini && (i_iter < num_iter || !noassimlastiter)
                inifiles = [joinpath(config["storagedir"], "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d_%03d", i_cycle, i_iter-1, i)).nc") for i in 1:num_ens_curr]
                #inifiles_new = [joinpath(rundirs[i], "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d_%03d", i_cycle, i_iter, i)).nc") for i in 1:num_ens_curr]
                inifiles_new = [joinpath(rundirs[i], "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d", i_cycle, i_iter)).nc") for i in 1:num_ens_curr]
                if ! skip_da
                    for i in 1:num_ens_curr
                        @info "copying $(inifiles[i]) => $(inifiles_new[i])"
                        cp(inifiles[i], inifiles_new[i], force=true)
                    end
                end
                inifiles = inifiles_new
            end

            if skip_da
                println("Skipping this data assimilation cycle.")
            elseif config["estimationtype"] == "bioparameter"
                if (i_iter < num_iter) || !noassimlastiter
                    parameter_values = config["assimilation_function"](i_cycle, i_iter, rpi, mod_values[:, obsindex_sub], obs_values_sub, parameter_values, obs_info_sub)
                    if size(parameter_values) ≠ (num_ens, length(parameter_names))
                        error("The output of $(config["assimilation_function"]) has invalid shape (expected $((num_ens, length(parameter_names))) but obtained $(size(parameter_values))).")
                    end
                end

                if display_function !== nothing
                    display_function(rpi, parameter_values)
                end

                if any(parameter_values .< 0.0)
                    @warn "Negative parameter(s) found!"
                    for (ip, p) in enumerate(parameter_values)
                        if p < 0.0
                            @warn "   $(parameter_names[ip]): $(p*parameter_values[ip])"
                        end
                    end
                    if config["negative_parameter_response"] == "adjust"
                        # TODO-low allow setting of minval
                        parameter_values .= max.(parameter_values, config["negative_parameter_minval"])
                    elseif config["negative_parameter_response"] ≠ "warn"
                        error("Negative parameter(s).")
                    end
                end
            elseif config["estimationtype"] == "state"
                if (i_iter < num_iter) || !noassimlastiter
                    testinfo = nothing
                    if config["assimilation_function_usefileinput"]
                        # TODO-low remove these checks
                        ncheck = length(glob(joinpath("*", "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d", i_cycle, i_iter)).nc"), config["rundir"]))
                        if ncheck ≠ num_ens
                            throw(ROMSError("Incorrect number of mod-files (found $ncheck, expected $num_ens)."))
                        end
                        ncheck = length(glob("$(config["file_prefix"])_mod_$(suffix)_*.nc", config["storagedir"]))
                        if ncheck ≠ num_ens
                            throw(ROMSError("Incorrect number of ini-files (found $ncheck, expected $num_ens)."))
                        end
                        output_raw = config["assimilation_function"](i_cycle, i_iter, config["obsfile"],
                                            config["storagedir"], "$(config["file_prefix"])_mod_$(suffix)_*.nc",
                                            config["rundir"], joinpath("*", "$(config["file_prefix"])_ini_$(@sprintf("%03d_%03d", i_cycle, i_iter)).nc"),
                                            findfirst(obsindex), findlast(obsindex))
                    else
                        output_raw = config["assimilation_function"](i_iter, num_iter, mod_values[:, obsindex_sub], obs_values_sub, inifiles, obs_info_sub)
                    end
                    if isa(output_raw, Tuple)
                        # for 4DEnOI
                        output, testinfo = output_raw
                    else
                        output = output_raw
                    end
                    if output !== nothing
                        if length(output) > 1 && isfile(output[1])

                            num_testfiles = length(output)

                            @info "Found $(num_testfiles) output files."

                            #
                            # perform runs for test files
                            #

                            rundirs_test = [r * "_test" for r in rundirs[1:num_testfiles]]
                            rfm.rundirs = rundirs_test
                            suffix = @sprintf("test_%03d_%03d", i_cycle, i_iter-1) # note that first run starts at 000
                            rfm.logfilename = "roms_$(suffix).log"

                            set_output_names!(config, suffix)
                            create_files(rfm)

                            set!(rfm, "oceanin", "ININAME", output)

                            start_jobs(rs, rfm)
                            @info "moving model output to \"$(config["storagedir"])\""
                            move_output(rfm, config["storagedir"], filenamemask="*_[mri]*_*.nc")

                            #
                            # evaluate results
                            #

                            modfiles_test = glob("$(config["file_prefix"])_mod_$(suffix)_*.nc", config["storagedir"])
                            rstfiles_test = glob("$(config["file_prefix"])_rst_$(suffix)_*.nc", config["storagedir"])

                            if length(modfiles_test) ≠ num_testfiles
                                throw(ROMSError("Expected $num_testfiles mod files but found $(length(modfiles_test))."))
                            end
                            if length(rstfiles_test) ≠ num_testfiles
                                throw(ROMSError("Expected $num_testfiles rst files but found $(length(rstfiles_test))."))
                            end
                            mod_values_test = fill(NaN, (num_testfiles, sum(obsindex)))

                            @info "extracting $(sum(obsindex)) mod file values"
                            for (imod, modfile) in enumerate(modfiles_test)
                                @info("opening \"$(modfile)\"")
                                NCDatasets.Dataset(modfile) do nc
                                    ## apply minval
                                    #mod_values_test[imod, :] = max.(minval, nc["NLmodel_value"][obsindex])
                                    # if there is a MethodError: Cannot `convert` an object of type Missing to an object of type Float64
                                    # this may indicate a ROMS problem like a blowup
                                    if any(nc["NLmodel_value"][obsindex] .=== missing)
                                        @warn("Found missing values in NLmodel_value, this indicates a blowup.")
                                    else
                                        mod_values_test[imod, :] = nc["NLmodel_value"][obsindex]
                                        tmp = nc["obs_scale"][obsindex]
                                        if !all(obs_scale .=== tmp)
                                            @warn """Found differing number of values obs_scale = 0 in testing run
                                            in reference from regular run: $(sum(obs_scale .=== missing))"
                                            in file \"$(basename(modfile))\":        $(sum(tmp .=== missing))"""
                                            obsindex_sub .&= tmp .!== missing
                                        end
                                        #@assert(all(obs_scale .== tmp), "Problem with mod file \"$(modfile)\", obs_scale values do not agree with reference.")
                                    end
                                end
                            end

                            stats = Dict{Int, Dict{String, Any}}()
                            obs_type_unique = sort(unique(obs_info_sub["obs_type"]))
                            for ot in obs_type_unique
                                cind = obs_info_sub["obs_type"] .== ot
                                stats[ot] = Dict{String, Any}("numobs" => sum(cind),
                                    "rmse" => [sqrt(mean((mod_values_test[i, obsindex_sub][cind].-obs_values_sub[cind]).^2)) for i in 1:num_testfiles],
                                    "bias" => [mean(mod_values_test[i, obsindex_sub][cind])-mean(obs_values_sub[cind]) for i in 1:num_testfiles],
                                    "normcost" => [mean(((mod_values_test[i, obsindex_sub][cind].-obs_values_sub[cind]).^2)./obs_info_sub["obs_error"][cind]) for i in 1:num_testfiles],
                                )
                            end
                            stats[0] = Dict{String, Array{Float64}}("numobs" => length(obs_values_sub),
                                "rmse" => [sqrt(mean((mod_values_test[i, obsindex_sub].-obs_values_sub).^2)) for i in 1:num_testfiles],
                                "bias" => [mean(mod_values_test[i, obsindex_sub])-mean(obs_values_sub) for i in 1:num_testfiles],
                                "normcost" => [mean(((mod_values_test[i, obsindex_sub].-obs_values_sub).^2)./obs_info_sub["obs_error"]) for i in 1:num_testfiles],
                            )

                            # NOTE:
                            # use a non-logtransfomed file (including non-logtransfomed error)

                            cind = obs_info_sub["obs_type"] .<= 7
                            stats[0]["normcost_biolog"] = zeros(num_testfiles)
                            for i in 1:num_testfiles
                                stats[0]["normcost_biolog"][i] = sum(((mod_values_test[i, obsindex_sub][cind].-obs_values_sub[cind]).^2)./obs_info_sub["obs_error"][cind])
                            end
                            cind = cind .== false #(!).cind
                            for i in 1:num_testfiles
                                stats[0]["normcost_biolog"][i] += sum(((log.(max.(minval, mod_values_test[i, obsindex_sub][cind])).-log.(max.(minval, obs_values_sub[cind]))).^2)./obs_info_sub["obs_error"][cind])
                                stats[0]["normcost_biolog"][i] /= length(cind)
                            end

                            for name in keys(stats[0])
                                if name != "numobs"
                                    stats[0][name][stats[0][name] .=== NaN] .= Inf
                                end
                            end

                            i_opt = argmin(stats[0][config["costtype"]::String])

                            @info begin
                                s = "summary:\n"
                                s *= "cost type = " * config["costtype"] * ", "
                                s *= "i_opt = $(i_opt), num_sigma = $(testinfo["numstd"])\n"
                                s *= "number of observations: $(length(obs_values_sub))\n"
                                s *= " run     "
                                s *= "    stepsize "
                                s *= " RMSE: total "
                                for ot in obs_type_unique
                                    s *= @sprintf(" obs type %2d ", ot)
                                end
                                s *= " lcost total:"
                                s *= " cost: total "
                                for ot in obs_type_unique
                                    s *= @sprintf(" obs type %2d ", ot)
                                end
                                s *= "\n"
                                for i in 1:num_testfiles
                                    s *= @sprintf(" test %2d ", i)
                                    s *= @sprintf("%12.6f ", testinfo["stepsize"][i])
                                    s *= @sprintf("%12.6f ", stats[0]["rmse"][i])
                                    for ot in obs_type_unique
                                        s *= @sprintf("%12.6f ", stats[ot]["rmse"][i])
                                    end
                                    s *= @sprintf("%12.6f ", stats[0]["normcost_biolog"][i])
                                    s *= @sprintf("%12.6f ", stats[0]["normcost"][i])
                                    for ot in obs_type_unique
                                        s *= @sprintf("%12.6f ", stats[ot]["normcost"][i])
                                    end
                                    s *= "\n"
                                end
                                s
                            end

                            #
                            # add in new rstfile and update inifile
                            #

                            #@info """copying \"$(rstfiles_test[i_opt])\" to \"$(rstfiles[1])\"
                            #copying \"$(output[i_opt])\" to \"$(inifiles[1])\" """
                            @info "copying \"$(output[i_opt])\" to \"$(inifiles[1])\" "
                            #cp(rstfiles_test[i_opt], rstfiles[1], force=true)
                            #error("Not doing that, terminating instead.")
                            cp(output[i_opt], inifiles[1], force=true)


                            #
                            # cleanup
                            #

                            # TODO: note that all test files in output are currently in the first run directory
                            for f in output
                                if isfile(f)
                                    @info "removing $(f)"
                                    rm(f, force=true)
                                end
                            end
                        else
                            error("Unknown or invalid output of assimilation function.")
                        end
                    end
                    #if refillensemble && length(inifiles) < length(inifiles_save)
                    #    @warn "refilling ini ensemble ($(length(inifiles))->$(length(inifiles_save))); new ensemble:"
                    #    inifiles_save[1:length(inifiles)] = inifiles
                    #    inifiles = inifiles_save
                    #    for i in 1:length(inifiles)
                    #        @info inifiles[i]
                    #    end
                    #end
                end
            end
            if refillensemble && length(rstfiles) < length(rstfiles_save)
                @warn "refilling rst ensemble ($(length(rstfiles))->$(length(rstfiles_save))); new ensemble:"
                rstfiles_save[1:length(rstfiles)] = rstfiles
                rstfiles = rstfiles_save
                for i in 1:length(rstfiles)
                    @info rstfiles[i]
                end
            end
        end
        #move_output(rfm, config["storagedir"], filenamemask="*_ini_*.nc")
    end
    if config["cleanup_on_success"]
        @info "Performing cleanup."
        rm(config["rundir"], recursive=true)
    end
    @info "Successfully completed last cycle."
end

if isinteractive()
    #
    # input args
    #

    println(" entering isinteractive")

    s = ArgParse.ArgParseSettings("Start an an ensemble of ROMS simulations " *
        "for data assimilation.",
        version = "version 0.1", # version info
        add_version = true)      # auto-add version option

    ArgParse.@add_arg_table! s begin
		"--config", "-c"
			default = "romsassim_config.yaml"
			arg_type = String
			help = "the romsassim configuration file"
#         "--debug", "-d"
#            action = :store_true
#            help = "Print some debugging output."
         "--test", "-t"
			arg_type = Int
            default = 0
            help = "Test/debug mode."
#        "--opt1"
#            nargs = '?'              # '?' means optional argument
#            arg_type = Int           # only Int arguments allowed
#            default = 0              # this is used when the option is not passed
#            constant = 1             # this is used if --opt1 is passed with no argument
#            help = "an option"
#        "--karma", "-k"
#            action = :count_invocations  # increase a counter each time the option is given
#            help = "increase karma"
#        "arg1"
#            nargs = 2                        # eats up two arguments; puts the result in a Vector
#            help = "first argument, two " *
#                   "entries at once"
#            required = true
#        "arg2"
#            nargs = '*'                            # eats up as many arguments as possible before an option
#            default = Any["no_arg_given"]          # since the result will be a Vector{Any}, the default must
#                                                   # also be (or it can be [] or nothing)
#            help = "second argument, eats up " *
#                   "as many items as possible " *
#                   "before an option"
    end

    args = ArgParse.parse_args(s)
    println("Parsed args:")
    for (key, val) in args
        println("  $key => $(repr(val))")
    end

    #
    # read config
    #

    if endswith(args["config"], ".config")
        c = ConfParser.ConfParse(args["config"])
        ConfParser.parse_http(c)
        config = Dict()
        for (key, val) in c._data
            if length(val) == 1
                if match(r"^[0-9]+$", val[1]) !== nothing
                    config[key] = parse(Int, val[1])
                elseif val[1] == "true"
                    config[key] = true
                elseif val[1] == "false"
                    config[key] = false
                else
                    config[key] = val[1]
                end
            else
                config[key] = val
            end
        end
    elseif endswith(args["config"], ".yaml") || endswith(args["config"], ".yml")
        config = YAML.load(open(args["config"]))
    else
        throw(ConfigurationError("Configuration file format not yet supported."))
    end

    for (key, val) in config
        @info "$key => $(repr(val))"
    end

    romsassim(config)
    nothing
end # isinteractive

end # module
