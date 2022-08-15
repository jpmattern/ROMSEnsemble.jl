#
# general ROMSStarter
#

abstract type ROMSStarter end

function start_jobs(rs::ROMSStarter, rfm::ROMSFileManager; numjobs::Int=-1, monitor::Bool=true, maxnumrestarts::Int=5, commandargs::Union{Nothing, String}=nothing)
    if numjobs < 1
        numjobs = length(rfm.rundirs)
    else
        @assert(numjobs<=length(rfm.rundirs), "Input \"numjobs\" cannot exceed the number of run directories ($(length(rfm.rundirs))).")
    end

    for i in 1:numjobs
        d = rfm.rundirs[i]
        @debug("directory $d")
        start_job(rs, i, d, rfm.localnames["oceanin"], rfm.logfilename, commandargs=commandargs)
    end

    if monitor
        numrestarts = zeros(Int, numjobs)
        i = 1
        while i <= numjobs
            wait_for(rs, i)
            @info("Run $(i) terminated.")
            if maxnumrestarts > 0
                d = rfm.rundirs[i]
                romslogfile = joinpath(d, rfm.logfilename)
                @debug "   checking $romslogfile"
                if run_finished(romslogfile)
                    i += 1
                else
                    if numrestarts[i] < maxnumrestarts
                        @warn "Run $(i) terminated with error -- attempting restart."
                        start_job(rs, i, d, rfm.localnames["oceanin"], rfm.logfilename, commandargs=commandargs)
                        numrestarts[i] += 1
                        sleep(10)
                    else
                        @error "Run $(i) terminated with error -- maximum number of restarts ($(maxnumrestarts)) reached."
                        i += 1
                        break
                    end
                end
            end
        end
    end
    nothing
end


#
# helper functions
#

function run_finished(romslogfile::String)
    # look for "ROMS/TOMS: DONE..."
    if !isfile(romslogfile)
        return false
    end
    # 2020-03-29: there was a Julia issue where
    # a return statement in the loop would not exit the function
    # I am now programming around this
    result = false
    open(romslogfile) do f
        for line in eachline(f)
            if occursin("ROMS/TOMS: DONE", line)
                result = true
                break
            elseif occursin("Bus error", line)
                @error("Bus error detected:\n$(line)")
                break
            elseif occursin("WRT_RST", line) # TODO-low maybe remove this. At the moment needed because some runs may fail past end of ROMS output
                #@info("Found final WRT_RST.")
                result = true
                break
            end
        end
    end
    return result
end

#
# SrunROMSStarter
#

mutable struct SrunROMSStarter <: ROMSStarter
    exe :: String
    ntasks :: Int
    sleeptime :: Float64
    pids :: Array{Any, 1}
    function SrunROMSStarter(numruns::Int, exe::String, ntasks::Int; sleeptime::Float64=0.0)
        @assert(numruns > 0, "Number of runs must be greater than 0.")
        @assert(isfile(exe), "ROMS executable must exist.")
        @assert(ntasks > 0, "Number of tasks must be greater than 0.")
        rs = new()
        rs.exe = exe
        rs.ntasks = ntasks
        rs.sleeptime = sleeptime
        rs.pids = Array{Any, 1}(undef, numruns)
        rs
    end
end

function Base.show(io::IO, rs::SrunROMSStarter)
    print(io, "SrunROMSStarter")
end

function start_job(rs::SrunROMSStarter, irun::Int, rundir, oceaninfilename, logfilename, commandargs::Union{Nothing, String}=nothing)
    cmd = `srun --chdir $(rundir) --ntasks $(rs.ntasks) $(rs.exe) $(oceaninfilename)`
    @debug "cmd=$(cmd)"
    rs.pids[irun] = run(pipeline(cmd, stdout=joinpath(rundir, logfilename), stderr=joinpath(rundir, "stderr.txt") ), wait=false)
    sleep(rs.sleeptime)
    nothing
end

function wait_for(rs::SrunROMSStarter, irun::Int)
    Base.wait(rs.pids[irun])
end

#
# SBatchROMSStarter
#

mutable struct SBatchROMSStarter <: ROMSStarter
    exe :: String
    ntasks :: Int
    sleeptime :: Float64
    jids :: Array{Int, 1}
    sbatchtemplatefile :: String
    jobnameprefix :: String
    regex_matchjid :: Regex
    function SBatchROMSStarter(numruns::Int, exe::String, ntasks::Int, sbatchtemplatefile::String; sleeptime::Float64=0.0, jobnameprefix::String="SBatchROMSStarter", regex_matchjid::Regex=r".*job +([0-9]+)")
        @assert(numruns > 0, "Number of runs must be greater than 0.")
        @assert(isfile(exe), "ROMS executable must exist.")
        @assert(ntasks > 0, "Number of tasks must be greater than 0.")
        @assert(isfile(sbatchtemplatefile), "The sbatch template (\"$(sbatchtemplatefile)\") file must exist.")
        rs = new()
        rs.exe = exe
        rs.ntasks = ntasks
        rs.sbatchtemplatefile = sbatchtemplatefile
        rs.sleeptime = sleeptime
        rs.jids = Array{Int, 1}(undef, numruns)
        rs.jobnameprefix = jobnameprefix
        rs.regex_matchjid = regex_matchjid
        rs
    end
end

function Base.show(io::IO, rs::SBatchROMSStarter)
    print(io, "SBatchROMSStarter")
end

function start_job(rs::SBatchROMSStarter, irun::Int, rundir::String, oceaninfilename::String, logfilename::String; commandargs::Union{Nothing, String}=nothing)
    sbatchjobname = rs.jobnameprefix*"_"*@sprintf("%03d", irun)

    # copy template file to rundir
    sbatchfile = joinpath(rundir, "sbatch.bash")
    cp(rs.sbatchtemplatefile, sbatchfile, force=true)

    # temporarily added --exclude=node27, node30, node31, node32
    # temporarily added --exclude=node01, node08, node09, node10, node11
    # temporarily added --partition=nrt
    if commandargs === nothing
        cmd = `sbatch --chdir $(rundir) --job-name="$sbatchjobname" --export="ALL,ROMS_EXECUTABLE=$(rs.exe),ROMS_INFILE=$(oceaninfilename)" --ntasks $(rs.ntasks) --output="$(logfilename)" $(sbatchfile)`
    else
        cmd = `sbatch --chdir $(rundir) --job-name="$sbatchjobname" --export="ALL,ROMS_EXECUTABLE=$(rs.exe),ROMS_INFILE=$(oceaninfilename)" --ntasks $(rs.ntasks) --output="$(logfilename)" $(commandargs) $(sbatchfile)`
    end
    @info "submitting job \"$(sbatchjobname)\""
    tmp = read(pipeline(cmd, stderr=joinpath(rundir, "sbatcherr.txt") ), String)
    try
        # expecting something like "Submitted batch job 1234"
        rs.jids[irun] = parse(Int, match(rs.regex_matchjid, tmp).captures[1])
    catch
        @warn "Failed to obtain job ID for submitted job."
        rs.jids[irun] = -1
    end
    sleep(rs.sleeptime)
    nothing
end

function wait_for(rs::SBatchROMSStarter, irun::Int)
    cmd = `squeue -j $(rs.jids[irun]) --noheader`
    try
        while length(read(pipeline(cmd))) > 0
            sleep(5)
        end
    catch e
        # slurm_load_jobs error: Invalid job id specified
        if isa(e, ProcessFailedException)
            # older job IDs may be forgotten
            return
        else
            rethrow(e)
        end
    end
    return
end

