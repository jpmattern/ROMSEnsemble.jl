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

"""
Read ROMS log file and return whether simulation finished without errors.
"""
function run_finished(romslogfile::String) :: Bool
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

"""
    SrunROMSStarter(numruns::Int, exe::String, ntasks::Int; sleeptime::Float64=0.0)

A ROMSStarter using the workload manager Slurm's `srun` command to start ROMS.
"""
mutable struct SrunROMSStarter <: ROMSStarter
    exe :: String
    ntasks :: Int
    sleeptime :: Float64
    pids :: Array{Any, 1}
    function SrunROMSStarter(numruns::Int, exe::String, ntasks::Int; sleeptime::Float64=0.0)
        if numruns ≤ 0
            error("Number of runs must be greater than 0.")
        end
        if ! isfile(exe)
            error("ROMS executable must exist.")
        end
        if ntasks ≤ 0
            error("Number of tasks must be greater than 0.")
        end
        return new(exe, ntasks, sleeptime, Array{Any, 1}(undef, numruns))
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
# SbatchROMSStarter
#

"""
    SbatchROMSStarter(numruns::Int, exe::String, ntasks::Int, sbatchtemplatefile::String; sleeptime::Float64=0.0, jobnameprefix::String="SbatchROMSStarter", regex_matchjid::Regex=r".*job +([0-9]+)")

A ROMSStarter using the workload manager Slurm's `sbatch` command to start ROMS.
"""
mutable struct SbatchROMSStarter <: ROMSStarter
    exe :: String
    ntasks :: Int
    sleeptime :: Float64
    jids :: Array{Int, 1}
    sbatchtemplatefile :: String
    jobnameprefix :: String
    regex_matchjid :: Regex
    function SbatchROMSStarter(numruns::Int, exe::String, ntasks::Int, sbatchtemplatefile::String; sleeptime::Float64=0.0, jobnameprefix::String="SbatchROMSStarter", regex_matchjid::Regex=r".*job +([0-9]+)")
        if numruns ≤ 0
            error("Number of runs must be greater than 0.")
        end
        if ! isfile(exe)
            error("ROMS executable must exist.")
        end
        if ntasks ≤ 0
            error("Number of tasks must be greater than 0.")
        end
        if ! isfile(sbatchtemplatefile)
            error("The sbatch template (\"$(sbatchtemplatefile)\") file must exist.")
        end
        rs = new()
        rs.exe = exe
        rs.ntasks = ntasks
        rs.sbatchtemplatefile = sbatchtemplatefile
        rs.sleeptime = sleeptime
        rs.jids = Array{Int, 1}(undef, numruns)
        rs.jobnameprefix = jobnameprefix
        rs.regex_matchjid = regex_matchjid
        return rs
    end
end

function Base.show(io::IO, rs::SbatchROMSStarter)
    print(io, "SbatchROMSStarter")
end

function start_job(rs::SbatchROMSStarter, irun::Int, rundir::String, oceaninfilename::String, logfilename::String; commandargs::Union{Nothing, String}=nothing)
    sbatchjobname = rs.jobnameprefix*"_"*@sprintf("%03d", irun)

    # copy template file to rundir
    sbatchfile = joinpath(rundir, "sbatch.bash")
    cp(rs.sbatchtemplatefile, sbatchfile, force=true)

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
        @debug "job ID: $(rs.jids[irun])"
    catch
        @warn "Failed to obtain job ID for submitted job."
        rs.jids[irun] = -1
    end
    sleep(rs.sleeptime)
    nothing
end

function wait_for(rs::SbatchROMSStarter, irun::Int)
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
    nothing
end

