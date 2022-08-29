# ROMSEnsemble.jl

**ROMSEnsemble.jl** is used to run ocean model ensembles and perform ensemble-based data assimilation using the Regional Ocean Modeling System ([ROMS](https://www.myroms.org/)). It automates starting and stopping of ROMS simulations, modifying input files and extracting output, and it requires an existing ROMS setup and a computer with a workload manager, such as [Slurm](https://slurm.schedmd.com/). **ROMSEnsemble.jl** does not provide any data assimilation capabilities, that is, it extracts the necessary state and parameter information from a model ensemble and passe it to a user-supplied data assimilation function.

**ROMSEnsemble.jl** is based on and expands upon this [MATLAB code](https://github.com/bwang63/EnKF_3D_github).

## Getting started

**ROMSEnsemble.jl** requires a working ROMS setup, including a ROMS executable compiled with the `VERIFICATION` option, so that ROMS accepts ROMS observations files as input and so that it outputs the model state at the observation locations. With that in place,  **ROMSEnsemble.jl** is configured using a configuration YAML file (or, alternatively, a configuration contained in a Julia `Dict{String,Any}`) and a modified Slurm sbatch script if `sbatch` is used to run ROMS using Slurm.

Here is example configuration file:
```
%YAML 1.2
---
# A description of all configuration parameters is shown in the Wiki.

# The ROMS executable to use for the ensemble simulations.
executable: "$HOME/romsensemble/romsM"

# ROMS main text-based input file (previously referred to as the "ocean.in" file). 
# Note that this is a template file that will be copied, the original will remain unchanged by ROMSAssim.jl. 
# The `paramchanges_ocean` parameter (below) can be used to make changes to the ROMS parameters in all 
# copies.
ocean_in: "$HOME/romsensemble/romsinfiles/ocean.in"

# Changes to the template `ocean_in` file (optional). Each key must correspond to a parameter in the file 
# specified by `ocean_in`.
# These modifications are only applied to copies of the file specified by `ocean_in`, the original file will
# remain unchanged.
paramchanges_ocean:
    TITLE: "ROMSEnsemble run"
    NtileI: 4
    NtileJ: 5
    NHIS: 96
    DT: 300

# ROMS biological input file containing parameters for the biological model (optional: not required for 
# physics-only run without a coupled biological model).
# Note that this is a template file that will be copied, the original will remain unchanged by ROMSAssim.jl.
bio_in: "$HOME/romsensemble/romsinfiles/npzd_iron.in"

# ROMS data assimilation input file containing parameters for configuring data assimilation, including the
# observation file.
# Note that this is a template file that will be copied, the original will remain unchanged by ROMSAssim.jl.
s4dvar_in: "$HOME/romsensemble/romsinfiles/s4dvar.in"

# The ROMS observation file which includes observations for the length of the full data assimilation run.
# ROMSAssim.jl slices this file and provides only the observations belonging to the current data assimilation 
# cycle to the assimilation function.
# The original file will be compied and not modified by ROMSAssim.jl.
obsfile: "$HOME/romsensemble/data/observationfile.nc"

# The ROMSAssim.jl run directory, used to temporarily store copies of ROMS input files and ROMS output files. 
# Relevant output files will be copied to `storagedir` at the end of each assimilation cycle.
# Unless it is specified as an absolute path, `rundir` is considered to be leative to Julia's current working 
# directory. 
# This directory will be created if it does not exist and ROMSAssim.jl may delete and overwrite files in this 
# directory!
rundir: "run"

# The ROMSAssim.jl storage directory, used to store ROMS output files. Relevant output files will be copied
# from `rundir` at the end of each assimilation cycle.
# This directory will be created if it does not exist and ROMSAssim.jl may delete and overwrite files in this 
# directory!
storagedir: "output"

# A prefix that is appended to the filename of ROMS output files.
file_prefix: "romsens"

# The initial conditions used for the ROMSEnsemble.jl ensemble run.
# These can be specified as a single file name (all ensemble members start with the same initial conditions),
# a glob-expression to be expanded (e.g. "initial_conditions/ini*.nc"), or an array containing the individual 
# file names (e.g. ["ic/ini01.nc", "ic/ini02.nc", "ic/ini03.nc"]). Note that in the latter two cases, the 
# number of files must match the size of the ensemble `n_ens` (below).
initial_conditions: "$HOME/romsensemble/initial_conditions/ini*.nc"

# The length of a spinup in units of days, set to 0 for no spinup. 
spinup_days: 0

# The number of ensemble members.
n_ens: 5

# The function performing the data assimilation. See the function stubs 
# `ROMSEnsemble.af_stateestimation_stub`
# for reference.
assimilation_function: "ROMSEnsemble.af_stateestimation_stub"

# The length of each data assimilation cycle in days.
# Note that instead of specifying `cycle_length` and `num_cycles`, the stop dates between cycles can also 
# be specified with the `stopdates` option.
cycle_length: 4

# The number of data assimilation cycle to run.
num_cycles: 1

# The data assimilation type: use "state" for state estimation, and "parameter" for parameter estimation.
estimationtype: "state"

# The cycle setup, determining the number of iterations or "outer loops" to perform in each cycle and how 
# many of the ensemble members to run. This can be a `String` or a vector of integers.
# "EnKF" is equivalent to [`n_ens`, `n_ens`], i.e. run the full ensemble forward, perform data assimilation,
# then run the full (updated) ensemble forward to the end of the cycle.
# When `cycle_setup` is specified by a vector, the `noassimlastiter` variable determines if an assimilation
# step is performed in the last iteration (default: perform no assimilation in last iteration 
# `noassimlastiter=true`)
cycle_setup: "EnKF"

# The variables to pass from the observation file to the assimilation function. These provide additional 
# information about the observations, for example, the observation errors or observations locations, useful 
# for implementing localization.
obsfile_variables: ["obs_type", "obs_error", "obs_Xgrid", "obs_Ygrid", "obs_Zgrid", "obs_depth"]

# The way to start each ROMS simulation of the ensemble. The strings "sbatch" and "srun" indicate the use of
# Slurm's sbatch and srun commands (see https://slurm.schedmd.com/).
# Alternatively, an instance of type `ROMSEnsemble.ROMSStarter` can be added to a compiled configuration 
# dictionary in order to use a custom ROMSStarter.
ROMSStarter: "sbatch"

# The "sbatch" option for `ROMSStarter` requires the specification of a template sbatch script. The structure 
# for this file is described in the wiki.
batchscript: "$HOME/romsensemble/data/sbatch_script.template"
```

## Issues and Discussions

Please use the issues page to report issues and the discussions page for discussions. **ROMSEnsemble.jl** was initially created and is still being used for my personal research, and I would be happy to expand it, if there is interest.
