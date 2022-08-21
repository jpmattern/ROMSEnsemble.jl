# ROMSEnsemble.jl

**ROMSEnsemble.jl** is used to run ocean model ensembles and perform ensemble-based data assimilation using the Regional Ocean Modeling System ([ROMS](https://www.myroms.org/)). It automates starting and stopping of ROMS simulations, modifying input files and extracting output, and it requires an existing ROMS setup and a computer with a workload manager, such as [Slurm](https://slurm.schedmd.com/). **ROMSEnsemble.jl** does not provide any data assimilation capabilities, that is, it extracts the necessary state and parameter information from a model ensemble and passe it to a user-supplied data assimilation function.

**ROMSEnsemble.jl** is based on and expands upon this [MATLAB code](https://github.com/bwang63/EnKF_3D_github).

## Issues and Discussions

Please use the issues page to report issues and the discussions page for discussions. **ROMSEnsemble.jl** was initially created and is till being used for my personal research, and I would be willing expand it, if there is interest.
