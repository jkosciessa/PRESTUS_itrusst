# **ITRUSST benchmark simulations**

The goal of these analyses is to explore the impact of tissue medium parameter assumptions on relevant outcomes of acoustic intensity and heating.
This simulation uses PRESTUS, a wrapper for k-wave, to perform acoustic and heating simulations.
Different sequences and parameter setups are specified either via dedicated config files, or within the MATLAB call scripts.
Analyses are run on a GPU-supported HPC cluster hosted at the Donders Institute.

For an overview of the goals and ongoing results (*access restricted*), see [here](https://docs.google.com/document/d/16yBCTZDG966979RLcoQdo9XqpzFn1oSb-pxK2KgR4dM/edit?pli=1&tab=t.0).

There are three simulation goals here:

[1] Impact of individual parameter variation on outcomes (all else being equal)

Relevant scripts are called ```PRESTUS_paramvar```.

[2] Intercomparison of different tools 

Relevant scripts are called ```PRESTUS_pipeline```.

[3] Variations of skull thickness layer 

Relevant scripts are called ```PRESTUS_thickness```.

Shared functions live in `code/functions/`, side analyses in `code/explore/` and numerical checks in `code/debug/` (see [code/readme.md](code/readme.md)). Methods are documented in [doc/PRESTUS_itrusst_methods.md](doc/PRESTUS_itrusst_methods.md). PRESTUS is included as submodule `tools/PRESTUS` (cloned from `../v05/tools/PRESTUS`; update with `git -C tools/PRESTUS pull v05_local development`).

Benchmark phantoms (`code/createPhantom.m`) follow the PRESTUS examples, see [here](https://github.com/Donders-Institute/PRESTUS/blob/development/examples/createPhantom.m). The code uses PRESTUS' ```phantom``` implementation (see [here](https://github.com/Donders-Institute/PRESTUS/blob/development/documentation/doc_medium.md)). 

## Rerunning after a change of tissue parameters

Tissue properties live in `data/tissue_parameters.csv` (bounds) and `data/tissue_setups.csv` (which bound setup 1/2 use); the `config_medium_<setup>.yaml` files are generated from them. After editing the CSVs:

```matlab
addpath('code/functions');
itrusst_backup_outputs(pwd, '<yymmdd>');  % data/tussim, data/calibration -> *_<yymmdd>, figures/* -> figures/backup/figures_<yymmdd>/
cd code
itrusst_rerun_all('prepare')    % local: medium configs, thickness phantoms, range/setup figures
itrusst_rerun_all('calibrate')  % cluster, blocking: free-water calibration -> protocol configs
itrusst_rerun_all('simulate')   % cluster: submits all simulations
itrusst_rerun_all('plot')       % once all jobs have finished
```

The backup matters: several simulations use `overwrite_files: never` and would otherwise skip conditions with existing outputs. Side analyses in `code/explore/` and `code/debug/` are not part of the rerun.
