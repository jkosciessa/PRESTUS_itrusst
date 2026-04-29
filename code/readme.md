## model overview

sub-001 | soft-tissue only, 64 mm focus
sub-002 | soft-tissue + skull, 64 mm focus
sub-003 | soft-tissue + skull (thickness variations), 64 mm focus

sub-004 | soft-tissue + skull, 80 mm focus
sub-005 | soft-tissue + skull, 80 mm focus

## shared utilities

### `load_paramvar_results.m`

Loads parameter variation simulation results from CSV output files.

```matlab
[results, xlabels] = load_paramvar_results(data_dir, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values)
```

Reads `sub-XXX_phantom_output_table<affix>.csv` files produced by the simulation loop and returns:
- `results` — struct with fields `sweep`, `Ipa_target`, `maxT` (temperature rise from 37 °C), `Isppa_brain`, `riseT_brain`, `Isppa_skull`, `riseT_skull`, `real_focal_distance_mm`; each `[n_sweep × n_conds]`
- `xlabels` — cell array of condition labels (e.g. `'b_sound_speed'`, `'s_density'`)

Missing columns in a CSV are skipped with a warning; missing files are also warned and left as `NaN`.

Used by: `PRESTUS_paramvar_02_benchmark1_plot.m`, `benchmark2_plot.m`, `benchmark303_plot.m`, `benchmark_allvariants_plot.m`.

---

### `setup_itrusst_paths.m`

Initialises the PRESTUS toolbox path struct and adds all required directories to the MATLAB path.

```matlab
pn = setup_itrusst_paths(rootpath, user)
```

Sets up `pn.prestus`, `pn.prestus_fun`, `pn.prestus_tools`, `pn.kwave`, `pn.minimize`, `pn.nifti`, `pn.configs`, `pn.data`, `pn.path_seg`, `pn.code`, and (for known users) `pn.simnibs`. Must be called after `addpath(fullfile(rootpath, 'code'))` so MATLAB can find the function.

Used by: `PRESTUS_paramvar_01_benchmarks_sim.m`, `PRESTUS_pipeline_01_benchmark1_sim.m`.
