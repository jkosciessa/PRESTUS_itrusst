## layout

| folder | content |
|---|---|
| `code/` | main analysis chain (all stages: `itrusst_rerun_all.m`, see top-level README): phantoms, calibration, `PRESTUS_paramvar_00-04`, `PRESTUS_pipeline_01-02`, `PRESTUS_thickness_01-07`, setup/range figures |
| `code/functions/` | shared helpers called by the scripts (no script runs from here) |
| `code/explore/` | side analyses not (yet) used: `PRESTUS_paramvar_06-10/13`, skull sound-speed case study, sub-303 sweep, multiskull, per-voxel `plot_sub*` |
| `code/debug/` | numerical checks (0.25 mm / 0.5 mm / 3D resolution, 2D vs 3D, calibration 2D vs 3D, benchmark 2 mismatch) and stale scripts with hard-coded paths |

Scripts must be run as files (they resolve `rootpath` from their own location) and add `code/functions` to the path themselves.

Figures are written to `figures/<family>/` with the same split: `tissue_ranges/`, `setups/`, `pipeline/`, `paramvar/`, `thickness/`, `explore/<topic>/`, `debug/<topic>/`. Design tables of multi-run designs are written next to their simulation output in `data/tussim/<variant>/`.

## model overview

sub-001 | soft-tissue only, 64 mm focus (2D axisymmetric, 0.5 mm)
sub-002 | soft-tissue + skull, 64 mm focus (2D axisymmetric, 0.5 mm)
sub-003 | soft-tissue + skull with cortical/trabecular layers (multiskull, explore)

sub-101/102 | 3D counterparts of sub-001/002 (pipeline comparison, 3D benchmarks)
sub-201/202 | 0.25 mm counterparts of sub-001/002 (debug)
sub-301..312 | skull thickness 1..12 mm
sub-351..361 | skull thickness 1.5..11.5 mm (see `itrusst_thickness_subjects.m`)

## shared utilities (`code/functions/`)

### `setup_itrusst_paths.m`

```matlab
pn = setup_itrusst_paths(rootpath, user)
```

Sets up `pn.prestus` (`tools/PRESTUS`), `pn.kwave`, `pn.minimize`, `pn.nifti`, `pn.configs`, `pn.data`, `pn.path_seg`, `pn.code` and (for known users) `pn.simnibs`, and adds them to the MATLAB path. Call after `addpath(fullfile(rootpath, 'code', 'functions'))`.

### `load_paramvar_results.m`

```matlab
[results, xlabels] = load_paramvar_results(data_dir, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols)
```

Reads `sub-XXX_phantom<affix>.csv` per tissue x property x sweep level and returns `results` with fields `sweep`, `Ipa_target`, `maxT` (rise from 37 °C), `Isppa_brain`, `riseT_brain`, `Isppa_skull`, `riseT_skull`, `real_focal_distance_mm`, `focal_vol_3dB_mm3` (each `[n_sweep x n_conds]`) and condition labels `xlabels` (TeX symbols if `property_symbols` is given). Missing files/columns are warned and left as `NaN`.

### `load_focus_heating_metrics.m`

Same inputs as `load_paramvar_results`; reads the heatrise NIfTIs at the empirical focus voxel (`riseT_focus`) and at the fixed nominal target (`riseT_target`).

### `load_babelbrain_benchmark2.m`

```matlab
[results, ref, xlabels] = load_babelbrain_benchmark2(data_dir, ppw, tissue_labels, property_labels, property_symbols)
```

Reads the BabelBrain paramvar results of benchmark 2 (Sam Pichardo; newest `data/pichardo/Benchmark2-<date>-<ppw>PPW.csv`, 6 or 9 points per wavelength) in the layout of `load_paramvar_results` (rows LOW/HIGH = sweep level 1/2). Fields `Isppa_brain`, `MI_brain`, `riseT_brain`, `riseT_skull`, `riseT_target` (fixed target), `CEM43_brain`, `CEM43_skull`; `ref` holds the `Reference` row. PRESTUS vs BabelBrain 6/9 PPW differences: `PRESTUS_paramvar_02_benchmark2_babelbrain_plot.m` (`figures/paramvar/babelbrain/`).

### `plot_metric_with_patches.m`, `make_parallel_plot.m`

Range plot (min-max patch per condition, bottom/top levels as dots; brain grey, skull red) and parallel-coordinates plot of the setup property tables.

### `itrusst_calibrate.m`, `itrusst_apply_calibration.m`, `itrusst_thickness_subjects.m`

Free-water calibration of the annular transducer to a target intensity; copy of the calibrated drive (`annular: {...}` block of `data/calibration/PRESTUS_virtual_parameters/<setup>-F*mm-I*wpercm2.yaml`) into `config_<setup>.yaml`; mapping between skull thickness (mm) and thickness-phantom subject IDs.

### `itrusst_tissue_properties.m`, `itrusst_write_medium_configs.m`

Tissue properties of a setup (`setup1`, `setup2`, `reference`, `bottom`, `top`) from `data/tissue_parameters.csv` + `tissue_setups.csv`; generation of `data/configs/config_medium_<setup>.yaml` from them.

### `itrusst_backup_outputs.m`

Moves `data/tussim`, `data/calibration` and `figures/*` aside under a date tag before a rerun, keeping the I/O paths.
