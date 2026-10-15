## layout

| folder | content |
|---|---|
| `code/` | main analysis chain (all stages: `itrusst_rerun_all.m`, see top-level README): phantoms, calibration, `PRESTUS_paramvar_00-04`, `PRESTUS_pipeline_01-02`, `PRESTUS_thickness_01-07`, setup/range figures |
| `code/functions/` | shared helpers called by the scripts (no script runs from here) |
| `code/paper/` | panels of the paper figures (`paper_fig*`, `paper_supp_fig*`), stage `itrusst_rerun_all('plot-paper')`; output in `figures/paper/fig<N>/` |
| `code/explore/` | side analyses not (yet) used: `PRESTUS_paramvar_06-10/13`, skull sound-speed case study, sub-303 sweep, multiskull, per-voxel `plot_sub*` |
| `code/debug/` | numerical checks (0.25 mm / 0.5 mm / 3D resolution, 2D vs 3D, calibration 2D vs 3D, benchmark 2 mismatch) and stale scripts with hard-coded paths |

Scripts must be run as files (they resolve `rootpath` from their own location) and add `code/functions` to the path themselves.

Figures are written to `figures/<family>/` with the same split: `paper/`, `tissue_ranges/`, `setups/`, `pipeline/`, `paramvar/`, `thickness/`, `explore/<topic>/`, `debug/<topic>/`. Design tables of multi-run designs are written next to their simulation output in `data/tussim/<variant>/`.

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

### paper panels (`code/paper/`, helpers `paper_*.m` in `code/functions/`)

Panels at print size (8 pt Helvetica, axes boxes fixed in cm via `paper_style.m`), saved as vector PDF + 300 dpi PNG by `paper_save.m` for assembly in Illustrator. Run them all with `itrusst_rerun_all('plot-paper')` (independent of the slower diagnostic `plot` stage). Schematics, icons, legends and frames are added there. BabelBrain panels (benchmark 2 = sub-002, collaborator results in `data/pichardo/`, see `load_babelbrain_benchmark2.m`) are drawn next to their PRESTUS counterparts on a shared y-range; `paper_placeholder.m` (empty panel with the neighbouring axes) remains for collaborator panels without data.

| script | panels |
|---|---|
| `paper_fig1_property_ranges` | Fig. 1: property ranges per tissue (`data/tissue_parameters.csv`; attenuation at 500 kHz in dB/cm) |
| `paper_fig1_property_ranges_bone` | Fig. 1 variant: literature values per tissue class, skull split into cortical / trabecular / layer unclear (ITRUSST acoustic-properties table, `measurements.csv`, status 'included') |
| `paper_fig23_axial_profiles` | Figs. 2/3: vertical on-axis I_SPPA / heat-rise profiles, sub-001 vs sub-002, mean and range across paramvar runs |
| `paper_fig234_paramvar_panels` | Figs. 2-4: paramvar range panels with reference-medium lines (`paper_reference_metrics.m`) + BabelBrain counterparts of the sub-002 panels (`bb_ppw`: 9 PPW run by default, 6 PPW as convergence check) |
| `paper_fig3_heating_difference_maps` | Fig. 3: heat-rise maps top - bottom level for skull absorption fraction and attenuation (sub-002) |
| `paper_fig4_thickness` | Fig. 4: brain I_SPPA vs skull thickness |
| `paper_fig5_setups` | Fig. 5: parallel plots of Setup A/B (= setup1/setup2) and toolbox defaults; brain I_SPPA and temperature rise, 3D benchmarks |
| `paper_supp_fig1_focus` | Supp. Fig. 1: focal distance / -3 dB volume range panels; on-axis I_SPPA and -3 dB focal masks for low/high skull sound speed (sub-002) |
| `paper_supp_fig2_density_thickness` | Supp. Fig. 2: skull temperature rise vs thickness for low/high skull density (`itrusst_protocol1a_density_thickness`) |

Two inputs are not tracked in the repository: the BabelBrain CSVs (default `data/pichardo/`, override with `setenv('ITRUSST_BABELBRAIN_DIR', '<folder>')`) and the ITRUSST acoustic-properties table (default `../itrusst_acousticproperties/literature_web/data/measurements.csv`, override with `setenv('ITRUSST_MEASUREMENTS_CSV', '<file>')`). Scripts that need them stop with an error naming the missing path and the variable (`itrusst_input_path.m`).

Thermal reference lines need a reference-medium (`_baseline`) run, which `PRESTUS_paramvar_00_baseline_sim.m` provides for sub-002 only; for sub-001 and sub-352 only the acoustic reference is drawn.
