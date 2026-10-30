# ITRUSST benchmark simulations

Simulations of how assumptions about tissue acoustic and thermal properties change predicted transcranial ultrasound intensity and heating. They use [PRESTUS](https://github.com/Donders-Institute/PRESTUS) (k-Wave acoustics and Pennes bioheat thermal modelling) on digital phantoms: tissue-label volumes generated in code (`code/createPhantom.m`, `code/PRESTUS_thickness_01_createPhantoms.m`), not physical phantoms or imaging data. Subject IDs (`sub-XXX`) label phantom geometries.

## Requirements

- MATLAB with the Image Processing Toolbox (`niftiread`).
- PRESTUS as submodule `tools/PRESTUS` (k-Wave included): clone with `git clone --recurse-submodules`, or run `git submodule update --init --recursive` after cloning.
- Optional: a SLURM or qsub scheduler (used automatically when found) and an NVIDIA GPU for the 3D setup comparison.
- Supp. Fig. 3a reads the ITRUSST acoustic-properties table (`measurements.csv`), which is not part of this repository; point `ITRUSST_MEASUREMENTS_CSV` to a copy.

## Layout

| path | content |
|---|---|
| `code/` | simulation and plot scripts of the three analyses, `itrusst_rerun_all.m` (all stages) |
| `code/functions/` | shared functions |
| `code/paper/` | paper figure panels and the supplementary table (see [code/readme.md](code/readme.md)) |
| `data/tissue_parameters.csv` | low / reference / high value per tissue and property (source: `data/tissue_parameters_source.md`) |
| `data/tissue_setups.csv` | bound per property for each literature setup |
| `data/configs/` | protocol config and medium configs (`config_medium_<setup>.yaml`; generated for the literature setups, hand-written for the toolbox defaults) |
| `data/literature_skull-attenuation_values.csv` | reported skull attenuation vs frequency (Supp. Fig. 3) |

Simulation outputs (`data/tussim/`), phantoms (`data/simnibs/`) and figures (`figures/`) are generated and not tracked.

## Analyses

| analysis | phantoms | scripts |
|---|---|---|
| 1: one-at-a-time (OAAT) sensitivity | sub-001, sub-002, sub-352 (2D) | `PRESTUS_paramvar_00-04` |
| 2: skull thickness | sub-301–312, sub-351–361 (2D) | `PRESTUS_thickness_01-07` |
| 3: setup and toolbox default comparison | sub-101, sub-102 (3D) | `PRESTUS_pipeline_01-02` |

**1. OAAT sensitivity.** Each of seven properties per tissue (sound speed, density, attenuation, thermal conductivity, specific heat capacity, perfusion, absorption) is set to its low or high value while all others stay at reference. Benchmark 1 (water / soft tissue, sub-001) varies brain properties only; Benchmark 2 (water / skull / brain, sub-002) and the most transmissive skull thickness (2.5 mm, sub-352) vary brain and skull properties. Output: `data/tussim/itrusst_protocol1a/`.

**2. Skull thickness.** 23 skull phantoms of uniform thickness from 1 to 12 mm in 0.5 mm steps, all properties at reference, acoustic simulation only. Two follow-ups repeat the thickness series at low, reference and high skull density, and sweep skull sound speed on the 3, 6.5 and 9 mm phantoms (both with a shortened 20 s thermal protocol).

**3. Setup and toolbox default comparison.** Four literature setups built from the OAAT directions: Temp-H / Temp-L take, per property, the bound that maximises / minimises brain heating; Trans-H / Trans-L take, per acoustic property, the bound that maximises / minimises brain intensity and are simulated without heating. They are compared with the defaults of PRESTUS, BabelBrain and k-Plan on 3D versions of both benchmarks (axisymmetric acoustics on GPU, thermal simulation on the radially expanded 3D field). The properties of all setups are listed by `code/paper/paper_supp_table_setups.m`.

## Shared settings

- Transducer: single-element focused bowl (64 mm aperture and radius of curvature), 500 kHz, calibrated in free water to I_SPPA = 60 W/cm² (`itrusst_calibrate_all.m`).
- Grid: 0.5 mm, 70 × 120 mm active domain plus a 20 mm PML; the Benchmark 2 skull (sub-002) is a spherical cap (outer apex at 30 mm depth, 6.5 mm thick on axis).
- Heating: deposition from the absorption coefficient β (`thermal.absorption_model: coefficient`). Brain β is specified as a fraction of attenuation and converted with the attenuation it is combined with (`itrusst_resolve_absorption.m`); skull β is a coefficient capped at the attenuation.
- Protocol (`config_itrusst_protocol.yaml`): 20 ms pulses every 200 ms for 80 s, then 80 s cooling.
- CPU (`cpp_cpu`, single precision) for calibration and Analyses 1–2; GPU (`cpp_gpu`) for Analysis 3.

## Rerunning

After editing `data/tissue_parameters.csv` or `data/tissue_setups.csv`:

```matlab
addpath('code/functions');
itrusst_backup_outputs(pwd, '<yymmdd>');  % move previous outputs aside (several runs never overwrite)
cd code
itrusst_rerun_all('prepare')    % medium configs, thickness phantoms, range/setup figures
itrusst_rerun_all('calibrate')  % free-water calibration -> protocol config
itrusst_rerun_all('simulate')   % all simulations
itrusst_rerun_all('plot')       % diagnostic figures, once all jobs have finished
itrusst_rerun_all('plot-paper') % paper figure panels (figures/paper/)
```

Without a scheduler the simulations run one after another in MATLAB (`setenv('ITRUSST_PLATFORM', 'matlab')` forces this on a cluster node). A CPU run with heating took 70–105 min on the cluster nodes; without an NVIDIA GPU, set `code_type = 'cpp_cpu'` in `PRESTUS_pipeline_01_benchmark1_sim.m`.
