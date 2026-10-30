# **ITRUSST benchmark simulations**

The goal of these analyses is to explore the impact of tissue medium parameter assumptions on relevant outcomes of acoustic intensity and heating.
This simulation uses PRESTUS, a wrapper for k-wave, to perform acoustic and heating simulations.
Different sequences and parameter setups are specified either via dedicated config files, or within the MATLAB call scripts.
Analyses are run on a GPU-supported HPC cluster hosted at the Donders Institute.

For an overview of the goals and ongoing results (*access restricted*), see [here](https://docs.google.com/document/d/16yBCTZDG966979RLcoQdo9XqpzFn1oSb-pxK2KgR4dM/edit?pli=1&tab=t.0).

The analyses are described under [Methods](#methods); script families map onto them as follows:

| analysis | scripts |
|---|---|
| 1: tool intercomparison (tissue property configurations) | `PRESTUS_pipeline_*` |
| 2: one-at-a-time parameter variation | `PRESTUS_paramvar_00-04` |
| 3: skull thickness | `PRESTUS_thickness_*` |
| 4: interactions between tissue properties (exploratory) | `code/explore/PRESTUS_paramvar_06-10`, `_13` (kept locally, not tracked) |

Shared functions live in `code/functions/` (see [code/readme.md](code/readme.md)); side analyses (`code/explore/`) and numerical checks (`code/debug/`) are kept locally and not tracked. PRESTUS is included as submodule `tools/PRESTUS` (Donders-Institute/PRESTUS); clone with `git clone --recurse-submodules` or run `git submodule update --init --recursive` after cloning.

All simulations use digital phantoms: tissue-label volumes generated in code (`code/createPhantom.m`, `code/PRESTUS_thickness_01_createPhantoms.m`), not physical phantoms and not imaging data of participants. Subject IDs (`sub-XXX`) label phantom geometries. The phantoms follow the PRESTUS examples, see [here](https://github.com/Donders-Institute/PRESTUS/blob/development/examples/createPhantom.m). The code uses PRESTUS' ```phantom``` implementation (see [here](https://github.com/Donders-Institute/PRESTUS/blob/development/documentation/doc_medium.md)). 

## Rerunning after a change of tissue parameters

Tissue properties live in `data/tissue_parameters.csv` (bounds) and `data/tissue_setups.csv` (which bound each setup uses: heat_max, heat_min, trans_max, trans_min); the `config_medium_<setup>.yaml` files are generated from them. After editing the CSVs:

```matlab
addpath('code/functions');
itrusst_backup_outputs(pwd, '<yymmdd>');  % data/tussim, data/calibration -> *_<yymmdd>, figures/* -> figures/backup/figures_<yymmdd>/
cd code
itrusst_rerun_all('prepare')    % local: medium configs, thickness phantoms, range/setup figures
itrusst_rerun_all('calibrate')  % blocking: free-water calibration -> protocol configs
itrusst_rerun_all('simulate')   % submits all simulations (runs them in turn without a scheduler)
itrusst_rerun_all('plot')       % once all jobs have finished
itrusst_rerun_all('plot-paper') % paper figure panels only (figures/paper/)
```

The backup matters: several simulations use `overwrite_files: never` and would otherwise skip conditions with existing outputs.

Without a cluster, the stages run in MATLAB: `setup_itrusst_paths` submits to SLURM/qsub when `sbatch`/`qsub` is on the PATH and otherwise sets the PRESTUS platform to `matlab` (`setenv('ITRUSST_PLATFORM', 'matlab')` forces this on a cluster node). Local runs are slow: the sweeps run one simulation after another, and a CPU run with heating took 70–105 min on the cluster nodes; and `PRESTUS_pipeline_01_benchmark1_sim.m` sets `code_type = 'cpp_gpu'` (set `'cpp_cpu'` without an NVIDIA GPU). `user = 'julkos'` in the scripts is not used: the digital phantoms need no SimNIBS installation. Side analyses (`code/explore/`) and numerical checks (`code/debug/`) are not part of the rerun and not tracked in this repository.

## Methods

### Overview

Transcranial focused ultrasound (TUS) simulations were conducted using PRESTUS (PREprocessing & Simulations for Transcranial Ultrasound Stimulation; v05), an open-source MATLAB toolbox that integrates k-Wave pseudospectral acoustic simulation with Pennes bioheat equation thermal modelling. All analyses were performed within the PRESTUS_itrusst analysis environment, which extends PRESTUS with a structured benchmarking framework to evaluate the sensitivity of simulation outputs to assumptions about tissue acoustic and thermal properties. Three main analyses were performed: (1) a pipeline comparison across eight tissue property configurations, (2) a systematic one-at-a-time sensitivity analysis of individual tissue parameters, and (3) an evaluation of the effect of skull thickness on transcranial transmission. A fourth, exploratory analysis tests interactions between tissue properties.

---

### Simulation Framework

#### Digital phantom geometry

Two-dimensional axisymmetric digital phantoms (tissue-label volumes generated in code) were constructed at 0.5 mm isotropic resolution and stored as NIfTI volumes in the SimNIBS m2m-folder layout that PRESTUS reads (`final_tissues.nii.gz`); no SimNIBS segmentation was run. The active simulation domain spanned 70 × 120 mm (140 × 240 voxels laterally and axially), with a 20 mm perfectly-matched layer (PML; 40 voxels at 0.5 mm resolution) appended externally by k-Wave to absorb outgoing waves.

Two standard benchmarks were used. **Benchmark 1** modelled a water–soft tissue interface, with soft tissue (label 1) beginning at 30 mm depth from the proximal domain edge. **Benchmark 2** modelled a transcranial geometry comprising a water column, a spherical skull cap (outer and inner radius of curvature 75 and 68.5 mm about a common centre on the brain side; outer apex at 30 mm depth, 6.5 mm thick on the beam axis), and brain tissue filling the interior. Both phantoms were treated as 2D axisymmetric (cylindrical symmetry) to reduce computation time. For the pipeline comparison (Analysis 1), three-dimensional phantoms of both benchmarks (sub-101, sub-102; 140 × 140 × 240 voxels at 0.5 mm) were used. Their acoustic simulations were also axisymmetric (140 × 240 half-plane grid, C++ GPU code); the pressure field was expanded radially to 3D, and the thermal simulation ran on the 140 × 140 × 240 grid.

#### Transducer

A single-element annular transducer was modelled (outer diameter 64 mm, radius of curvature 64 mm, inner diameter 0 mm) operating at a nominal focal depth of 64 mm. The transducer bowl was placed at the proximal edge of the active domain (voxel coordinate [70, 1] at 0.5 mm resolution); the geometric focus is computed from the bowl focal distance as voxel [70, 129] (64 mm axial depth).

#### Transducer calibration

The source amplitude was calibrated once per protocol by a free-water acoustic simulation (`code/itrusst_calibrate_all.m` → `itrusst_calibrate`); calibration is transducer-specific and independent of the phantom. The peak spatial-peak pulse-average intensity (Isppa) along the axial profile was extracted from the simulated pressure field as:

> Isppa = p²_max / (2ρc) × 10⁻⁴  [W/cm²]

Because intensity scales with the square of the source pressure in the linear simulation, the source pressure amplitude (`elem_p_amp`) was multiplied by √(60 W/cm² / Isppa) to reach the target intensity of 60 W/cm². The transducer has a single element, so its phase has no effect and was kept from the config. The calibrated amplitude and phase are written to `data/calibration/PRESTUS_virtual_parameters/` (`save_optimized_values`), then copied into the `annular` block of `data/configs/config_itrusst_protocol.yaml` (`itrusst_apply_calibration`) and used by all subsequent tissue simulations.

#### Acoustic simulation

Acoustic pressure fields were computed using k-Wave (v1.4, bundled with PRESTUS), a pseudospectral time-domain solver implementing the linearised pressure–velocity equations in heterogeneous media. Simulations used a Courant-Friedrichs-Lewy (CFL) stability number of 0.15 and were run with the k-Wave C++ binaries: CPU (`cpp_cpu`) for calibration, the reference-medium baseline, the one-at-a-time sweep and the thickness sweeps, GPU (`cpp_gpu`) for the pipeline comparison. The output of each acoustic simulation was the spatial maximum pressure field p_max(x) across all time steps, from which the following metrics were derived per tissue compartment:

- **Isppa** [W/cm²]: spatial peak pulse-average intensity = p²_max / (2ρc)
- **Ipa_target** [W/cm²]: mean intensity averaged over a circular target region
- **Psptp** [Pa]: spatial peak temporal peak pressure
- **Mechanical Index (MI)**: p_neg / (√f × 10⁶), a dimensionless cavitation risk metric
- **Real focal distance** [mm]: depth of the peak intensity from the transducer exit plane

#### Thermal simulation

Thermal simulations were performed using k-Wave's `kWaveDiffusion` solver, which implements the Pennes bioheat transfer equation:

> ρ c_p (∂T/∂t) = ∇(k∇T) − ρ_b c_b ω (T − T_blood) + Q

where ρ is density, c_p is specific heat capacity, k is thermal conductivity, ω is the blood perfusion coefficient [s⁻¹], T_blood = 37°C, and Q is the acoustic heat deposition rate:

> Q = α_abs × p²_max / (ρ × c)  [W/m³]

where c is the sound speed and α_abs is the absorption coefficient [Np/m] at the source frequency, passed to PRESTUS per tissue as `absorption_coeff` [dB/(MHz^y·cm)] with the tissue's attenuation power law y (`thermal.absorption_model: coefficient`, set in `config_itrusst_protocol.yaml` and the medium configs). Brain absorption is specified as a fraction of attenuation (`absorption_fraction`; 29 / 77 / 100 %), skull absorption as a coefficient. `itrusst_resolve_absorption.m` converts the brain fraction with the attenuation it is combined with (α_abs = fraction × α), because PRESTUS applies one absorption model per run.

The sonication timing followed a pulse-train paradigm: pulse duration (PD) 20 ms, pulse repetition interval (PRI) 200 ms, pulse train duration (PTD) 200 ms, pulse train repetition interval (PTRI) 200 ms, and total pulse-train repetition duration (PTRD) 80 s, followed by an 80 s post-PTRI cooldown. Thermal modelling time steps were 20 ms during the pulse and 1 s post-pulse. Thermal outputs included maximum temperature rise (maxT, riseT), end-of-sonication temperature (endT), and thermal dose expressed as Cumulative Equivalent Minutes at 43°C (CEM43) computed both by the k-Wave native method and per ISO 13 standards. All metrics were extracted per tissue compartment (brain, skull) and written to per-subject CSV summary tables.

---

### Analysis 1: Pipeline Benchmark — Tissue Property Configuration Comparison

To quantify how assumptions about tissue acoustic and thermal properties propagate into simulation outputs, eight property configurations were evaluated: four literature setups, the literature reference and the defaults of three TUS simulation tools.

| Configuration | Rationale |
|---|---|
| **Temp-H** (`heat_max`) | Per property, the literature bound that maximises predicted brain heating (the kPlan defaults exceed Temp-H in skull absorption) |
| **Temp-L** (`heat_min`) | Per property, the bound that minimises predicted brain heating |
| **Trans-H** (`trans_max`) | Per acoustic property (sound speed, density, attenuation), the bound that maximises brain intensity; acoustic simulation only |
| **Trans-L** (`trans_min`) | Per acoustic property, the bound that minimises brain intensity; acoustic simulation only |
| **reference** | Literature reference values |
| **PRESTUS** | PRESTUS toolbox defaults |
| **BabelBrain** | BabelBrain defaults |
| **kPlan** | kPlan defaults (perfusion not reported; simulated as 0) |

The literature values (bottom/reference/top per tissue and property) are kept in `data/tissue_parameters.csv` (source table: `data/tissue_parameters_source.md`), and the bound each setup uses per property in `data/tissue_setups.csv` (one column per setup; the Trans setups leave the thermal properties unspecified, which are then held at the reference). The bound directions follow the one-at-a-time effects (Analysis 2). `config_medium_heat_max/heat_min/trans_max/trans_min/reference.yaml` are generated from these files (`code/functions/itrusst_write_medium_configs.m`); the tool defaults are hand-written `config_medium_<tool>.yaml` files. Water is identical across configurations (c = 1482 m/s, ρ = 994 kg/m³, α = 0.00055 dB/cm at 500 kHz with y = 2).

For each configuration, simulations were run on the 3D phantoms of Benchmark 1 (soft tissue only, sub-101) and Benchmark 2 (skull + brain, sub-102) at a target intensity of 60 W/cm² (`PRESTUS_pipeline_01_benchmark1_sim.m`; axisymmetric acoustics, thermal simulation on the radially expanded 3D field; Trans-H/L without the thermal simulation). `PRESTUS_pipeline_02_benchmark1_plot.m` compares the temperature rise at the focus and the maximum brain and skull temperature rise (riseT_brain, riseT_skull) across Temp-H, Temp-L and the tools, and the maximum brain intensity (Isppa_brain) across Trans-H, Trans-L and the tools.

Attenuation and absorption are specified at 500 kHz and converted to PRESTUS' 1 MHz coefficients α₀ = α(500 kHz) / 0.5^y and α_abs,0 = α_abs(500 kHz) / 0.5^y (dB/(MHz^y·cm); y = 1 for brain and skull). Both are listed below. PRESTUS uses its v05 default absorption coefficients (`config/config_default.yaml`; brain = attenuation, a conservative and provisional default; skull 3.78 dB/(MHz·cm) after Pinton et al., 2012) and its default brain power law y = 1.2. BabelBrain and kPlan specify absorption as a fraction of attenuation; their absorption coefficients are that fraction × the tool's own α₀ (BabelBrain: brain 0.85, skull 0.163; kPlan: 1.0, 1.0). All tool configurations are simulated with `thermal.absorption_model: coefficient` like the literature setups.

**Brain properties by configuration** (Trans-H/L: thermal properties and absorption not used, –):

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | α_abs @ 500 kHz (dB/cm) | α_abs,0 (dB/(MHz·cm)) |
|---|---|---|---|---|---|---|---|---|---|
| Temp-H | 1573 | 1040 | 0.493 | 0.986 | 0.48 | 3525 | 196 | 0.493 | 0.986 |
| Temp-L | 1520 | 1050 | 0.154 | 0.308 | 0.55 | 3718 | 809 | 0.045 | 0.089 |
| Trans-H | 1520 | 1050 | 0.154 | 0.308 | – | – | – | – | – |
| Trans-L | 1573 | 1040 | 0.493 | 0.986 | – | – | – | – | – |
| reference | 1560 | 1040 | 0.30 | 0.60 | 0.51 | 3630 | 559 | 0.231 | 0.462 |
| PRESTUS | 1546 | 1046 | 0.257 | 0.59 (y = 1.2) | 0.51 | 3630 | 559 | 0.257 | 0.59 (y = 1.2) |
| BabelBrain | 1562 | 1041 | 0.24 | 0.48 | 0.51 | 3630 | 559 | 0.204 | 0.408 |
| kPlan | 1550 | 1045 | 0.26 | 0.52 | 0.51 | 3600 | 0 | 0.26 | 0.52 |

**Skull properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | α_abs @ 500 kHz (dB/cm) | α_abs,0 (dB/(MHz·cm)) |
|---|---|---|---|---|---|---|---|---|---|
| Temp-H | 2871 | 1080 | 4.02 | 8.04 | 0.30 | 826 | 0 | 4.02 | 8.04 |
| Temp-L | 2190 | 1991 | 11.7 | 23.4 | 0.36 | 2524 | 30 | 1.05 | 2.10 |
| Trans-H | 2190 | 1080 | 4.02 | 8.04 | – | – | – | – | – |
| Trans-L | 2871 | 1991 | 11.7 | 23.4 | – | – | – | – | – |
| reference | 2550 | 1775 | 6.0 | 12.0 | 0.316 | 1793 | 20 | 1.40 | 2.80 |
| PRESTUS | 2800 | 1850 | 6.65 | 13.3 | 0.32 | 1793 | 20 | 1.89 | 3.78 |
| BabelBrain | 2800 | 1850 | 6.65 | 13.3 | 0.315 | 1793 | 20 | 1.08 | 2.17 |
| kPlan | 2627 | 1850 | 6.65 | 13.3 | 0.32 | 1300 | 0 | 6.65 | 13.3 |

In Temp-H the skull absorbs all attenuated energy: its top absorption coefficient (8.04 dB/(MHz·cm)) equals the bottom attenuation that Temp-H uses.

The skull reference sound speed, density and attenuation are the mean of the cortical and trabecular bone values of Aubry et al. (2022) rather than literature-table values (`data/tissue_parameters_source.md`).

---

### Analysis 2: Sensitivity Analysis — One-at-a-Time Parameter Perturbation

To isolate the contribution of individual tissue parameters to simulation outputs, a one-at-a-time (OAT) sensitivity analysis was performed. Each parameter was independently set to its literature-reported lower or upper bound while all other parameters were held at their reference values. Seven acoustic and thermal parameters were perturbed for each of two tissue types (brain, skull):

| Parameter | Symbol | Unit |
|---|---|---|
| Sound speed | c | m/s |
| Density | ρ | kg/m³ |
| Acoustic attenuation coefficient | α₀ | dB/(MHz^y·cm) |
| Thermal conductivity | k | W/m/K |
| Specific heat capacity | c_p | J/kg/K |
| Blood perfusion | ω | mL/min/kg |
| Absorption coefficient | α_abs,0 | dB/(MHz^y·cm) |

Parameter bounds were drawn from a literature-informed tissue property table (`tissue_parameters.csv`). Brain absorption is a fraction of the attenuation it is combined with (`itrusst_resolve_absorption.m`): with attenuation at reference, its bounds are 0.087 and 0.30 dB/cm at 500 kHz (29 and 100 %), and the brain attenuation bounds scale the reference absorption (77 %) to 0.119 and 0.380 dB/cm. Skull absorption is a coefficient, capped at the attenuation. Brain parameters were perturbed for Benchmark 1 (sub-001); both brain and skull parameters were perturbed for Benchmark 2 (sub-002) and for the 2.5 mm thickness phantom (sub-352, the most transmissive skull thickness; see Analysis 3). For each perturbation, a separate simulation was run with a unique output filename affix encoding the tissue, parameter, and bound level (e.g., `_b_c_1` = brain sound speed lower bound; `_s_a0_2` = skull attenuation upper bound). One protocol variant was evaluated (`itrusst_protocol1a`: `config_itrusst_protocol.yaml`, 500 kHz, 60 W/cm²; outputs in `data/tussim/itrusst_protocol1a/`).

Output metrics extracted for each perturbation included Ipa_target, maxT, Isppa_brain, riseT_brain, and real focal distance. All simulations used C++ CPU code in single-precision floating point.

---

### Analysis 3: Skull Thickness Sensitivity

To characterise the relationship between skull thickness and transcranial acoustic transmission, 23 digital skull phantoms were constructed with uniform cortical bone thickness from 1 to 12 mm in 0.5 mm steps (one voxel per step at 0.5 mm resolution; subjects sub-301 to sub-312 for 1–12 mm and sub-351 to sub-361 for 1.5–11.5 mm, see `code/functions/itrusst_thickness_subjects.m`). For each phantom, the spherical skull cap was placed at a fixed outer radius of curvature of 75 mm (outer apex at 30 mm depth), with the inner radius adjusted to achieve the target thickness. All tissue properties were set to their reference values. Simulations were acoustic only (`PRESTUS_thickness_02_sim.m`) and used C++ CPU code at 0.5 mm resolution. The primary outcomes were Ipa_target (at the fixed target) and Isppa_brain (maximum in brain) as a function of skull thickness, quantifying acoustic insertion loss attributable to bone. Transmission is not monotonic in thickness: Isppa_brain oscillates with a period of about half a wavelength in bone (2.55 mm at 2550 m/s) on top of the attenuation decay. With the current tissue parameters the most transmissive thickness is 2.5 mm (sub-352, Isppa_brain 21.9 W/cm²), which therefore serves as the maximal-transmission benchmark in the one-at-a-time analysis. In the previous parameter set (`data/tussim_261008`) the peak lay at 3 mm (sub-303, now 20.9 W/cm², 4.5 % below the peak), which was the earlier benchmark.

As a follow-up, skull density was varied across its bottom, reference and top levels on every thickness phantom (all other properties at their reference values, heating enabled, shortened 20 s thermal protocol; `PRESTUS_thickness_04_density_sim.m` / `_05_density_plot.m`), to test whether the effect of skull density depends on skull thickness.

Skull sound speed was swept at 9 evenly spaced values between its bottom and top bound plus the reference value (2190–2871 m/s) on the 3, 6.5 and 9 mm phantoms (sub-303/356/309; all other properties at their reference values, heating enabled, shortened 20 s thermal protocol; `PRESTUS_thickness_06_soundspeed_sim.m` / `_07_soundspeed_plot.m`). Motivation: on sub-002 the reference value transmits less intensity to the brain than either bound, in PRESTUS and BabelBrain, consistent with a skull-layer resonance (transmission minimum at an odd multiple of λ/4 in bone, maximum at multiples of λ/2). A plane-wave layer estimate (normal incidence, reference density and attenuation, 500 kHz) predicts transmission rising monotonically with sound speed at 3 mm, a minimum near 2680 m/s at 6.5 mm, and a maximum near 2230 m/s at 9 mm; the plot overlays this estimate on brain I_SPPA.

---

### Analysis 4: Interactions Between Tissue Properties (exploratory)

*Status: exploratory, not part of the main analysis; the scripts (`code/explore/`) are kept locally and not included in this repository.*

The pipeline-comparison configurations (Analysis 1) and the one-at-a-time sensitivity sweep (Analysis 2) each probe the tissue-property space in a limited way: the former evaluates only two literature-bound "constellations" (all seven properties per tissue simultaneously at their Temp-H or Temp-L bound), and the latter perturbs one property at a time from the reference values, holding all others fixed. Neither directly establishes whether combining every property's individually assigned bound actually produces the expected combined (maximal/minimal) outcome, since one-at-a-time effects need not sum linearly when parameters are varied jointly.

A preliminary check compared, for each combined "constellation", the actual simulated outcome against the sum of the seven (Benchmark 1) or fourteen (Benchmark 2) individual one-at-a-time effects from Analysis 2 relative to the reference configuration. Every constellation moved each outcome (riseT_brain, riseT_skull, Isppa_brain, Ipa_target) in the direction predicted by the additive sum, but the magnitude deviated substantially: heating outcomes were super-additive (actual change 1.6–2.0× the additive prediction, i.e. combining bounds heats more than one-at-a-time effects would suggest), while intensity outcomes were close to additive with no skull present but sub-additive with skull present (actual change 0.65–0.79× predicted), consistent with attenuation acting multiplicatively (exponentially) rather than additively along the propagation path. This motivated testing interactions between properties directly rather than relying on the additive approximation.

Both designs below work on Benchmark 2 (sub-002) with the 14 brain and skull properties at their literature bottom/top bounds. Every evaluation is a full, independent acoustic + thermal simulation (no acoustic cache is reused; see the script headers for why), with the thermal protocol shortened to 20 s stimulation and no post-stimulation cool-off, since only riseT_brain/riseT_skull/Isppa_brain/Ipa_target are needed. This override is scoped to these scripts and is not applied to `config_itrusst_protocol.yaml`. Each evaluation is submitted as its own job, without dependencies between jobs; NIfTI creation, reports and matrix saving are disabled.

**Skull-density crossover.** `PRESTUS_paramvar_08_constellation_factorial_sim.m` crosses skull density (bottom/top) with each of the other 13 properties (bottom/top) in a 2×2 grid, all other properties at their reference values (13 pairs × 4 = 52 evaluations). The interaction per pair is [outcome(ρ_skull top, factor top) − outcome(ρ_skull bottom, factor top)] − [outcome(ρ_skull top, factor bottom) − outcome(ρ_skull bottom, factor bottom)]; it is analysed by `PRESTUS_paramvar_13_density_crossover_plot.m`.

**Morris elementary effects (interaction screening).** `PRESTUS_paramvar_09_morris_interactions_sim.m` screens all properties at once via a Morris (1991) elementary-effects screening design on the same 14-property, 2-level (bottom/top) space: r independent trajectories, each a "staircase" walk of 15 points (1 random baseline + all 14 properties flipped one at a time, in a random order grouped acoustic-factors-first, kept from an earlier caching scheme so that the design (rng seed 42) is unchanged), from which each property's mean absolute elementary effect (μ\*, overall importance) and standard deviation of elementary effects (σ, a direct per-property measure of how much its effect depends on the other 13 properties' current values, i.e. interaction) are computed. At r=10 (default), the design has 150 points (shared points are submitted once), each a full, independent acoustic + thermal simulation with the shortened 20 s thermal protocol, submitted as one job per evaluation without dependencies; cost (≈30–50 CPU job-hours at 6–29 min per evaluation) scales linearly with r rather than exponentially with the number of properties. No acoustic cache is reused: an earlier version reused each trajectory's acoustic field for its 8 thermal-only steps, but on the axisymmetric grid all cache-reusing jobs failed (cached 3D-expanded pressure field vs. 2D medium), leaving that run (`data/tussim_261001/itrusst_protocol1a_morris`, previous tissue ranges) with 47/150 results and no thermal elementary effects. The trajectory design (including every step's output_affix) is written to `data/tussim/itrusst_protocol1a_morris/morris_design.csv`; `PRESTUS_paramvar_10_morris_analysis_plot.m` joins it with the per-evaluation result CSVs and computes μ, μ\*, σ and n per property, tolerating a partially completed design.

---

### Outputs and Visualisation

Simulation results for each subject and configuration were written to CSV tables (`sub-XXX_phantom[_affix].csv`) located in per-setup output directories. These tables contain the full set of acoustic and thermal metrics enumerated above. For the pipeline comparison, bar plots per benchmark show the temperature rise at the focus, riseT_brain and riseT_skull for Temp-H, Temp-L and the tools, and Isppa_brain for Trans-H, Trans-L and the tools. For the sensitivity analysis, metrics were plotted as min–max ranges with individual data points overlaid, colour-coded by tissue type (brain: grey; skull: red). For the skull thickness analysis, Ipa_target and Isppa_brain were plotted as a continuous function of skull thickness in millimetres.
