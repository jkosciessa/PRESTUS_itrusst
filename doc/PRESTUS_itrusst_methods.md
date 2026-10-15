# Methods

## Overview

Transcranial focused ultrasound (TUS) simulations were conducted using PRESTUS (Pipeline for Running E-field and TUS Simulations; v05), an open-source MATLAB toolbox that integrates k-Wave pseudospectral acoustic simulation with Pennes bioheat equation thermal modelling. All analyses were performed within the PRESTUS_itrusst analysis environment, which extends PRESTUS with a structured benchmarking framework to evaluate the sensitivity of simulation outputs to assumptions about tissue acoustic and thermal properties. Three analyses were performed: (1) a pipeline comparison across six tissue property configurations, (2) a systematic one-at-a-time sensitivity analysis of individual tissue parameters, and (3) an evaluation of the effect of skull thickness on transcranial transmission.

---

## Simulation Framework

### Phantom geometry

Two-dimensional axisymmetric phantoms were constructed at 0.5 mm isotropic resolution and stored as NIfTI volumes in SimNIBS m2m-folder format (`final_tissues.nii.gz`). The active simulation domain spanned 70 × 120 mm (140 × 240 voxels laterally and axially), with a 20 mm perfectly-matched layer (PML; 40 voxels at 0.5 mm resolution) appended externally by k-Wave to absorb outgoing waves.

Two standard benchmarks were used. **Benchmark 1** modelled a water–soft tissue interface, with soft tissue (label 1) beginning at 30 mm depth from the proximal domain edge. **Benchmark 2** modelled a transcranial geometry comprising a water column, a hemispherical skull annulus (inner radius 68.5 mm, outer radius 75 mm from the geometric centre), and brain tissue filling the interior. Both phantoms were treated as 2D axisymmetric (cylindrical symmetry) to reduce computation time. For the pipeline comparison (Analysis 1), three-dimensional counterparts of both benchmarks (sub-101, sub-102; 140 × 140 × 240 voxels at 0.5 mm) were simulated with C++ GPU code.

### Transducer

A single-element annular transducer was modelled (outer diameter 64 mm, radius of curvature 64 mm, inner diameter 0 mm) operating at a nominal focal depth of 64 mm. The transducer bowl was placed at the proximal edge of the active domain (voxel coordinate [70, 1] at 0.5 mm resolution); the geometric focus is computed from the bowl focal distance as voxel [70, 129] (64 mm axial depth).

### Transducer calibration

The source amplitude was calibrated once per protocol by a free-water acoustic simulation (`code/itrusst_calibrate_all.m` → `itrusst_calibrate`); calibration is transducer-specific and independent of the phantom. The peak spatial-peak pulse-average intensity (Isppa) along the axial profile was extracted from the simulated pressure field as:

> Isppa = p²_max / (2ρc) × 10⁻⁴  [W/cm²]

A linear correction factor was then applied to scale the source pressure amplitude (`elem_p_amp`) to achieve the target intensity of 60 W/cm². The calibrated amplitude and phase are copied into the `annular` block of `data/configs/config_itrusst_protocol.yaml` (`itrusst_apply_calibration`) and used by all subsequent tissue simulations.

### Acoustic simulation

Acoustic pressure fields were computed using k-Wave (v1.4.1), a pseudospectral time-domain solver implementing the linearised pressure–velocity equations in heterogeneous media. Simulations used a Courant-Friedrichs-Lewy (CFL) stability number of 0.15 and were run with the k-Wave C++ binaries: CPU (`cpp_cpu`) for calibration, the one-at-a-time sweep and the thickness sweeps, GPU (`cpp_gpu`) for the reference-medium baseline and the 3D pipeline comparison. The output of each acoustic simulation was the spatial maximum pressure field p_max(x) across all time steps, from which the following metrics were derived per tissue compartment:

- **Isppa** [W/cm²]: spatial peak pulse-average intensity = p²_max / (2ρc)
- **Ipa_target** [W/cm²]: mean intensity averaged over a circular target region
- **Psptp** [Pa]: spatial peak temporal peak pressure
- **Mechanical Index (MI)**: p_neg / (√f × 10⁶), a dimensionless cavitation risk metric
- **Real focal distance** [mm]: depth of the peak intensity from the transducer exit plane

### Thermal simulation

Thermal simulations were performed using k-Wave's `kWaveDiffusion` solver, which implements the Pennes bioheat transfer equation:

> ρc (∂T/∂t) = ∇(k∇T) − ρ_b c_b ω (T − T_blood) + Q

where ρ is density, c is specific heat capacity, k is thermal conductivity, ω is the blood perfusion coefficient [s⁻¹], T_blood = 37°C, and Q is the acoustic heat deposition rate:

> Q = α × p²_max / (ρ × c) × absorption_fraction  [W/m³]

The sonication timing followed a pulse-train paradigm: pulse duration (PD) 20 ms, pulse repetition interval (PRI) 200 ms, pulse train duration (PTD) 200 ms, pulse train repetition interval (PTRI) 200 ms, and total pulse-train repetition duration (PTRD) 80 s, followed by an 80 s post-PTRI cooldown. Thermal modelling time steps were 20 ms during the pulse and 1 s post-pulse. Thermal outputs included maximum temperature rise (maxT, riseT), end-of-sonication temperature (endT), and thermal dose expressed as Cumulative Equivalent Minutes at 43°C (CEM43) computed both by the k-Wave native method and per ISO 13 standards. All metrics were extracted per tissue compartment (brain, skull, skin) and written to per-subject CSV summary tables.

---

## Analysis 1: Pipeline Benchmark — Tissue Property Configuration Comparison

To quantify how assumptions about tissue acoustic and thermal properties propagate into simulation outputs, six property configurations were evaluated: three literature setups and the defaults of three TUS simulation tools.

| Configuration | Rationale |
|---|---|
| **setup 1** | Per property, the literature bound that maximises predicted brain heating |
| **setup 2** | Per property, the opposite bound (minimises predicted heating) |
| **reference** | Literature reference values |
| **PRESTUS** | PRESTUS toolbox defaults |
| **BabelBrain** | BabelBrain defaults |
| **kPlan** | kPlan defaults (perfusion not reported; simulated as 0) |

The literature values (bottom/reference/top per tissue and property) are kept in `data/tissue_parameters.csv` (source table: `data/tissue_parameters_source.md`), and the bound used by setup 1 and setup 2 per property in `data/tissue_setups.csv`. The bound directions follow the one-at-a-time effect on brain heating (Analysis 2). `config_medium_setup1/reference/setup2.yaml` are generated from these files (`code/functions/itrusst_write_medium_configs.m`); the tool defaults are hand-written `config_medium_<tool>.yaml` files. Water is identical across configurations (c = 1482 m/s, ρ = 994 kg/m³, α = 0.00055 dB/cm at 500 kHz with y = 2).

For each configuration, acoustic and thermal simulations were run on the 3D versions of Benchmark 1 (soft tissue only, sub-101) and Benchmark 2 (skull + brain, sub-102) at a target intensity of 60 W/cm² (`PRESTUS_pipeline_01_benchmark1_sim.m`). Simulation outputs (maxT, Ipa_target) were extracted from per-subject CSV tables and compared across configurations.

Attenuation is specified at 500 kHz and converted to PRESTUS' 1 MHz coefficient α₀ = α(500 kHz) / 0.5^y (dB/(MHz^y·cm); y = 1 for brain and skull). Both are listed below.

**Brain properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | Absorption fraction |
|---|---|---|---|---|---|---|---|---|
| setup 1 | 1520 | 1040 | 0.493 | 0.986 | 0.48 | 3525 | 196 | 1.0 |
| setup 2 | 1568 | 1050 | 0.236 | 0.472 | 0.55 | 3718 | 809 | 0.73 |
| reference | 1560 | 1040 | 0.30 | 0.60 | 0.51 | 3630 | 559 | 0.865 |
| PRESTUS | 1546 | 1046 | 0.26 | 0.52 | 0.51 | 3630 | 559 | 1.0 |
| BabelBrain | 1562 | 1041 | 0.24 | 0.48 | 0.51 | 3630 | 559 | 0.85 |
| kPlan | 1550 | 1045 | 0.26 | 0.52 | 0.51 | 3600 | 0 | 1.0 |

**Skull properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | Absorption fraction |
|---|---|---|---|---|---|---|---|---|
| setup 1 | 2871 | 1080 | 11.7 | 23.4 | 0.30 | 826 | 0 | 1.0 |
| setup 2 | 2190 | 1991 | 3.22 | 6.44 | 0.36 | 2524 | 30 | 0.16 |
| reference | 2550 | 1775 | 6.0 | 12.0 | 0.315 | 1793 | 20 | 0.58 |
| PRESTUS | 2800 | 1850 | 6.65 | 13.3 | 0.32 | 1793 | 20 | 0.28 |
| BabelBrain | 2800 | 1850 | 6.65 | 13.3 | 0.315 | 1793 | 20 | 0.163 |
| kPlan | 2800 | 1850 | 6.65 | 13.3 | 0.32 | 1300 | 0 | 1.0 |

The skull reference sound speed, density and attenuation are the mean of the cortical and trabecular bone values of Aubry et al. (2022) rather than literature-table values (`data/tissue_parameters_source.md`).

---

## Analysis 2: Sensitivity Analysis — One-at-a-Time Parameter Perturbation

To isolate the contribution of individual tissue parameters to simulation outputs, a one-at-a-time (OAT) sensitivity analysis was performed. Each parameter was independently set to its literature-reported lower or upper bound while all other parameters were held at their reference values. Seven acoustic and thermal parameters were perturbed for each of two tissue types (brain, skull):

| Parameter | Symbol | Unit |
|---|---|---|
| Sound speed | c | m/s |
| Density | ρ | kg/m³ |
| Acoustic attenuation coefficient | α₀ | dB/(MHz^y·cm) |
| Thermal conductivity | k | W/m/K |
| Specific heat capacity | c_p | J/kg/K |
| Blood perfusion | ω | mL/min/kg |
| Absorption fraction | — | — |

Parameter bounds were drawn from a literature-informed tissue property table (`tissue_parameters.csv`). Brain parameters were perturbed for Benchmark 1 (sub-001); both brain and skull parameters were perturbed for Benchmark 2 (sub-002) and for the 2.5 mm thickness phantom (sub-352, the most transmissive skull thickness; see Analysis 3). For each perturbation, a separate simulation was run with a unique output filename affix encoding the tissue, parameter, and bound level (e.g., `_b_c_1` = brain sound speed lower bound; `_s_a0_2` = skull attenuation upper bound). One protocol variant was evaluated (`itrusst_protocol1a`: `config_itrusst_protocol.yaml`, 500 kHz, 60 W/cm²; outputs in `data/tussim/itrusst_protocol1a/`).

Output metrics extracted for each perturbation included Ipa_target, maxT, Isppa_brain, riseT_brain, and real focal distance. All simulations used C++ CPU code in single-precision floating point.

---

## Analysis 3: Skull Thickness Sensitivity

To characterise the relationship between skull thickness and transcranial acoustic transmission, 23 skull phantoms were constructed with uniform cortical bone thickness from 1 to 12 mm in 0.5 mm steps (one voxel per step at 0.5 mm resolution; subjects sub-301 to sub-312 for 1–12 mm and sub-351 to sub-361 for 1.5–11.5 mm, see `code/functions/itrusst_thickness_subjects.m`). For each phantom, the spherical skull cap was placed at a fixed outer radius of curvature of 75 mm (outer apex at 30 mm depth), with the inner radius adjusted to achieve the target thickness. All tissue properties were set to their reference values. Simulations used C++ CPU code at 0.5 mm resolution. The primary outcomes were Ipa_target (at the fixed target) and Isppa_brain (maximum in brain) as a function of skull thickness, quantifying acoustic insertion loss attributable to bone. Transmission is not monotonic in thickness: Isppa_brain oscillates with a period of about half a wavelength in bone (2.55 mm at 2550 m/s) on top of the attenuation decay. With the current tissue parameters the most transmissive thickness is 2.5 mm (sub-352, Isppa_brain 21.7 W/cm²), which therefore serves as the maximal-transmission benchmark in the one-at-a-time analysis. In the previous parameter set (`data/tussim_261008`) the peak lay at 3 mm (sub-303, now 20.8 W/cm², 4 % below the peak), which was the earlier benchmark.

As a follow-up, skull density was varied across its bottom, reference and top levels on every thickness phantom (all other properties at their reference values, heating enabled, shortened 20 s thermal protocol; `PRESTUS_thickness_04_density_sim.m` / `_05_density_plot.m`), to test whether the effect of skull density depends on skull thickness.

Skull sound speed was swept at 9 evenly spaced values between its bottom and top bound plus the reference value (2190–2871 m/s) on the 3, 6.5 and 9 mm phantoms (sub-303/356/309; all other properties at their reference values, heating enabled, shortened 20 s thermal protocol; `PRESTUS_thickness_06_soundspeed_sim.m` / `_07_soundspeed_plot.m`). Motivation: on sub-002 the reference value transmits less intensity to the brain than either bound, in PRESTUS and BabelBrain, consistent with a skull-layer resonance (transmission minimum at an odd multiple of λ/4 in bone, maximum at multiples of λ/2). A plane-wave layer estimate (normal incidence, reference density and attenuation, 500 kHz) predicts transmission rising monotonically with sound speed at 3 mm, a minimum near 2680 m/s at 6.5 mm, and a maximum near 2230 m/s at 9 mm; the plot overlays this estimate on brain I_SPPA.

---

## Analysis 4: Interactions Between Tissue Properties (exploratory)

*Status: exploratory, not part of the main analysis (scripts in `code/explore/`).*

The pipeline-comparison configurations (Analysis 1) and the one-at-a-time sensitivity sweep (Analysis 2) each probe the tissue-property space in a limited way: the former evaluates only two literature-bound "constellations" (all seven properties per tissue simultaneously at their setup 1 or setup 2 bound), and the latter perturbs one property at a time from the reference values, holding all others fixed. Neither directly establishes whether combining every property's individually assigned bound actually produces the expected combined (maximal/minimal) outcome, since one-at-a-time effects need not sum linearly when parameters are varied jointly.

A preliminary check compared, for each combined "constellation", the actual simulated outcome against the sum of the seven (Benchmark 1) or fourteen (Benchmark 2) individual one-at-a-time effects from Analysis 2 relative to the reference configuration. Every constellation moved each outcome (riseT_brain, riseT_skull, Isppa_brain, Ipa_target) in the direction predicted by the additive sum, but the magnitude deviated substantially: heating outcomes were super-additive (actual change 1.6–2.0× the additive prediction, i.e. combining bounds heats more than one-at-a-time effects would suggest), while intensity outcomes were close to additive with no skull present but sub-additive with skull present (actual change 0.65–0.79× predicted), consistent with attenuation acting multiplicatively (exponentially) rather than additively along the propagation path. This motivated testing interactions between properties directly rather than relying on the additive approximation.

Both designs below work on Benchmark 2 (sub-002) with the 14 brain and skull properties at their literature bottom/top bounds. Every evaluation is a full, independent acoustic + thermal simulation (no acoustic cache is reused; see the script headers for why), with the thermal protocol shortened to 20 s stimulation and no post-stimulation cool-off, since only riseT_brain/riseT_skull/Isppa_brain/Ipa_target are needed. This override is scoped to these scripts and is not applied to `config_itrusst_protocol.yaml`. Each evaluation is submitted as its own job, without dependencies between jobs; NIfTI creation, reports and matrix saving are disabled.

**Skull-density crossover.** `PRESTUS_paramvar_08_constellation_factorial_sim.m` crosses skull density (bottom/top) with each of the other 13 properties (bottom/top) in a 2×2 grid, all other properties at their reference values (13 pairs × 4 = 52 evaluations). The interaction per pair is [outcome(ρ_skull top, factor top) − outcome(ρ_skull bottom, factor top)] − [outcome(ρ_skull top, factor bottom) − outcome(ρ_skull bottom, factor bottom)]; it is analysed by `PRESTUS_paramvar_13_density_crossover_plot.m`.

**Morris elementary effects (interaction screening).** `PRESTUS_paramvar_09_morris_interactions_sim.m` screens all properties at once via a Morris (1991) elementary-effects screening design on the same 14-property, 2-level (bottom/top) space: r independent trajectories, each a "staircase" walk of 15 points (1 random baseline + all 14 properties flipped one at a time, in a random order grouped acoustic-factors-first, kept from an earlier caching scheme so that the design (rng seed 42) is unchanged), from which each property's mean absolute elementary effect (μ\*, overall importance) and standard deviation of elementary effects (σ, a direct per-property measure of how much its effect depends on the other 13 properties' current values, i.e. interaction) are computed. At r=10 (default), the design has 150 points (shared points are submitted once), each a full, independent acoustic + thermal simulation with the shortened 20 s thermal protocol, submitted as one job per evaluation without dependencies; cost (≈30–50 CPU job-hours at 6–29 min per evaluation) scales linearly with r rather than exponentially with the number of properties. No acoustic cache is reused: an earlier version reused each trajectory's acoustic field for its 8 thermal-only steps, but on the axisymmetric grid all cache-reusing jobs failed (cached 3D-expanded pressure field vs. 2D medium), leaving that run (`data/tussim_261001/itrusst_protocol1a_morris`, previous tissue ranges) with 47/150 results and no thermal elementary effects. The trajectory design (including every step's output_affix) is written to `data/tussim/itrusst_protocol1a_morris/morris_design.csv`; `PRESTUS_paramvar_10_morris_analysis_plot.m` joins it with the per-evaluation result CSVs and computes μ, μ\*, σ and n per property, tolerating a partially completed design.

---

## Outputs and Visualisation

Simulation results for each subject and configuration were written to CSV tables (`sub-XXX_phantom[_affix].csv`) located in per-setup output directories. These tables contain the full set of acoustic and thermal metrics enumerated above. For the pipeline comparison, bar plots of maxT and of Ipa_target alongside Isppa_brain were generated across the six configurations for each benchmark. For the sensitivity analysis, metrics were plotted as min–max ranges with individual data points overlaid, colour-coded by tissue type (brain: grey; skull: red). For the skull thickness analysis, Ipa_target and Isppa_brain were plotted as a continuous function of skull thickness in millimetres.
