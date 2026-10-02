# Methods

## Overview

Transcranial focused ultrasound (TUS) simulations were conducted using PRESTUS (Pipeline for Running E-field and TUS Simulations; v05), an open-source MATLAB toolbox that integrates k-Wave pseudospectral acoustic simulation with Pennes bioheat equation thermal modelling. All analyses were performed within the PRESTUS_itrusst analysis environment, which extends PRESTUS with a structured benchmarking framework to evaluate the sensitivity of simulation outputs to assumptions about tissue acoustic and thermal properties. Three analyses were performed: (1) a pipeline comparison across six tissue property configurations, (2) a systematic one-at-a-time sensitivity analysis of individual tissue parameters, and (3) an evaluation of the effect of skull thickness on transcranial transmission.

---

## Simulation Framework

### Phantom geometry

Two-dimensional axisymmetric phantoms were constructed at 0.5 mm isotropic resolution and stored as NIfTI volumes in SimNIBS m2m-folder format (`final_tissues.nii.gz`). The active simulation domain spanned 70 × 120 mm (140 × 240 voxels laterally and axially), with a 20 mm perfectly-matched layer (PML; 40 voxels at 0.5 mm resolution) appended externally by k-Wave to absorb outgoing waves.

Two standard benchmarks were used. **Benchmark 1** modelled a water–soft tissue interface, with soft tissue (label 1) beginning at 30 mm depth from the proximal domain edge. **Benchmark 2** modelled a transcranial geometry comprising a water column, a hemispherical skull annulus (inner radius 68.5 mm, outer radius 75 mm from the geometric centre), and brain tissue filling the interior. Both phantoms were treated as 2D axisymmetric (cylindrical symmetry) to reduce computation time. For a subset of analyses, fully three-dimensional simulations were run on a realistic head segmentation (sub-602) obtained with SimNIBS 4 CHARM segmentation from T1- and T2-weighted MRI, simulated at 0.75 mm resolution using C++ GPU code.

### Transducer

A single-element annular transducer was modelled (outer diameter 64 mm, radius of curvature 64 mm, inner diameter 0 mm) operating at a nominal focal depth of 64 mm. The transducer was positioned at the proximal edge of the active domain (1 mm offset; voxel coordinate [70, 2] at 0.5 mm resolution), with the geometric focus placed at voxel [70, 130] (64 mm axial depth).

### Transducer calibration

Prior to each phantom simulation, source amplitude was calibrated by running a free-water acoustic simulation using `itrusst_calibrate`. The peak spatial-peak pulse-average intensity (Isppa) along the axial profile was extracted from the simulated pressure field as:

> Isppa = p²_max / (2ρc) × 10⁻⁴  [W/cm²]

A linear correction factor was then applied to scale the source pressure amplitude (`elem_amp`) to achieve the desired target intensity (60 W/cm² or 30 W/cm² depending on protocol variant). The resulting calibrated amplitude and phase were used as transducer source parameters in the subsequent tissue simulation.

### Acoustic simulation

Acoustic pressure fields were computed using k-Wave (v1.4), a pseudospectral time-domain solver implementing the linearised pressure–velocity equations in heterogeneous media. Simulations used a Courant-Friedrichs-Lewy (CFL) stability number of 0.15 and were run in MATLAB CPU mode (2D phantoms) or C++ CPU/GPU mode (3D layered). The output of each acoustic simulation was the spatial maximum pressure field p_max(x) across all time steps, from which the following metrics were derived per tissue compartment:

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

For each configuration, acoustic and thermal simulations were run for Benchmark 1 (soft tissue only) and Benchmark 2 (skull + brain) at a target intensity of 60 W/cm². Simulation outputs (maxT, Ipa_target) were extracted from per-subject CSV tables and compared across configurations.

Attenuation is specified at 500 kHz and converted to PRESTUS' 1 MHz coefficient α₀ = α(500 kHz) / 0.5^y (dB/(MHz^y·cm); y = 1 for brain and skull). Both are listed below.

**Brain properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | Absorption fraction |
|---|---|---|---|---|---|---|---|---|
| setup 1 | 1520 | 1040 | 0.412 | 0.824 | 0.48 | 3525 | 196 | 1.0 |
| setup 2 | 1568 | 1050 | 0.206 | 0.412 | 0.55 | 3718 | 809 | 0.8 |
| reference | 1546 | 1045 | 0.24 | 0.48 | 0.51 | 3630 | 559 | 0.9 |
| PRESTUS | 1546 | 1046 | 0.26 | 0.52 | 0.51 | 3630 | 559 | 1.0 |
| BabelBrain | 1562 | 1041 | 0.24 | 0.48 | 0.51 | 3630 | 559 | 0.85 |
| kPlan | 1550 | 1045 | 0.26 | 0.52 | 0.51 | 3600 | 0 | 1.0 |

**Skull properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α @ 500 kHz (dB/cm) | α₀ (dB/(MHz·cm)) | k (W/m/K) | c_p (J/kg/K) | Perfusion (mL/min/kg) | Absorption fraction |
|---|---|---|---|---|---|---|---|---|
| setup 1 | 2871 | 1178 | 10.2 | 20.4 | 0.31 | 826 | 0 | 1.0 |
| setup 2 | 2190 | 2100 | 3.22 | 6.44 | 0.36 | 2524 | 30 | 0.16 |
| reference | 2800 | 1908 | 6.65 | 13.3 | 0.315 | 1793 | 20 | 0.58 |
| PRESTUS | 2800 | 1850 | 6.65 | 13.3 | 0.32 | 1793 | 20 | 0.28 |
| BabelBrain | 2800 | 1850 | 6.65 | 13.3 | 0.315 | 1793 | 20 | 0.163 |
| kPlan | 2800 | 1850 | 6.65 | 13.3 | 0.32 | 1300 | 0 | 1.0 |

---

## Analysis 2: Sensitivity Analysis — One-at-a-Time Parameter Perturbation

To isolate the contribution of individual tissue parameters to simulation outputs, a one-at-a-time (OAT) sensitivity analysis was performed. Each parameter was independently set to its literature-reported lower or upper bound while all other parameters were held at midpoint (baseline) values. Seven acoustic and thermal parameters were perturbed for each of two tissue types (brain, skull):

| Parameter | Symbol | Unit |
|---|---|---|
| Sound speed | c | m/s |
| Density | ρ | kg/m³ |
| Acoustic attenuation coefficient | α₀ | dB/(MHz^y·cm) |
| Thermal conductivity | k | W/m/K |
| Specific heat capacity | c_p | J/kg/K |
| Blood perfusion | ω | mL/min/kg |
| Absorption fraction | — | — |

Parameter bounds were drawn from a literature-informed tissue property table (`tissue_parameters.csv`). Brain parameters were perturbed for Benchmark 1 (sub-001); both brain and skull parameters were perturbed for Benchmark 2 (sub-002). For each perturbation, a separate simulation was run with a unique output filename affix encoding the tissue, parameter, and bound level (e.g., `_b_c_1` = brain sound speed lower bound; `_s_a0_2` = skull attenuation upper bound). Three protocol variants were evaluated in this sensitivity analysis:

| Variant | Frequency | Target intensity |
|---|---|---|
| itrusst_protocol1a | 500 kHz | 60 W/cm² |
| itrusst_protocol1b | 500 kHz | 30 W/cm² |
| itrusst_protocol1c | 250 kHz | 60 W/cm² |

Output metrics extracted for each perturbation included Ipa_target, maxT, Isppa_brain, riseT_brain, and real focal distance. All simulations used C++ CPU code in single-precision floating point.

---

## Analysis 3: Skull Thickness Sensitivity

To characterise the relationship between skull thickness and transcranial acoustic transmission, 23 skull phantoms were constructed with uniform cortical bone thickness from 1 to 12 mm in 0.5 mm steps (one voxel per step at 0.5 mm resolution; subjects sub-301 to sub-312 for 1–12 mm and sub-351 to sub-361 for 1.5–11.5 mm, see `code/functions/itrusst_thickness_subjects.m`). For each phantom, the spherical skull cap was placed at a fixed outer radius of curvature of 75 mm (outer apex at 30 mm depth), with the inner radius adjusted to achieve the target thickness. All tissue properties were set to midpoint values. Simulations used C++ CPU code at 0.5 mm resolution. The primary outcomes were Ipa_target (at the fixed target) and Isppa_brain (maximum in brain) as a function of skull thickness, quantifying acoustic insertion loss attributable to bone.

As a follow-up, skull density was varied across its bottom/mid/top levels on every thickness phantom (all other properties at midpoint, heating enabled, shortened 20 s thermal protocol; `PRESTUS_thickness_04_density_sim.m` / `_05_density_plot.m`), to test whether the effect of skull density depends on skull thickness.

---

## Analysis 4: Full-Factorial Constellation Verification via Acoustic/Thermal Decoupling

*Status: exploratory, not part of the main analysis (scripts in `code/explore/`). **The design described below (449 evaluations, acoustic/thermal caching, batched jobs) is not what the code implements.** `PRESTUS_paramvar_08_constellation_factorial_sim.m` currently runs a 52-evaluation 2×2 crossover of skull density with each of the other 13 properties (analysed by `PRESTUS_paramvar_13_density_crossover_plot.m`), without acoustic caching or batching. The Morris design (`PRESTUS_paramvar_09_morris_interactions_sim.m`) is analysed by `PRESTUS_paramvar_10_morris_analysis_plot.m`. This section needs rewriting before use.*

The pipeline-comparison configurations (Analysis 1) and the one-at-a-time sensitivity sweep (Analysis 2) each probe the tissue-property space in a limited way: the former evaluates only two literature-bound "constellations" (all seven properties per tissue simultaneously at their setup 1 or setup 2 bound), and the latter perturbs one property at a time from the midpoint, holding all others fixed. Neither directly establishes whether combining every property's individually assigned bound actually produces the expected combined (maximal/minimal) outcome, since one-at-a-time effects need not sum linearly when parameters are varied jointly.

A preliminary check compared, for each combined "constellation", the actual simulated outcome against the sum of the seven (Benchmark 1) or fourteen (Benchmark 2) individual one-at-a-time effects from Analysis 2 relative to the midpoint configuration. Every constellation moved each outcome (riseT_brain, riseT_skull, Isppa_brain, Ipa_target) in the direction predicted by the additive sum, but the magnitude deviated substantially: heating outcomes were super-additive (actual change 1.6–2.0× the additive prediction, i.e. combining bounds heats more than one-at-a-time effects would suggest), while intensity outcomes were close to additive with no skull present but sub-additive with skull present (actual change 0.65–0.79× predicted), consistent with attenuation acting multiplicatively (exponentially) rather than additively along the propagation path. This motivated evaluating the full combinatorial parameter space directly rather than relying on the additive approximation.

**Method.** A two-level (literature bottom/top bound) factorial design over all fourteen Benchmark 2 (brain+skull) tissue properties is run by exploiting a structural property of the PRESTUS pipeline: the acoustic (k-Wave) simulation depends only on sound speed, density and attenuation coefficient per tissue (confirmed both empirically, via zero measured effect of any other property on Isppa_brain/Ipa_target in Analysis 2, and in the toolbox source, where thermal-only properties are explicitly excluded from the medium struct passed to the acoustic solver), while the remaining four properties per tissue (thermal conductivity, specific heat capacity, perfusion, absorption fraction) affect only the subsequent bio-heat thermal simulation, which consumes the acoustic pressure field as a fixed input. The acoustic simulation is therefore run once per unique combination of the acoustic-relevant properties and cached; each combination of the remaining thermal-relevant properties is then evaluated by reloading the cached acoustic field and re-running only the thermal simulation with the new thermal property values, reproducing the exact joint outcome (not a statistical approximation) without repeating the acoustic solve. Benchmark 1 (brain only) is not part of this design - its smaller 2⁷ = 128-combination space was used earlier only to cross-check effect directions between Analyses 2 and 4.

A full 2¹⁴ = 16,384-combination design would require 16,448 simulations, far exceeding the study's job-count ceiling (≤1,000, ideally ≤128) even before scheduler queue limits are considered. This is not an exhaustive design: it instead exhaustively covers two *marginal* spaces at representative anchor points of the other dimension, in three complementary parts totalling 449 evaluations: (A) all 2⁶ = 64 acoustic combinations, each run with heating enabled at midpoint thermal properties - since these 64 acoustic simulations are required for caching regardless, running heating on them adds no extra evaluations and yields the exact 6-factor acoustic interaction structure at representative thermal settings; (B) one additional "mid-acoustic" anchor (all six acoustic properties at their literature midpoint) crossed with the full 2⁸ = 256-combination thermal factorial, giving the exact 8-factor thermal interaction structure at a representative acoustic point; and (C) the two already-established worst-case (setup 1) and best-case (setup 2) thermal bound patterns (Analysis 1; `data/tissue_setups.csv`) applied to each of the 64 part-(A) acoustic combinations, testing whether that worst/best-case assignment holds regardless of acoustic conditions. The full 64 × 256 acoustic × thermal joint grid is not covered - a gap the preliminary additive-vs-actual check above already bounds as a quantified, modest deviation rather than leaving unconstrained.

**Batching.** PRESTUS's own job dispatcher (`hpc_submit_job.m`) submits exactly one HPC job per parameter set, with no native array-job or batching support. To keep the job count within budget, the 449 evaluations are grouped into batches of `combos_per_job` (4 by default) sequential, in-process evaluations per submitted job, calling `prestus_pipeline()` directly rather than the per-call dispatcher. Because every thermal-only evaluation depends on its acoustic combination's cache already being written, batches run in two sequential phases: phase 1 (the 65 acoustic-defining/anchor evaluations of parts A and B, 17 batches at the default batch size) must complete before phase 2 (the 384 thermal-only evaluations of parts B and C, 96 batches) is submitted - 113 jobs total at the default batch size. `modules.run_nifti_creation` and `modules.generate_report` (both enabled by default) are disabled for every evaluation, since only the numeric per-combination CSV output table is used; `io.save_thermal_matrices` is likewise disabled throughout, and `io.save_acoustic_matrices` is enabled only for the 65 phase-1 evaluations whose cache phase 2 depends on.

To keep each evaluation minimal, the thermal simulation protocol is shortened from itrusst_protocol's default 80 s pulse-train duration + 80 s post-stimulation cool-off (4,080 thermal timesteps) to a 20 s pulse-train duration with no post-stimulation cool-off (1,000 timesteps, ≈24.5% of the original), since this design only requires riseT_brain/riseT_skull/Isppa_brain/Ipa_target - reached at or near the end of the stimulation period - not a full CEM43 dose calculation over an observed cool-down; this override is scoped to this analysis only and is not applied to the shared `config_itrusst_protocol.yaml` used by Analyses 1-3. Combined with the acoustic/thermal decoupling, per-evaluation cost is estimated (grounded in the one available fresh, current-grid timing reference: ≈12 s acoustic, ≈109 s thermal per combination at the original 80 s/80 s protocol, linearly rescaled by timestep count) at ≈51 s (acoustic-only) and ≈65 s (acoustic+heating or thermal-only), giving ≈8.1 h of total compute across all 449 evaluations. At the default batch size, each of the 113 submitted jobs runs ≈4 evaluations sequentially (≈4.3 min/job); if all batches within a phase run concurrently, wall-clock is on the order of ≈4-5 min per phase (≈10 min total for both phases in sequence) at full cluster parallelism, or correspondingly longer under limited concurrency (e.g. ≈30 min at 20-way parallelism). These are order-of-magnitude projections from a single timing reference and a linear timestep-scaling assumption, not a benchmarked measurement of this specific design. Before committing the full design, a small pilot subset (`pilot_mode = true`) must be run and cross-checked against an independently run, non-decoupled, non-cached, non-batched full simulation of the same property combinations (same shortened protocol, `io.acoustic_cache_affix` unset), to confirm the caching mechanism reproduces identical results.

**Complementary design: Morris elementary effects (interaction screening).** The design above does not cover the full acoustic × thermal joint grid; `PRESTUS_paramvar_09_morris_interactions_sim.m` targets that gap directly via a Morris (1991) elementary-effects screening design on the same 14-property, 2-level (bottom/top) space: r independent trajectories, each a "staircase" walk of 15 points (1 random baseline + all 14 properties flipped one at a time, in a random order grouped acoustic-factors-first to preserve the acoustic-caching benefit), from which each property's mean absolute elementary effect (μ\*, overall importance) and standard deviation of elementary effects (σ, a direct per-property measure of how much its effect depends on the other 13 properties' current values, i.e. interaction) are computed. At r=10 (default, 150 evaluations, ≈2.7 h total compute), cost scales linearly with r rather than exponentially with the number of properties, unlike a full factorial. Unlike the design above, no cross-evaluation acoustic-cache reuse is attempted beyond a trajectory's own 8 trailing thermal-only steps, since batches run in parallel with no ordering guarantee and a premature cache-reuse attempt fails silently (see the script's header for the specific mechanism) rather than erroring - every acoustic-cost evaluation therefore always requests a fresh solve, opportunistically benefiting from another evaluation's already-written cache when timing happens to favour it, but never depending on it. The trajectory design (including every step's output_affix) is written to `data/tussim/itrusst_protocol1a_morris/morris_design.csv`; the μ\*/σ analysis itself has not yet been run or scripted, since no simulation results exist yet.

---

## Outputs and Visualisation

Simulation results for each subject and configuration were written to CSV tables (`sub-XXX_phantom[_affix].csv`) located in per-setup output directories. These tables contain the full set of acoustic and thermal metrics enumerated above. For the pipeline comparison, bar plots of maxT and of Ipa_target alongside Isppa_brain were generated across the six configurations for each benchmark. For the sensitivity analysis, metrics were plotted as min–max ranges with individual data points overlaid, colour-coded by tissue type (brain: grey; skull: red). For the skull thickness analysis, Ipa_target and Isppa_brain were plotted as a continuous function of skull thickness in millimetres.
