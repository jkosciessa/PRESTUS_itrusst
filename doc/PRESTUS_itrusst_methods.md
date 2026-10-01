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

To quantify how assumptions about tissue acoustic and thermal properties propagate into simulation outputs, six property configurations were evaluated representing the defaults used by five TUS simulation tools/frameworks as well as a conservative (maximum attenuation) parameter set:

| Configuration | Rationale |
|---|---|
| **conservative** | High brain attenuation (α = 0.94 dB/MHz/cm), low perfusion (196 s⁻¹); maximises predicted heating |
| **liberal** | Low brain attenuation (α = 0.36 dB/MHz/cm), high perfusion (809 s⁻¹); minimises predicted heating |
| **mid** | Midpoint values across the literature range (α = 0.48, perfusion = 559 s⁻¹) |
| **PRESTUS** | PRESTUS toolbox defaults (α = 0.52, perfusion = 559 s⁻¹, absorption fraction = 1.0) |
| **BabelBrain** | BabelBrain defaults (sound speed = 1562 m/s, absorption fraction = 0.85) |
| **kPlan** | kPlan defaults (perfusion = 0, absorption fraction = 1.0) |

For each configuration, acoustic and thermal simulations were run for Benchmark 1 (soft tissue only; sub-001) and Benchmark 2 (skull + brain; sub-002) at a target intensity of 60 W/cm². Simulation outputs (maxT, Ipa_target) were extracted from per-subject CSV tables and compared across configurations.

The full tissue property values used in each configuration are summarised below.

**Brain properties by configuration:**

| Setup | Sound speed (m/s) | Density (kg/m³) | α (dB/MHz/cm) | k (W/m/K) | c (J/kg/K) | Perfusion (s⁻¹) | Absorption fraction |
|---|---|---|---|---|---|---|---|
| conservative | 1520 | 1046 | 0.94 | 0.48 | 3583 | 196 | 1.0 |
| liberal | 1568 | 1046 | 0.36 | 0.55 | 3696 | 809 | 0.8 |
| mid | 1546 | 1045 | 0.48 | 0.51 | 3630 | 559 | 0.9 |
| PRESTUS | 1546 | 1046 | 0.52 | 0.51 | 3630 | 559 | 1.0 |
| BabelBrain | 1562 | 1041 | 0.48 | 0.51 | 3630 | 559 | 0.85 |
| kPlan | 1550 | 1045 | 0.52 | 0.51 | 3600 | 0 | 1.0 |

Brain absorption fraction for conservative/liberal was corrected from the originally assigned literature bottom/top bound (0.8/1.0) to 1.0/0.8: the one-at-a-time sensitivity sweep (Analysis 2) shows a *higher* absorption fraction increases predicted brain heating (riseT_brain +25% for the top vs. bottom bound, Benchmark 1), so the bottom bound was the wrong choice for "conservative" (maximises heating). The corresponding skull property table (used in PRESTUS_itrusst_setups_skull.m) has an analogous correction for skull density, which was similarly assigned the wrong bound (a *lower* skull density increases heating by 25%, but conservative had been assigned the top bound).

---

## Analysis 2: Sensitivity Analysis — One-at-a-Time Parameter Perturbation

To isolate the contribution of individual tissue parameters to simulation outputs, a one-at-a-time (OAT) sensitivity analysis was performed. Each parameter was independently set to its literature-reported lower or upper bound while all other parameters were held at midpoint (baseline) values. Seven acoustic and thermal parameters were perturbed for each of two tissue types (brain, skull):

| Parameter | Symbol | Unit |
|---|---|---|
| Sound speed | c | m/s |
| Density | ρ | kg/m³ |
| Acoustic attenuation coefficient | α | dB/MHz/cm |
| Thermal conductivity | k | W/m/K |
| Specific heat capacity | c_p | J/kg/K |
| Blood perfusion | ω | s⁻¹ |
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

To characterise the relationship between skull thickness and transcranial acoustic transmission, twelve cylindrical skull phantoms were constructed with uniform cortical bone thickness ranging from 1 to 12 mm (subjects sub-301 to sub-312). For each phantom, the skull annulus was placed at a fixed outer radius of 75 mm from the geometric centre, with the inner radius adjusted to achieve the target thickness. All tissue properties were set to midpoint values. Simulations used C++ CPU code at 0.5 mm resolution. The primary outcome was Ipa_target at the brain focus as a function of skull thickness, quantifying acoustic insertion loss attributable to bone.

---

## Outputs and Visualisation

Simulation results for each subject and configuration were written to CSV tables (`sub-XXX_phantom_output_table[_affix].csv`) located in per-setup output directories. These tables contain the full set of acoustic and thermal metrics enumerated above. For the pipeline comparison, bar plots of maxT and Ipa_target were generated across the six configurations for each benchmark. For the sensitivity analysis, metrics were plotted as min–max ranges with individual data points overlaid, colour-coded by tissue type (brain: grey; skull: red). For the skull thickness analysis, Ipa_target was plotted as a continuous function of skull thickness in millimetres.
