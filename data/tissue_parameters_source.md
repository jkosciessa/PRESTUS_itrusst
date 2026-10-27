# Tissue parameter source table

Source of `tissue_parameters.csv` (updated table, 2026-10-27; previous revisions 2026-10-16, 2026-10-08, 2026-10-01). Attenuation and absorption are given at 500 kHz in dB/cm;
`tissue_parameters.csv` stores PRESTUS' 1 MHz coefficients `alpha_coeff = alpha_500kHz / 0.5^alpha_power` and `absorption_coeff = alpha_abs_500kHz / 0.5^alpha_power`
(dB/(MHz^y cm); y = 2 for water, y = 1 for brain and skull). Levels map LOW/reference/HIGH -> bottom/mid/top.

| Tissue | Setting | c [m/s] | ρ [kg/m³] | α @ 500 kHz [dB/cm] | κ [W/m/K] | c_p [J/(kg K)] | ω [mL/min/kg] | α_abs @ 500 kHz [dB/cm] | Absorption fraction (deprecated) |
|---|---|---|---|---|---|---|---|---|---|
| Water | uniform | 1482 | 994 | 0.00055 | 0.6 | 4178 | 0 | 0.00055 | 1 |
| Skull | LOW | 2190 | 1080 † | 4.02 | 0.30 † | 826 † | 0 | 1.05 | 0.16 |
| Skull | reference | 2550 | 1775 | 6.0 | 0.315 | 1793 | 20 | 1.40 † | 0.58 |
| Skull | HIGH | 2871 | 1991 † | 11.7 | 0.36 | 2524 † | 30 | 4.02 | 1.0 |
| Brain | LOW | 1520 † | 1040 | 0.236 | 0.48 | 3525 | 196 | 0.045 † | 0.73 † |
| Brain | reference | 1560 | 1040 | 0.30 | 0.51 | 3630 | 559 | 0.093 † | 0.865 |
| Brain | HIGH | 1573 | 1050 | 0.493 | 0.55 | 3718 | 809 | 0.231 † | 1.0 |

† as marked in the source table (2026-10-16).

Heat deposition uses the absorption coefficient α_abs (PRESTUS `thermal.absorption_model: coefficient`, `absorption_coeff` in dB/(MHz^y cm), y as for attenuation). The absorption fraction column is deprecated: it is listed for comparison with earlier revisions and with tools that specify absorption as a fraction, and is not used in the simulations. It does not equal α_abs / α; the fractions implied by the coefficients are 0.26 / 0.23 / 0.34 (skull LOW / reference / HIGH) and 0.19 / 0.31 / 0.47 (brain).

Each HIGH absorption coefficient is at most the LOW attenuation of the same tissue (skull 4.02 = 4.02; brain 0.231 < 0.236), so absorption never exceeds attenuation in any one-at-a-time or setup combination. Skull LOW attenuation (4.02) is the low end of whole-skull attenuation in White et al. (2006), which also serves as a proxy for the high end of skull absorption (HIGH, 4.02). White et al. report 53 Np/m (4.60 dB/cm) at 0.548 MHz and a least-squares power-law exponent of 1.47 across 0.272 / 0.548 / 0.840 MHz; scaling the 0.548 MHz value with b = 1.47 gives 4.02 dB/cm at 500 kHz. Simulations run at 500 kHz with y = 1 for skull, so the 500 kHz value is kept and the 1 MHz coefficient is rescaled with y = 1 (`alpha_coeff` = `absorption_coeff` = 4.02 / 0.5 = 8.04 dB/(MHz cm)) rather than with y = 1.47. An earlier version of the 2026-10-16 table gave 3.22 dB/cm, which does not match this conversion, and the 2026-10-16 table gave 4.03 dB/cm (b = 1.43 for attenuation, 1.46 for absorption). Skull LOW (1.05) and reference (1.40) are the bottom and top of White et al. (1978; 2.1–2.8 dB/cm at 1 MHz). Skull LOW was 1.35 (longitudinal-wave absorption of Pinton et al., 2012; 2.7 dB/cm at 1 MHz) until 2026-10-19. Brain reference (0.093) matches Goss et al. (1979; via Duck 1990). Source of brain LOW (0.045) and HIGH (0.231): to be documented.

Skull reference sound speed, density and attenuation (2550 m/s, 1775 kg/m³, 6.0 dB/cm) are a project decision: the mean of the cortical (2800 m/s, 1850 kg/m³, 4 dB/cm) and trabecular (2300 m/s, 1700 kg/m³, 8 dB/cm) bone values of Aubry et al. (2022). The 2026-10-08 table used the cortical values alone; the other skull reference values follow the table.
