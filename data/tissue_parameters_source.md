# Tissue parameter source table

Source of `tissue_parameters.csv` (updated table, 2026-10-08; previous revision 2026-10-01). Attenuation is given at 500 kHz in dB/cm;
`tissue_parameters.csv` stores PRESTUS' 1 MHz coefficient `alpha_coeff = alpha_500kHz / 0.5^alpha_power`
(dB/(MHz^y cm); y = 2 for water, y = 1 for brain and skull). Levels map LOW/reference/HIGH -> bottom/mid/top.

| Tissue | Setting | c [m/s] | ρ [kg/m³] | α @ 500 kHz [dB/cm] | κ [W/m/K] | c_p [J/(kg K)] | ω [mL/min/kg] | Absorption fraction |
|---|---|---|---|---|---|---|---|---|
| Water | uniform | 1482 | 994 | 0.00055 | 0.6 | 4178 | 0 | 1 |
| Skull | LOW | 2190 | 1080 † | 3.22 | 0.30 † | 826 † | 0 | 0.16 |
| Skull | reference | 2800 | 1850 | 4.0 | 0.315 | 1793 | 20 | 0.58 |
| Skull | HIGH | 2871 | 1991 † | 11.7 | 0.36 | 2524 † | 30 | 1.0 |
| Brain | LOW | 1520 † | 1040 | 0.236 | 0.48 | 3525 | 196 | 0.73 † |
| Brain | reference | 1560 | 1040 | 0.30 | 0.51 | 3630 | 559 | 0.865 |
| Brain | HIGH | 1568 | 1050 | 0.493 | 0.55 | 3718 | 809 | 1.0 |

† as marked in the source table (2026-10-08).
