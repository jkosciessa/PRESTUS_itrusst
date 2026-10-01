% PRESTUS_ITRUSST_SETUPS_BRAIN
%
% Parallel-coordinates plot of the soft-tissue (brain) property comparison
% table (Setup A / Setup B / PRESTUS / BabelBrain / kPlan) across the
% seven reported properties (sound speed, density, attenuation, thermal
% conductivity, specific heat, perfusion, absorption fraction). Data taken
% from data/configs/config_medium_*.yaml (the ground truth for what was
% actually modelled), cross-checked against the "Brain properties by
% configuration" table in doc/PRESTUS_itrusst_methods.md (Analysis 1:
% Pipeline Benchmark). Setup A = conservative (maximises predicted heating)
% and Setup B = liberal (minimises predicted heating) - mirroring the
% Setup A / Setup B naming used for the skull properties in
% PRESTUS_itrusst_setups_skull.m. The "mid" configuration is omitted here
% (it is the midpoint of the literature range, not a distinct tool/
% benchmark setup).
%
% NOTE: two properties were corrected here (and in data/configs/
% config_medium_conservative.yaml / config_medium_liberal.yaml /
% doc/PRESTUS_itrusst_methods.md) after checking the OAT sensitivity sweep
% (data/tussim/itrusst_protocol1a/sub-00{1,2}, via PRESTUS_paramvar_06_
% impact_direction_plot.m), because the literature bottom/top bound had
% been assigned to conservative/liberal WITHOUT checking whether that
% direction actually maximises/minimises heating for that property:
%   - Density: Setup A uses the bottom bound (1041 kg/m^3), not the
%     1046 kg/m^3 originally used for "conservative". Raising density from
%     1041 to 1046 DECREASES riseT_brain (-0.0023 degC in Benchmark 1,
%     brain only; -0.0015 degC in Benchmark 2, brain+skull), so 1041 is the
%     heating-maximising ("conservative") value - a small effect (<1%).
%   - Absorption fraction: Setup A now uses the TOP bound (1.0), Setup B
%     the bottom bound (0.8) - the reverse of the originally assigned
%     0.8/1.0. Raising absorption fraction from 0.8 to 1.0 INCREASES
%     riseT_brain by +25% (Benchmark 1) / +0.9% (Benchmark 2), so 1.0 is
%     the heating-maximising value, not 0.8.
% Density's and absorption fraction's effect on Ipa_target were not used to
% judge direction: density's is negligible and flips sign between
% benchmarks (-0.008 vs +0.015); Ipa_target is driven by acoustic
% properties (sound speed, attenuation), not absorption fraction, which
% only affects the heat source term.
%
% Att. is the raw config alpha_coeff (dB/(MHz cm)), unscaled by frequency.
% This is NOT the same unit as PRESTUS_itrusst_setups_skull.m's Att.
% column, which instead plots the effective attenuation in dB/cm at the
% itrusst protocol1 frequency (500 kHz) - i.e. alpha_coeff x freq_MHz.
%
% Each property gets its own vertical axis, scaled independently to its
% own min/max across setups (classic parallel-coordinates convention),
% so that a comparison of a fast-varying property (e.g. perfusion) is not
% swamped by a slow-varying one (e.g. sound speed) on a shared axis.
% Connecting lines trace each setup across all properties.
%
% Missing values (k-Plan perfusion, reported as "/") are skipped: the line
% breaks across that property instead of being interpolated. This mirrors
% the skull table's treatment of k-Plan* perfusion: k-Plan's own defaults
% do not document a perfusion value; config_medium_kPlan.yaml substitutes 0
% (brain and skull alike) so PRESTUS has a numeric value to simulate with,
% but that substitution is not shown here. Specific heat (1300 for k-Plan
% skull) is the analogous case in the skull table where a tool default
% falls just outside the Setup A/B literature range.
%
% Two figures are produced:
%   itrusst_table_parallel_brain          - all axes oriented min (bottom) to max (top)
%   itrusst_table_parallel_brain_flipped  - same data, but the Sound Speed,
%                                            Density, Thermal conduct.,
%                                            Specific Heat and Perf. axes are
%                                            flipped (max at bottom, min at
%                                            top), so that Setup A and Setup
%                                            B - the two extremes of the
%                                            literature range - trace out a
%                                            clean envelope (A on top, B on
%                                            bottom) without crossing over
%                                            each other
%
% Output: figures/table_comparison/itrusst_table_parallel_brain*.{png,eps}

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures', 'table_comparison');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.code    = fullfile(rootpath, 'code'); addpath(fullfile(pn.code, 'functions'))

%% ── Table data ──────────────────────────────────────────────────────────────
% Columns: Sound Speed, Density, Att. (dB/MHz/cm), Thermal conduct.,
%          Specific Heat, Perf., Abs. %

property_labels = {'Sound Speed', 'Density', 'Att. (dB/MHz/cm)', 'Thermal conduct.', ...
                    'Specific Heat', 'Perf.', 'Abs. %'};

setup_labels = {'Setup A', 'Setup B', 'PRESTUS', 'BabelBrain', 'k-Plan'};

data = [ ...
    1520, 1041, 0.94, 0.48, 3583, 196, 1.00; ... % Setup A (conservative; density & absorption fraction corrected, see NOTE above)
    1568, 1046, 0.36, 0.55, 3696, 809, 0.80; ... % Setup B (liberal; density & absorption fraction corrected, see NOTE above)
    1546, 1046, 0.52, 0.51, 3630, 559, 1.00; ... % PRESTUS
    1562, 1041, 0.48, 0.51, 3630, 559, 0.85; ... % BabelBrain
    1550, 1045, 0.52, 0.51, 3600, NaN, 1.00  ... % k-Plan  (perfusion not reported)
];

setup_colors = { ...
    [0.70 0.00 0.00], ... % Setup A (conservative)
    [0.00 0.60 0.50], ... % Setup B (liberal)
    [0.00 0.45 0.74], ... % PRESTUS
    [0.49 0.18 0.56], ... % BabelBrain
    [0.47 0.67 0.19]  ... % k-Plan
};

%% ── Plot 1: all axes min-to-max ─────────────────────────────────────────────

make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
    [], pn, 'itrusst_table_parallel_brain');

%% ── Plot 2: Setup A / Setup B envelope, no crossover ────────────────────────
% Setup A is the heating-maximising bound on every property (top for Att.
% and Abs. %, bottom for everything else), so Sound Speed, Density, Thermal
% conduct., Specific Heat and Perf. are flipped to put Setup A consistently
% on top and Setup B consistently on the bottom across all seven properties
% (Att. and Abs. % are already unflipped-top for Setup A, so are left as-is).

flip_cols = find(ismember(property_labels, ...
    {'Sound Speed', 'Density', 'Thermal conduct.', 'Specific Heat', 'Perf.'}));

make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
    flip_cols, pn, 'itrusst_table_parallel_brain_flipped');
