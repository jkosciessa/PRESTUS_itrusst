% PRESTUS_ITRUSST_SETUPS_SKULL
%
% Parallel-coordinates plot of the tissue-property comparison table
% (Setup A / Setup B / PRESTUS / k-Plan* / BabelBrain*) across the seven
% reported properties (sound speed, density, attenuation, thermal
% conductivity, specific heat, perfusion, absorption fraction).
%
% Att. is the EFFECTIVE attenuation in dB/cm at the itrusst protocol1
% frequency (500 kHz), i.e. the raw config alpha_coeff (dB/(MHz cm), see
% data/configs/config_medium_*.yaml) scaled by freq_MHz^alpha_power = 0.5^1.
% This is why the Att. values here (12.2, 3.45, 6.65, ...) are exactly half
% of the configs' alpha_coeff (24.4, 6.9, 13.3, ...). PRESTUS_itrusst_setups_
% brain.m instead plots the raw, unscaled alpha_coeff (dB/(MHz cm)) for
% brain Att. - the two Att. columns are NOT in the same units.
%
% Each property gets its own vertical axis, scaled independently to its
% own min/max across setups (classic parallel-coordinates convention),
% so that a comparison of a fast-varying property (e.g. perfusion) is not
% swamped by a slow-varying one (e.g. sound speed) on a shared axis.
% Connecting lines trace each setup across all properties.
%
% Missing values (k-Plan* perfusion, reported as "/") are skipped: the
% line breaks across that property instead of being interpolated. This
% reflects k-Plan's own defaults, which do not document a perfusion value;
% config_medium_kPlan.yaml substitutes 0 so PRESTUS has a numeric value to
% simulate with, but that substitution is not shown here.
%
% NOTE: Setup A/B's density was corrected (here and in data/configs/
% config_medium_conservative.yaml / config_medium_liberal.yaml) after
% checking the OAT sensitivity sweep (data/tussim/itrusst_protocol1a/
% sub-002, via PRESTUS_paramvar_06_impact_direction_plot.m). Setup A had
% originally been assigned the top bound (1908 kg/m^3) and Setup B the
% bottom bound (1178 kg/m^3), following the literature range mechanically -
% but a LOWER skull density INCREASES riseT_brain by 25% (bottom vs. top
% bound), so 1178 is actually the heating-maximising ("conservative")
% value. Setup A now uses 1178 and Setup B 1908 (the reverse of the
% original assignment) - this is the largest of the two directional
% corrections found across brain and skull (the brain-side one, in
% PRESTUS_itrusst_setups_brain.m, is absorption fraction).
%
% Two figures are produced:
%   itrusst_table_parallel          - all axes oriented min (bottom) to max (top)
%   itrusst_table_parallel_flipped  - same data, but the Density, Thermal
%                                      conduct., Specific Heat and Perf. axes
%                                      are flipped (max at bottom, min at top)
%
% Output: figures/table_comparison/itrusst_table_parallel*.{png,eps}

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures', 'table_comparison');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.code    = fullfile(rootpath, 'code'); addpath(pn.code)

%% ── Table data ──────────────────────────────────────────────────────────────
% Columns: Sound Speed, Density, Att. (dB/cm, 500 kHz), Thermal conduct.,
%          Specific Heat, Perf., Abs. %

property_labels = {'Sound Speed', 'Density', 'Att. (dB/cm)', 'Thermal conduct.', ...
                    'Specific Heat', 'Perf.', 'Abs. %'};

setup_labels = {'Setup A', 'Setup B', 'PRESTUS', 'k-Plan*', 'BabelBrain*'};

data = [ ...
    2822, 1178, 12.2,  0.31,  1313, 10, 1.0;  ... % Setup A (density corrected, see NOTE above)
    2257, 1908, 3.45,  0.32,  2626, 30, 0.16; ... % Setup B (density corrected, see NOTE above)
    2800, 1850, 6.65,  0.32,  1793, 20, 0.28; ... % PRESTUS
    2800, 1850, 6.65,  0.32,  1300, NaN, 1.0; ... % k-Plan*  (perfusion not reported)
    2800, 1850, 6.65,  0.315, 1793, 20, 0.163 ... % BabelBrain*
];

setup_colors = { ...
    [0 0 0],           ... % Setup A     - black, distinguished from Setup B by marker
    [0 0 0],           ... % Setup B     - black, distinguished from Setup A by marker
    [0.00 0.45 0.74], ... % PRESTUS
    [0.47 0.67 0.19], ... % k-Plan*
    [0.49 0.18 0.56]  ... % BabelBrain*
};
%% ── Plot 1: all axes min-to-max ─────────────────────────────────────────────

make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
    [], pn, 'itrusst_table_parallel');

%% ── Plot 2: Density / Thermal conduct. / Specific Heat / Perf. flipped ─────

flip_cols = find(ismember(property_labels, {'Density', 'Thermal conduct.', 'Specific Heat', 'Perf.'}));

make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
    flip_cols, pn, 'itrusst_table_parallel_flipped');
