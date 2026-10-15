% PRESTUS_PARAMVAR_06_IMPACT_DIRECTION_PLOT
%
% For every one-at-a-time (OAT) tissue-property perturbation from Analysis 2
% (see README.md, Methods), determines whether moving a
% parameter from its literature "bottom" bound to its "top" bound (per
% data/tissue_parameters.csv) INCREASES or DECREASES these outcomes:
%   - Isppa_brain      : intensity (spatial-peak pulse-average), wherever it
%                        actually peaks within the brain mask
%   - Ipa_target       : target intensity, read at the FIXED nominal target
%                        voxel (kept for cross-reference; tracks Isppa_brain
%                        closely in this dataset)
%   - riseT_skull      : max. skull temperature rise (deg C) — independent
%                        whole-skull-mask maximum
%   - riseT_brain      : max. brain temperature rise (deg C) — an
%                        INDEPENDENT whole-brain-mask maximum. This is NOT
%                        pinned to the acoustic focus: the thermal hotspot
%                        can land on a different voxel than Isppa_brain's.
%   - riseT_focus      : temperature rise AT the empirical acoustic focus
%                        voxel (Ix_brain_vox/Iz_brain_vox — wherever Isppa
%                        actually peaked for that run), via
%                        load_focus_heating_metrics.m
%   - riseT_target     : temperature rise AT the FIXED nominal target voxel
%                        (focus_pos_vox_1/_3 — same location for every
%                        condition), via load_focus_heating_metrics.m
%   - real_focal_distance_mm : distance of the empirical acoustic focus from
%                        the transducer (apex-referenced)
%   - focal_vol_3dB_mm3      : ITRUSST −3 dB / half-intensity focal volume
%                        (i.e. the volumetric FWHM), over the whole brain mask
%
% riseT_brain / riseT_focus / riseT_target let us tell apart three distinct
% notions of "heating at the target": an independent whole-brain maximum,
% heating at wherever the beam actually focused, and heating at the fixed
% coordinate the treatment was planned for — these can diverge once tissue
% properties shift the beam.
%
% Uses the same per-subject OAT sweep CSVs (and, for riseT_focus/riseT_target,
% heatrise NIfTIs) as PRESTUS_paramvar_02_benchmark1_plot.m /
% PRESTUS_paramvar_02_benchmark2_plot.m
% (data/tussim/itrusst_protocol1a/sub-00{1,2}), loaded via
% load_paramvar_results.m / load_focus_heating_metrics.m:
%   sub-001 = Benchmark 1, brain only (no skull layer, so riseT_skull is NaN)
%   sub-002 = Benchmark 2, brain + skull
% Brain properties are therefore checked in BOTH benchmarks (useful for
% sanity-checking whether a direction is robust, or benchmark-dependent);
% skull properties, and riseT_skull, only exist for sub-002.
%
% "Direction" is defined as sign(top - bottom): positive means the outcome
% is higher at the top bound than at the bottom bound.
%
% Output:
%   figures/explore/impact_direction/impact_direction_Isppa_brain.{png,eps}
%   figures/explore/impact_direction/impact_direction_riseT_skull.{png,eps}
%   figures/explore/impact_direction/impact_direction_riseT_brain.{png,eps}
%   figures/explore/impact_direction/impact_direction_riseT_focus.{png,eps}
%   figures/explore/impact_direction/impact_direction_riseT_target.{png,eps}
%   figures/explore/impact_direction/impact_direction_real_focal_distance_mm.{png,eps}
%   figures/explore/impact_direction/impact_direction_focal_vol_3dB_mm3.{png,eps}
%   figures/explore/impact_direction/paramvar_impact_direction.csv  (full results table,
%                                                     all outcomes)

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
pn.figures = fullfile(rootpath, 'figures', 'explore', 'impact_direction');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.code    = fullfile(rootpath, 'code'); addpath(fullfile(pn.code, 'functions'))

%% ── OAT sweep definitions (matches PRESTUS_paramvar_01_benchmarks_sim.m) ────

property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
sweep_values      = 1:2; % 1 = bottom bound, 2 = top bound (tissue_parameters.csv)

%% ── Load OAT sweep results ───────────────────────────────────────────────

% sub-001: Benchmark 1, brain only
[res_b1, labels_b1] = load_paramvar_results(fullfile(pn.data, 'sub-001'), 1, ...
    {'brain'}, {'_b'}, property_labels, property_suffixes, sweep_values, property_symbols);
[focus_b1, ~] = load_focus_heating_metrics(fullfile(pn.data, 'sub-001'), 1, ...
    {'brain'}, {'_b'}, property_labels, property_suffixes, sweep_values, property_symbols);

% sub-002: Benchmark 2, brain + skull
[res_b2, labels_b2] = load_paramvar_results(fullfile(pn.data, 'sub-002'), 2, ...
    {'brain', 'skull'}, {'_b', '_s'}, property_labels, property_suffixes, sweep_values, property_symbols);
[focus_b2, ~] = load_focus_heating_metrics(fullfile(pn.data, 'sub-002'), 2, ...
    {'brain', 'skull'}, {'_b', '_s'}, property_labels, property_suffixes, sweep_values, property_symbols);

%% ── Build a combined results table (one row per condition) ─────────────────

conditions = [ ...
    build_rows(res_b1, focus_b1, labels_b1, {'brain'},         property_labels, 'sub-001 (brain only)'); ...
    build_rows(res_b2, focus_b2, labels_b2, {'brain', 'skull'}, property_labels, 'sub-002 (brain + skull)') ...
];

fprintf('\n%-24s %-10s | %9s %9s | %9s %9s %9s %9s | %9s %9s\n', ...
    'Parameter', 'Benchmark', 'dIsppa', 'dIpa', 'dRiseTsk', 'dRiseTbr', 'dRiseTfoc', 'dRiseTtgt', 'dDist', 'dFWHMvol');
for r = 1:numel(conditions)
    c = conditions(r);
    fprintf('%-24s %-10s | %+9.4f %+9.4f | %+9.4f %+9.4f %+9.4f %+9.4f | %+9.4f %+9.4f\n', ...
        c.label, regexp(c.benchmark, 'sub-\d+', 'match', 'once'), ...
        c.Isppa_delta, c.Ipa_delta, c.riseT_skull_delta, c.riseT_brain_delta, ...
        c.riseT_focus_delta, c.riseT_target_delta, c.focal_distance_delta, c.focal_vol_delta);
end

fprintf('\nSystematic pattern check (sign consistency across conditions where a property is tested more than once):\n');
for p = 1:numel(property_labels)
    mask = strcmp({conditions.property}, property_labels{p}) & strcmp({conditions.tissue}, 'brain');
    deltas = [conditions(mask).riseT_brain_delta];
    if numel(deltas) > 1 && all(~isnan(deltas))
        signs = sign(deltas);
        fprintf('  brain %-24s riseT_brain deltas: %s -> %s\n', property_labels{p}, ...
            sprintf('%+.4f ', deltas), ternary(all(signs == signs(1)), 'CONSISTENT', 'SIGN FLIP (check magnitude - may be noise)'));
    end
end

%% ── Write results table to CSV ───────────────────────────────────────────

T = struct2table(conditions);
writetable(T, fullfile(pn.figures, 'paramvar_impact_direction.csv'));

%% ── Tornado plots: one per outcome, sorted by |delta|, coloured by sign ────

plot_impact_tornado(conditions, 'Isppa_delta', 'Isppa\_brain \Delta (top - bottom)', ...
    pn, 'impact_direction_Isppa_brain');

plot_impact_tornado(conditions, 'riseT_skull_delta', 'riseT\_skull \Delta (top - bottom), ^{\circ}C', ...
    pn, 'impact_direction_riseT_skull');

plot_impact_tornado(conditions, 'riseT_brain_delta', 'riseT\_brain (whole-brain max) \Delta (top - bottom), ^{\circ}C', ...
    pn, 'impact_direction_riseT_brain');

plot_impact_tornado(conditions, 'riseT_focus_delta', 'riseT\_focus (at empirical acoustic focus) \Delta (top - bottom), ^{\circ}C', ...
    pn, 'impact_direction_riseT_focus');

plot_impact_tornado(conditions, 'riseT_target_delta', 'riseT\_target (at fixed nominal target) \Delta (top - bottom), ^{\circ}C', ...
    pn, 'impact_direction_riseT_target');

plot_impact_tornado(conditions, 'focal_distance_delta', 'real focal distance \Delta (top - bottom), mm', ...
    pn, 'impact_direction_real_focal_distance_mm');

plot_impact_tornado(conditions, 'focal_vol_delta', 'focal\_vol\_3dB (FWHM volume) \Delta (top - bottom), mm^3', ...
    pn, 'impact_direction_focal_vol_3dB_mm3');

%% ── Local functions ──────────────────────────────────────────────────────

function rows = build_rows(results, focus_results, xlabels, tissue_labels, property_labels, benchmark_name)
% Reshape a load_paramvar_results() struct + a load_focus_heating_metrics()
% struct (both 2 x n_cond arrays, row 1 = bottom bound, row 2 = top bound)
% into one struct-array row per condition. Iterates tissue (outer) x
% property (inner), matching the exact order load_paramvar_results.m used
% to build xlabels/the result columns.
    n_cond = numel(tissue_labels) * numel(property_labels);
    rows(n_cond, 1) = struct('benchmark', '', 'label', '', 'tissue', '', 'property', '', ...
        'Isppa_bottom', NaN, 'Isppa_top', NaN, 'Isppa_delta', NaN, ...
        'Ipa_bottom', NaN, 'Ipa_top', NaN, 'Ipa_delta', NaN, ...
        'riseT_skull_bottom', NaN, 'riseT_skull_top', NaN, 'riseT_skull_delta', NaN, ...
        'riseT_brain_bottom', NaN, 'riseT_brain_top', NaN, 'riseT_brain_delta', NaN, ...
        'riseT_focus_bottom', NaN, 'riseT_focus_top', NaN, 'riseT_focus_delta', NaN, ...
        'riseT_target_bottom', NaN, 'riseT_target_top', NaN, 'riseT_target_delta', NaN, ...
        'focal_distance_bottom', NaN, 'focal_distance_top', NaN, 'focal_distance_delta', NaN, ...
        'focal_vol_bottom', NaN, 'focal_vol_top', NaN, 'focal_vol_delta', NaN);

    i_cond = 0;
    for t = 1:numel(tissue_labels)
        for p = 1:numel(property_labels)
            i_cond = i_cond + 1;
            rows(i_cond).benchmark = benchmark_name;
            rows(i_cond).label     = xlabels{i_cond};
            rows(i_cond).tissue    = tissue_labels{t};
            rows(i_cond).property  = property_labels{p};

            rows(i_cond).Isppa_bottom = results.Isppa_brain(1, i_cond);
            rows(i_cond).Isppa_top    = results.Isppa_brain(2, i_cond);
            rows(i_cond).Isppa_delta  = results.Isppa_brain(2, i_cond) - results.Isppa_brain(1, i_cond);

            rows(i_cond).Ipa_bottom = results.Ipa_target(1, i_cond);
            rows(i_cond).Ipa_top    = results.Ipa_target(2, i_cond);
            rows(i_cond).Ipa_delta  = results.Ipa_target(2, i_cond) - results.Ipa_target(1, i_cond);

            rows(i_cond).riseT_skull_bottom = results.riseT_skull(1, i_cond);
            rows(i_cond).riseT_skull_top    = results.riseT_skull(2, i_cond);
            rows(i_cond).riseT_skull_delta  = results.riseT_skull(2, i_cond) - results.riseT_skull(1, i_cond);

            rows(i_cond).riseT_brain_bottom = results.riseT_brain(1, i_cond);
            rows(i_cond).riseT_brain_top    = results.riseT_brain(2, i_cond);
            rows(i_cond).riseT_brain_delta  = results.riseT_brain(2, i_cond) - results.riseT_brain(1, i_cond);

            rows(i_cond).riseT_focus_bottom = focus_results.riseT_focus(1, i_cond);
            rows(i_cond).riseT_focus_top    = focus_results.riseT_focus(2, i_cond);
            rows(i_cond).riseT_focus_delta  = focus_results.riseT_focus(2, i_cond) - focus_results.riseT_focus(1, i_cond);

            rows(i_cond).riseT_target_bottom = focus_results.riseT_target(1, i_cond);
            rows(i_cond).riseT_target_top    = focus_results.riseT_target(2, i_cond);
            rows(i_cond).riseT_target_delta  = focus_results.riseT_target(2, i_cond) - focus_results.riseT_target(1, i_cond);

            rows(i_cond).focal_distance_bottom = results.real_focal_distance_mm(1, i_cond);
            rows(i_cond).focal_distance_top    = results.real_focal_distance_mm(2, i_cond);
            rows(i_cond).focal_distance_delta  = results.real_focal_distance_mm(2, i_cond) - results.real_focal_distance_mm(1, i_cond);

            rows(i_cond).focal_vol_bottom = results.focal_vol_3dB_mm3(1, i_cond);
            rows(i_cond).focal_vol_top    = results.focal_vol_3dB_mm3(2, i_cond);
            rows(i_cond).focal_vol_delta  = results.focal_vol_3dB_mm3(2, i_cond) - results.focal_vol_3dB_mm3(1, i_cond);
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function plot_impact_tornado(conditions, field, x_label, pn, fig_name)
% Horizontal "tornado" bar chart of conditions(:).(field), sorted by
% |value| (largest at top), coloured red (increases outcome) / blue
% (decreases outcome). Rows with NaN (e.g. riseT_skull for sub-001, which
% has no skull layer) are dropped.
    values = [conditions.(field)]';
    keep = ~isnan(values);
    conditions = conditions(keep);
    values = values(keep);

    labels = arrayfun(@(c) sprintf('%s (%s)', c.label, regexp(c.benchmark, 'sub-\d+', 'match', 'once')), ...
        conditions, 'UniformOutput', false);

    [~, order] = sort(abs(values), 'ascend'); % ascend so barh plots largest at top
    values = values(order);
    labels = labels(order);

    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.2 0.1 0.35 0.8]);
    hold on;

    n = numel(values);
    colors = repmat([0.00 0.45 0.74], n, 1); % blue = decreases outcome
    colors(values > 0, :) = repmat([0.85 0.10 0.10], sum(values > 0), 1); % red = increases outcome

    for i = 1:n
        barh(i, values(i), 'FaceColor', colors(i, :), 'EdgeColor', 'none', 'BarWidth', 0.7);
    end

    plot([0 0], [0.3 n + 0.7], '-k', 'LineWidth', 1);
    set(gca, 'YTick', 1:n, 'YTickLabel', labels, 'TickLabelInterpreter', 'tex', 'FontSize', 10);
    ylim([0.3 n + 0.7]);
    xlabel(x_label, 'Interpreter', 'tex');
    set(gca, 'Layer', 'top', 'Box', 'off');
    title('bottom \rightarrow top bound: red = outcome increases, blue = outcome decreases', 'FontSize', 10);

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
