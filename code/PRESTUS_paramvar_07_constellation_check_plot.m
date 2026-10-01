% PRESTUS_PARAMVAR_07_CONSTELLATION_CHECK_PLOT
%
% Checks whether combining every tissue property at its individually
% assigned "conservative"/"liberal" bound (data/configs/config_medium_
% conservative.yaml / config_medium_liberal.yaml) actually produces the
% expected combined outcome, or whether interactions between parameters
% make the true combined simulation deviate from what a naive sum of the
% one-at-a-time (OAT) effects would predict.
%
% Method: the OAT sweep (Analysis 2; data/tussim/itrusst_protocol1a/
% sub-00{1,2}) gives the effect of moving ONE property from mid to its
% bottom or top bound, with every other property held at mid. Summing these
% single-property effects - using whichever bound "conservative"/"liberal"
% assigns to each property - gives an ADDITIVE PREDICTION for the fully
% combined config. This is compared against the ACTUAL combined simulation
% (Analysis 1 pipeline comparison; data/tussim/medium_conservative,
% medium_liberal, medium_mid).
%
%   actual/predicted ~ 1  : additive - OAT sensitivity fully explains the
%                           combined outcome
%   actual/predicted > 1  : synergistic - combining the assigned bounds is
%                           MORE extreme than OAT would predict
%   actual/predicted < 1  : sub-additive - combining bounds is LESS extreme
%                           than OAT would predict (e.g. expected for
%                           attenuation-like effects, which are inherently
%                           multiplicative/exponential, not additive)
%   sign(actual) ~= sign(predicted) : the constellation does not even move
%                           the outcome in the intended direction - would
%                           mean Setup A/B ("conservative"/"liberal") do
%                           NOT reliably act as a worst/best case
%
% Bound assignment (bottom/top) per property is read directly from
% data/configs/config_medium_conservative.yaml / config_medium_liberal.yaml
% and matched against data/tissue_parameters.csv, so this check stays
% correct as those configs are corrected/updated - it does not hardcode
% which bound each property currently uses.
%
% NOTE: this script always reads the bound assignment from the CURRENT
% config_medium_conservative.yaml / config_medium_liberal.yaml, but the
% "actual" combined outcome comes from whatever simulation is currently on
% disk under data/tussim/medium_conservative / medium_liberal. As of the
% skull density / brain absorption fraction correction (see
% PRESTUS_itrusst_setups_brain.m / _skull.m), those two are OUT OF SYNC:
% the configs already reflect the corrected bounds, but the simulation
% outputs on disk still reflect the OLD (pre-correction) bounds for those
% two properties - so running this script now gives a mismatched,
% not-directly-interpretable comparison for skull density / brain
% absorption fraction specifically (the other 12 property/tissue
% combinations are unaffected and still valid). Re-run PRESTUS_pipeline_01
% _benchmark1_sim.m for 'conservative' and 'liberal' first so config and
% simulation agree again, then re-run this script.
%
% Output:
%   figures/paramvar/constellation_check.{png,eps}
%   figures/paramvar/paramvar_constellation_check.csv

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a');
pn.tussim  = fullfile(rootpath, 'data', 'tussim');
pn.configs = fullfile(rootpath, 'data', 'configs');
pn.figures = fullfile(rootpath, 'figures', 'paramvar');
if ~isfolder(pn.figures); mkdir(pn.figures); end
pn.code    = fullfile(rootpath, 'code'); addpath(pn.code)

%% ── OAT sweep definitions (matches PRESTUS_paramvar_01_benchmarks_sim.m) ────

property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
sweep_values      = 1:2; % 1 = bottom bound, 2 = top bound (tissue_parameters.csv)
metrics           = {'riseT_brain', 'riseT_skull', 'Isppa_brain', 'Ipa_target'};

params_table = readtable(fullfile(rootpath, 'data', 'tissue_parameters.csv'));

%% ── Load OAT sweep results (single-property effects, others held at mid) ───

[res_b1, ~] = load_paramvar_results(fullfile(pn.data, 'sub-001'), 1, ...
    {'brain'}, {'_b'}, property_labels, property_suffixes, sweep_values, property_symbols);
[res_b2, ~] = load_paramvar_results(fullfile(pn.data, 'sub-002'), 2, ...
    {'brain'}, {'_b'}, property_labels, property_suffixes, sweep_values, property_symbols);
[res_s2, ~] = load_paramvar_results(fullfile(pn.data, 'sub-002'), 2, ...
    {'skull'}, {'_s'}, property_labels, property_suffixes, sweep_values, property_symbols);

%% ── Read conservative/liberal bound assignment straight from the configs ───

setup_names = {'conservative', 'liberal'};
bound = struct();
for si = 1:numel(setup_names)
    yaml_path = fullfile(pn.configs, ['config_medium_' setup_names{si} '.yaml']);
    for p = 1:numel(property_labels)
        bound.(setup_names{si}).brain.(property_labels{p}) = ...
            config_bound_level(yaml_path, 'brain', property_labels{p}, params_table);
        bound.(setup_names{si}).skull.(property_labels{p}) = ...
            config_bound_level(yaml_path, 'skull', property_labels{p}, params_table);
    end
end

%% ── Benchmarks: OAT subject <-> combined-simulation subject ────────────────

benchmarks(1) = struct('name', 'Benchmark 1 (brain only)',  'res_brain', res_b1, 'res_skull', [],    'sub_combined', 101, 'has_skull', false);
benchmarks(2) = struct('name', 'Benchmark 2 (brain+skull)', 'res_brain', res_b2, 'res_skull', res_s2, 'sub_combined', 102, 'has_skull', true);

%% ── Compute predicted (additive) vs. actual combined delta ─────────────────

results = struct('benchmark', {}, 'setup', {}, 'metric', {}, ...
    'mid', {}, 'actual', {}, 'actual_delta', {}, 'predicted_delta', {}, 'ratio', {});

for bi = 1:numel(benchmarks)
    bm = benchmarks(bi);
    mid_vals = load_combined_outcome(pn.tussim, 'medium_mid', bm.sub_combined, metrics);

    for si = 1:numel(setup_names)
        setup_name  = setup_names{si};
        actual_vals = load_combined_outcome(pn.tussim, ['medium_' setup_name], bm.sub_combined, metrics);

        for mi = 1:numel(metrics)
            metric = metrics{mi};
            if strcmp(metric, 'riseT_skull') && ~bm.has_skull
                continue; % no skull layer in Benchmark 1
            end

            mid_val      = mid_vals.(metric);
            actual_delta = actual_vals.(metric) - mid_val;

            pred_delta = 0;
            for p = 1:numel(property_labels)
                lvl = bound.(setup_name).brain.(property_labels{p});
                if ~isnan(lvl)
                    pred_delta = pred_delta + (bm.res_brain.(metric)(lvl, p) - mid_val);
                end
            end
            if bm.has_skull
                for p = 1:numel(property_labels)
                    lvl = bound.(setup_name).skull.(property_labels{p});
                    if ~isnan(lvl)
                        pred_delta = pred_delta + (bm.res_skull.(metric)(lvl, p) - mid_val);
                    end
                end
            end

            row = struct('benchmark', bm.name, 'setup', setup_name, 'metric', metric, ...
                'mid', mid_val, 'actual', actual_vals.(metric), 'actual_delta', actual_delta, ...
                'predicted_delta', pred_delta, 'ratio', actual_delta / pred_delta);
            results(end+1) = row; %#ok<AGROW>
        end
    end
end

%% ── Print + save results ────────────────────────────────────────────────

fprintf('\n%-28s %-14s %-14s %10s %10s %10s\n', ...
    'Benchmark', 'Setup', 'Metric', 'actual_d', 'predict_d', 'ratio');
for r = 1:numel(results)
    c = results(r);
    fprintf('%-28s %-14s %-14s %+10.4f %+10.4f %10.2f\n', ...
        c.benchmark, c.setup, c.metric, c.actual_delta, c.predicted_delta, c.ratio);
end

sign_mismatch = arrayfun(@(c) sign(c.actual_delta) ~= sign(c.predicted_delta), results);
if any(sign_mismatch)
    fprintf('\nWARNING: %d condition(s) show a SIGN MISMATCH (constellation does not move the\n', sum(sign_mismatch));
    fprintf('outcome in the OAT-predicted direction) - see rows above with opposite-sign actual/predicted.\n');
else
    fprintf('\nNo sign mismatches: every constellation moves each outcome in the OAT-predicted\n');
    fprintf('direction. Magnitudes still deviate from additive where ratio != 1 (see above).\n');
end

T = struct2table(results);
writetable(T, fullfile(pn.figures, 'paramvar_constellation_check.csv'));

%% ── Scatter plot: predicted (additive) vs. actual combined delta ───────────

plot_constellation_scatter(results, pn, 'constellation_check');

%% ── Local functions ──────────────────────────────────────────────────────

function level = config_bound_level(yaml_path, tissue, property_name, params_table)
% Reads <property_name> for <tissue> from a config_medium_*.yaml file and
% returns 1 if it matches the literature "bottom" bound (tissue_parameters
% .csv), 2 if it matches "top", or NaN (with a warning) if it matches
% neither (e.g. a value shared with "mid", or a genuine mismatch).
    val = read_yaml_scalar(yaml_path, tissue, property_name);

    tissue_mask = strcmp(params_table.tissue, tissue);
    bottom_val  = params_table.(property_name)(tissue_mask & strcmp(params_table.level, 'bottom'));
    top_val     = params_table.(property_name)(tissue_mask & strcmp(params_table.level, 'top'));

    if abs(val - bottom_val) < 1e-6
        level = 1;
    elseif abs(val - top_val) < 1e-6
        level = 2;
    else
        level = NaN;
        warning('config_bound_level:noMatch', ...
            '%s: %s.%s = %.6g matches neither bottom (%.6g) nor top (%.6g) - excluded from the additive prediction.', ...
            yaml_path, tissue, property_name, val, bottom_val, top_val);
    end
end

function val = read_yaml_scalar(yaml_path, section, key)
% Minimal flat-YAML scalar reader for config_medium_*.yaml's fixed
% structure (medium_properties: > {water,brain,skull}: > key: value, at
% 2-space then 4-space indentation). Finds the <section> block (e.g.
% "brain:") and returns the numeric value of <key> within it.
    txt   = fileread(yaml_path);
    lines = regexp(txt, '\r\n|\r|\n', 'split');
    indent  = cellfun(@(l) numel(regexp(l, '^ *', 'match', 'once')), lines);
    trimmed = strtrim(lines);

    idx_start = find(indent == 2 & strcmp(trimmed, [section ':']), 1);
    if isempty(idx_start)
        error('read_yaml_scalar:noSection', 'Section "%s:" not found in %s', section, yaml_path);
    end

    rel_end = find(indent(idx_start+1:end) == 2, 1);
    if isempty(rel_end)
        idx_end = numel(lines) + 1;
    else
        idx_end = idx_start + rel_end;
    end

    val = NaN;
    for i = (idx_start+1):(idx_end-1)
        if indent(i) == 4 && startsWith(trimmed{i}, [key ':'])
            parts = strsplit(trimmed{i}, ':');
            val = str2double(strtrim(parts{2}));
            return;
        end
    end
    error('read_yaml_scalar:noKey', 'Key "%s" not found in section "%s:" of %s', key, section, yaml_path);
end

function vals = load_combined_outcome(pn_tussim, setup_folder, sub_id, metrics)
% Reads the single-row sub-XXX_phantom.csv from a combined (all-properties-
% at-once) simulation folder, e.g. data/tussim/medium_conservative/sub-101.
    fn = fullfile(pn_tussim, setup_folder, sprintf('sub-%03d', sub_id), ...
        sprintf('sub-%03d_phantom.csv', sub_id));
    Tb = readtable(fn);
    vals = struct();
    for m = 1:numel(metrics)
        vals.(metrics{m}) = Tb.(metrics{m})(1);
    end
end

function plot_constellation_scatter(results, pn, fig_name)
% Predicted (x) vs. actual (y) combined delta, one point per benchmark x
% setup x metric. Points on the dashed 1:1 line are perfectly additive;
% above it (in the direction away from zero) means the real combined
% effect is stronger than OAT sensitivity alone would predict (synergy);
% below it means weaker (sub-additive / saturating).
    color_map = struct('riseT_brain', [0.85 0.10 0.10], 'riseT_skull', [0.90 0.50 0.10], ...
        'Isppa_brain', [0.00 0.45 0.74], 'Ipa_target', [0.30 0.30 0.30]);

    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.25 0.2 0.34 0.45]);
    hold on;

    pred_vals   = [results.predicted_delta];
    actual_vals = [results.actual_delta];
    lim = max(abs([pred_vals, actual_vals])) * 1.2;

    plot([-lim lim], [-lim lim], '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([-lim lim], [0 0], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot([0 0], [-lim lim], '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5, 'HandleVisibility', 'off');

    plotted_metrics = {};
    for i = 1:numel(results)
        r = results(i);
        marker = 'o'; if strcmp(r.setup, 'liberal'), marker = 's'; end
        scatter(r.predicted_delta, r.actual_delta, 80, color_map.(r.metric), marker, ...
            'filled', 'MarkerEdgeColor', 'k');
        if ~ismember(r.metric, plotted_metrics)
            plotted_metrics{end+1} = r.metric; %#ok<AGROW>
        end
    end

    xlim([-lim lim]); ylim([-lim lim]);
    axis square;
    xlabel('Predicted \Delta (sum of OAT single-property effects)', 'Interpreter', 'tex');
    ylabel('Actual \Delta (combined simulation - mid)', 'Interpreter', 'tex');
    title({'Constellation check: additive OAT prediction vs. actual combined simulation', ...
        'dashed = 1:1 (additive)   o = conservative   \Box = liberal'}, 'FontSize', 10, 'Interpreter', 'tex');
    set(gca, 'Layer', 'top', 'Box', 'off');

    legend_h = gobjects(1, numel(plotted_metrics));
    for i = 1:numel(plotted_metrics)
        legend_h(i) = plot(nan, nan, 'o', 'MarkerFaceColor', color_map.(plotted_metrics{i}), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    legend(legend_h, plotted_metrics, 'Location', 'eastoutside', 'Box', 'off', 'Interpreter', 'none');

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
