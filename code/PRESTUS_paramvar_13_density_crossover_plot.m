% PRESTUS_PARAMVAR_13_DENSITY_CROSSOVER_PLOT
%
% Analyses the design built and run by PRESTUS_paramvar_08_constellation_
% factorial_sim.m (skull density x each of the other 13 tissue-property
% factors, 2x2 grid per pair, everything else at midpoint).
%
% Joins figures/paramvar/density_crossover_design.csv against each row's
% result CSV (data/tussim/itrusst_protocol1a_density_crossover/sub-002/
% sub-002_phantom<output_affix>.csv), then for every "other" factor
% computes the exact 2x2 interaction:
%   effect_at_other_bottom = outcome(density=top,other=bottom) - outcome(density=bottom,other=bottom)
%   effect_at_other_top    = outcome(density=top,other=top)    - outcome(density=bottom,other=top)
%   interaction             = effect_at_other_top - effect_at_other_bottom
% interaction ~ 0 means density's effect does not depend on that factor;
% a large interaction means it does (and effect_at_other_bottom/_top show
% the two effect sizes directly, not just their difference).
%
% acoustic_anchor rows (density-only, everything else mid) are loaded for
% completeness reporting but are not part of any pair's 2x2 grid, so they
% do not feed into the interaction calculation.
%
% Tolerates a partially completed run (missing/incomplete result files are
% reported and leave that cell/pair as NaN, not imputed).
%
% Output:
%   figures/paramvar/density_crossover_effects.csv (effect_at_other_bottom,
%   effect_at_other_top, interaction, n_available per factor per outcome)
%   figures/paramvar/density_crossover_<outcome>.{png,eps} - one row per
%   factor, showing density's effect on the outcome when that factor is at
%   bottom (blue dot) vs. top (red dot), connected by a line: horizontal
%   distance from 0 is the main effect of density, and the gap between the
%   two dots is the interaction (the quantity the tornado plot used to show
%   in isolation, now shown alongside the effect sizes it is a difference
%   of). Rows sorted by mean |effect| (largest overall density effect at
%   top), so importance and interaction are both visible in the same view.

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

variant_name = 'itrusst_protocol1a_density_crossover';   % must match PRESTUS_paramvar_08
subject_id   = 2;

pn.results = fullfile(rootpath, 'data', 'tussim', variant_name, sprintf('sub-%03d', subject_id));
pn.figures = fullfile(rootpath, 'figures', 'paramvar');

outcomes = {'Isppa_brain', 'Ipa_target', 'riseT_brain', 'riseT_skull'};

%% ── Load the design table and join results ─────────────────────────────────

% Delimiter forced to ',' - readtable's auto-detection can misfire on this
% file, since both the headers (other_tissue, other_property, ...) and the
% output_affix values are saturated with underscores, which it can mistake
% for an additional delimiter alongside the comma.
design = readtable(fullfile(pn.figures, 'density_crossover_design.csv'), 'Delimiter', ',', 'TextType', 'string');
n_rows = height(design);

for o = 1:numel(outcomes)
    design.(outcomes{o}) = NaN(n_rows, 1);
end
design.available = false(n_rows, 1);

n_missing = 0;
for i = 1:n_rows
    result_file = fullfile(pn.results, sprintf('sub-%03d_phantom%s.csv', subject_id, design.output_affix(i)));
    if ~isfile(result_file)
        n_missing = n_missing + 1;
        continue
    end
    R = readtable(result_file, 'Delimiter', ',');
    if ~all(ismember(outcomes, R.Properties.VariableNames))
        n_missing = n_missing + 1;   % ran but incomplete (e.g. killed mid-thermal-stage)
        continue
    end
    design.available(i) = true;
    for o = 1:numel(outcomes)
        design.(outcomes{o})(i) = R.(outcomes{o})(1);
    end
end

fprintf('Density crossover: %d/%d evaluations available (%d missing/incomplete).\n', ...
    n_rows - n_missing, n_rows, n_missing);

%% ── Compute the 2x2 interaction per "other" factor ─────────────────────────

pair_mask = design.eval_type ~= "acoustic_anchor";
pairs = design(pair_mask, :);
factor_keys = unique([pairs.other_tissue, pairs.other_property], 'rows', 'stable');
n_factors = size(factor_keys, 1);

effects = table('Size', [n_factors, 2 + 4 * numel(outcomes)], ...
    'VariableTypes', [{'string', 'string'}, repmat({'double'}, 1, 4 * numel(outcomes))], ...
    'VariableNames', [{'tissue', 'property'}, reshape(strcat(repmat(outcomes, 4, 1), ...
        repmat({'_effect_at_bottom'; '_effect_at_top'; '_interaction'; '_n'}, 1, numel(outcomes)))', 1, [])]);

for f = 1:n_factors
    tissue_f   = factor_keys(f, 1);
    property_f = factor_keys(f, 2);
    mask = pairs.other_tissue == tissue_f & pairs.other_property == property_f;
    rows = pairs(mask, :);

    effects.tissue(f) = tissue_f;
    effects.property(f) = property_f;

    for o = 1:numel(outcomes)
        v = get_cell(rows, outcomes{o});   % v(focal_level, other_level), NaN if missing/unavailable

        effect_at_bottom = v(2, 1) - v(1, 1);
        effect_at_top     = v(2, 2) - v(1, 2);
        interaction       = effect_at_top - effect_at_bottom;
        n_available = sum(rows.available);

        effects.(strcat(outcomes{o}, '_effect_at_bottom'))(f) = effect_at_bottom;
        effects.(strcat(outcomes{o}, '_effect_at_top'))(f)     = effect_at_top;
        effects.(strcat(outcomes{o}, '_interaction'))(f)       = interaction;
        effects.(strcat(outcomes{o}, '_n'))(f)                 = n_available;
    end
end

writetable(effects, fullfile(pn.figures, 'density_crossover_effects.csv'));

fprintf('\n%-8s %-24s | %10s %10s\n', 'Tissue', 'Property', 'interaction(Isppa)', 'interaction(riseT_brain)');
for f = 1:n_factors
    fprintf('%-8s %-24s | %+18.4f %+18.6f\n', effects.tissue(f), effects.property(f), ...
        effects.Isppa_brain_interaction(f), effects.riseT_brain_interaction(f));
end

%% ── Overview plots: main effect + interaction together, one per outcome ────

for o = 1:numel(outcomes)
    plot_density_effect_overview(effects, outcomes{o}, pn, sprintf('density_crossover_%s', outcomes{o}));
end

%% ── Local functions ──────────────────────────────────────────────────────

function v = get_cell(rows, outcome)
% v(focal_level, other_level) for focal_level,other_level in {1,2}; NaN if
% that specific combination is missing from rows or its evaluation failed.
    v = NaN(2, 2);
    for r = 1:height(rows)
        if ~rows.available(r), continue; end
        fl = rows.focal_level(r);
        ol = rows.other_level(r);
        v(fl, ol) = rows.(outcome)(r);
    end
end

function plot_density_effect_overview(effects, outcome, pn, fig_name)
% Dumbbell plot: one row per factor, two dots (density's effect on
% `outcome` when the paired factor is at bottom [blue] vs. top [red]),
% connected by a line. Distance from x=0 is the main effect of density;
% the horizontal gap between the two dots is the interaction. Rows with
% either effect missing (incomplete evaluation) are dropped, not imputed.
    effect_bottom = effects.(strcat(outcome, '_effect_at_bottom'));
    effect_top     = effects.(strcat(outcome, '_effect_at_top'));
    labels_all = strcat(effects.tissue, {' '}, effects.property);

    keep = ~isnan(effect_bottom) & ~isnan(effect_top);
    effect_bottom = effect_bottom(keep);
    effect_top     = effect_top(keep);
    labels = labels_all(keep);

    if isempty(effect_bottom)
        fprintf('No complete factor pairs for %s yet - skipping plot.\n', outcome);
        return
    end

    mean_abs_effect = (abs(effect_bottom) + abs(effect_top)) / 2;
    [~, order] = sort(mean_abs_effect, 'ascend'); % ascend so largest overall effect plots at top
    effect_bottom = effect_bottom(order);
    effect_top     = effect_top(order);
    labels = labels(order);

    n = numel(effect_bottom);
    color_bottom = [0.00 0.45 0.74]; % blue = other factor at its bottom bound
    color_top    = [0.85 0.10 0.10]; % red  = other factor at its top bound

    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.2 0.1 0.4 0.6]);
    hold on;

    x_all = [effect_bottom; effect_top; 0];
    xpad = max(range(x_all) * 0.08, eps);

    for i = 1:n
        plot([effect_bottom(i) effect_top(i)], [i i], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    end
    scatter(effect_bottom, 1:n, 55, 'filled', 'MarkerFaceColor', color_bottom);
    scatter(effect_top, 1:n, 55, 'filled', 'MarkerFaceColor', color_top);

    plot([0 0], [0.3 n + 0.7], '--k', 'LineWidth', 1);
    % 'none' interpreter: labels are plain tissue/property names (e.g.
    % "sound_speed") - a 'tex' interpreter would render the underscore as a
    % subscript command instead of a literal character.
    set(gca, 'YTick', 1:n, 'YTickLabel', labels, 'TickLabelInterpreter', 'none', 'FontSize', 10);
    ylim([0.3 n + 0.7]);
    xlim([min(x_all) - xpad, max(x_all) + xpad]);
    xlabel(sprintf('skull density effect on %s (bottom \\rightarrow top)', strrep(outcome, '_', '\_')), ...
        'Interpreter', 'tex');
    set(gca, 'Layer', 'top', 'Box', 'off');
    title({sprintf('%s: dot = density''s effect with the paired factor at bottom (blue) / top (red)', ...
        strrep(outcome, '_', '\_')), 'gap between dots = interaction; distance from 0 = main effect'}, ...
        'FontSize', 9, 'Interpreter', 'tex');

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
