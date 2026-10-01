% PRESTUS_PARAMVAR_10_MORRIS_ANALYSIS_PLOT
%
% Analyses the Morris elementary-effects design built and (partially or
% fully) run by PRESTUS_paramvar_09_morris_interactions_sim.m. Joins the
% trajectory design table (figures/paramvar/morris_design.csv) against
% each step's result CSV (data/tussim/<variant_name>/sub-002/sub-002_
% phantom<output_affix>.csv) to get, per design-table row:
%   Isppa_brain, Ipa_target, riseT_brain, riseT_skull
%
% For every trajectory and every step>0 (i.e. every flipped factor), the
% ELEMENTARY EFFECT is outcome(step) - outcome(step-1), sign-flipped to
% always read as "effect of bottom -> top" (see PRESTUS_paramvar_09's
% header for the full method description). Per factor (tissue x
% property), across all trajectories with both endpoints available:
%   mu     = mean(elementary effects)          - signed average effect
%   mu_star = mean(|elementary effects|)        - overall importance
%   sigma  = std(elementary effects)            - interaction / context-
%            dependence: large sigma means the factor's effect depends on
%            where the other 13 factors currently are
%   n      = number of trajectories contributing a usable elementary
%            effect for this factor (missing/not-yet-run evaluations on
%            either side of a step are dropped, not imputed)
%
% This script tolerates a PARTIALLY completed design (e.g. only the
% acoustic-cost phase-1 evaluations done so far): any step whose result
% CSV is missing is simply excluded from mu*/sigma for that step's
% factor, and completeness is reported explicitly so a small n is never
% silently mistaken for a stable estimate.
%
% Output:
%   figures/paramvar/morris_effects.csv   (mu, mu_star, sigma, n per
%                                           factor per outcome)
%   figures/paramvar/morris_plot_<outcome>.{png,eps}  (mu* vs sigma
%                                           scatter, one point per factor)

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

variant_name = 'itrusst_protocol1a_morris';   % must match PRESTUS_paramvar_09
subject_id   = 2;

pn.results = fullfile(rootpath, 'data', 'tussim', variant_name, sprintf('sub-%03d', subject_id));
pn.figures = fullfile(rootpath, 'figures', 'paramvar');
if ~isfolder(pn.figures); mkdir(pn.figures); end

outcomes = {'Isppa_brain', 'Ipa_target', 'riseT_brain', 'riseT_skull'};

%% ── Load the trajectory design table ─────────────────────────────────────

% Delimiter forced to ',' - readtable's auto-detection can misfire on
% underscore-saturated headers/values (see PRESTUS_paramvar_13's fix for
% the concrete failure this caused on a similarly-shaped design table).
design = readtable(fullfile(pn.figures, 'morris_design.csv'), 'Delimiter', ',', 'TextType', 'string');
n_rows = height(design);

for o = 1:numel(outcomes)
    design.(outcomes{o}) = NaN(n_rows, 1);
end
design.available = false(n_rows, 1);

%% ── Join each design row against its result CSV (tolerating gaps) ─────────

n_missing = 0;   % no result CSV at all (not yet run / still queued)
n_partial = 0;   % result CSV exists but lacks the thermal columns (heating
                  % stage did not finish - e.g. killed by the HPC walltime
                  % limit partway through; acoustic-only fields were still
                  % written, so the file exists but is incomplete)
missing_affixes = strings(0, 1);
partial_affixes = strings(0, 1);
for i = 1:n_rows
    affix = design.output_affix(i);
    result_file = fullfile(pn.results, sprintf('sub-%03d_phantom%s.csv', subject_id, affix));
    if ~isfile(result_file)
        n_missing = n_missing + 1;
        missing_affixes(end+1) = affix; %#ok<AGROW>
        continue
    end
    R = readtable(result_file, 'Delimiter', ',');
    if ~all(ismember(outcomes, R.Properties.VariableNames))
        n_partial = n_partial + 1;
        partial_affixes(end+1) = affix; %#ok<AGROW>
        continue
    end
    design.available(i) = true;
    for o = 1:numel(outcomes)
        design.(outcomes{o})(i) = R.(outcomes{o})(1);
    end
end

fprintf('Morris design: %d/%d evaluations have complete results.\n', n_rows - n_missing - n_partial, n_rows);
if n_missing > 0
    fprintf('  %d not yet run / no result file at all.\n', n_missing);
end
if n_partial > 0
    fprintf(['  %d ran the acoustic stage but never wrote thermal output (riseT_*) - likely killed \n' ...
             '  mid-heating-stage (check hpc.timelimit vs. actual thermal runtime). Affected output_affixes:\n'], n_partial);
    for k = 1:numel(partial_affixes)
        fprintf('    %s\n', partial_affixes(k));
    end
end
if n_missing > 0 || n_partial > 0
    fprintf('All incomplete evaluations above are excluded from mu*/sigma below (not imputed).\n');
end

%% ── Compute elementary effects per trajectory/step, per factor ────────────

trajectories = unique(design.trajectory);
% One struct array entry per (trajectory, flipped-factor) pair; besides
% tissue/property/trajectory it carries one field per outcome, holding
% that outcome's signed elementary effect (NaN if either endpoint of the
% step is missing/incomplete).
ee_fields = [{'tissue', 'property', 'trajectory'}, outcomes];
ee_args = cell(1, 2 * numel(ee_fields));
ee_args(1:2:end) = ee_fields;
ee_args(2:2:end) = {{}};
ee = struct(ee_args{:});

for t = 1:numel(trajectories)
    traj_rows = find(design.trajectory == trajectories(t));
    [~, order] = sort(design.step(traj_rows));
    traj_rows = traj_rows(order); % ensure step-ascending within this trajectory

    for k = 2:numel(traj_rows)
        row_prev = traj_rows(k - 1);
        row_cur  = traj_rows(k);

        entry.tissue     = design.tissue(row_cur);
        entry.property   = design.property(row_cur);
        entry.trajectory = trajectories(t);
        for o = 1:numel(outcomes)
            if design.available(row_prev) && design.available(row_cur)
                d = design.(outcomes{o})(row_cur) - design.(outcomes{o})(row_prev);
                if design.direction(row_cur) == "top_to_bottom"
                    d = -d;
                end
            else
                d = NaN;
            end
            entry.(outcomes{o}) = d;
        end
        ee(end+1) = entry; %#ok<AGROW>
    end
end

%% ── Aggregate elementary effects into mu / mu_star / sigma per factor ─────

factor_keys = unique([[ee.tissue]', [ee.property]'], 'rows', 'stable');
n_factors = size(factor_keys, 1);

effects = table('Size', [n_factors, 2 + 4 * numel(outcomes)], ...
    'VariableTypes', [{'string', 'string'}, repmat({'double'}, 1, 4 * numel(outcomes))], ...
    'VariableNames', [{'tissue', 'property'}, reshape(strcat(repmat(outcomes, 4, 1), ...
        repmat({'_mu'; '_mu_star'; '_sigma'; '_n'}, 1, numel(outcomes)))', 1, [])]);

for f = 1:n_factors
    tissue_f = factor_keys(f, 1);
    property_f = factor_keys(f, 2);
    mask = [ee.tissue]' == tissue_f & [ee.property]' == property_f;

    effects.tissue(f) = tissue_f;
    effects.property(f) = property_f;

    for o = 1:numel(outcomes)
        vals = [ee(mask).(outcomes{o})]';
        vals = vals(~isnan(vals));
        effects.(strcat(outcomes{o}, '_mu'))(f)      = mean(vals, 'omitnan');
        effects.(strcat(outcomes{o}, '_mu_star'))(f) = mean(abs(vals), 'omitnan');
        if numel(vals) >= 2
            effects.(strcat(outcomes{o}, '_sigma'))(f) = std(vals);
        else
            effects.(strcat(outcomes{o}, '_sigma'))(f) = NaN; % std undefined for n<2
        end
        effects.(strcat(outcomes{o}, '_n'))(f)       = numel(vals);
    end
end

writetable(effects, fullfile(pn.figures, 'morris_effects.csv'));

fprintf('\n%-8s %-24s | %8s %8s %8s %3s | %8s %8s %8s %3s\n', ...
    'Tissue', 'Property', 'muISPPA', 'muISPPA*', 'sigISPPA', 'n', 'muRiseTb', 'muRiseTb*', 'sigRiseTb', 'n');
for f = 1:n_factors
    fprintf('%-8s %-24s | %+8.3f %8.3f %8.3f %3d | %+8.4f %8.4f %8.4f %3d\n', ...
        effects.tissue(f), effects.property(f), ...
        effects.Isppa_brain_mu(f), effects.Isppa_brain_mu_star(f), effects.Isppa_brain_sigma(f), effects.Isppa_brain_n(f), ...
        effects.riseT_brain_mu(f), effects.riseT_brain_mu_star(f), effects.riseT_brain_sigma(f), effects.riseT_brain_n(f));
end

%% ── Morris plots: mu* (importance) vs sigma (interaction), per outcome ────

for o = 1:numel(outcomes)
    plot_morris(effects, outcomes{o}, pn, sprintf('morris_plot_%s', outcomes{o}));
end

%% ── Local functions ──────────────────────────────────────────────────────

function plot_morris(effects, outcome, pn, fig_name)
% Scatter of mu* (x, overall importance) vs sigma (y, interaction /
% context-dependence), one point per factor, labelled by tissue_property.
% Factors with n=0 (no usable elementary effect yet) are shown greyed out
% at the origin so incompleteness is visible rather than silently absent.
    mu_star = effects.(strcat(outcome, '_mu_star'));
    sigma   = effects.(strcat(outcome, '_sigma'));
    n       = effects.(strcat(outcome, '_n'));
    labels  = strcat(effects.tissue, {' '}, effects.property);

    has_data = n > 0;
    % NaN mu*/sigma (no usable elementary effect yet) are shown as grey
    % markers at the origin rather than left un-rendered, so incompleteness
    % is visible in the plot rather than silently absent.
    mu_star_plot = mu_star; mu_star_plot(~has_data) = 0;
    sigma_plot   = sigma;   sigma_plot(~has_data | isnan(sigma)) = 0;

    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.2 0.2 0.5 0.6]);
    hold on;

    scatter(mu_star_plot(has_data), sigma_plot(has_data), 60, 'filled', 'MarkerFaceColor', [0.85 0.10 0.10]);
    scatter(mu_star_plot(~has_data), sigma_plot(~has_data), 40, 'MarkerEdgeColor', [0.6 0.6 0.6]);

    text(mu_star_plot + range([0; mu_star_plot]) * 0.01, sigma_plot, labels, 'FontSize', 8, 'Interpreter', 'none');

    xlabel(sprintf('\\mu* (%s)  -  overall importance', strrep(outcome, '_', '\_')), 'Interpreter', 'tex');
    ylabel('\sigma  -  interaction / context-dependence', 'Interpreter', 'tex');
    title(sprintf('Morris screening: %s  (grey = no usable elementary effect yet)', strrep(outcome, '_', '\_')), ...
        'FontSize', 10, 'Interpreter', 'tex');
    box off;

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
