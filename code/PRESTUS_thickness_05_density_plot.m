% PRESTUS_THICKNESS_05_DENSITY_PLOT
%
% Plots the skull density x thickness design of
% PRESTUS_thickness_04_density_sim.m: each outcome (Isppa_brain, Ipa_target,
% riseT_brain, riseT_skull) vs. skull thickness, one line per density level.
% Parallel lines = no interaction; a widening/narrowing gap = interaction.
% Missing or incomplete result files are reported and left as gaps.
%
% Output: figures/paramvar/density_thickness_<outcome>.{png,eps}

%% ── Paths ───────────────────────────────────────────────────────────────────

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'));
rootpath = pwd;

variant_name = 'itrusst_protocol1a_density_thickness';   % must match PRESTUS_thickness_04
pn.results = fullfile(rootpath, 'data', 'tussim', variant_name);
pn.figures = fullfile(rootpath, 'figures', 'paramvar');

outcomes = {'Isppa_brain', 'Ipa_target', 'riseT_brain', 'riseT_skull'};
density_levels_name = {'bottom', 'mid', 'top'};
level_colors = [0.00 0.45 0.74; 0.50 0.50 0.50; 0.85 0.10 0.10];   % blue / grey / red

%% ── Load the design table and join results ─────────────────────────────────

% force ',' delimiter: auto-detection can misfire on underscore-heavy values
design = readtable(fullfile(pn.figures, 'density_thickness_design.csv'), 'Delimiter', ',', 'TextType', 'string');
n_rows = height(design);

for o = 1:numel(outcomes)
    design.(outcomes{o}) = NaN(n_rows, 1);
end
design.available = false(n_rows, 1);

n_missing = 0;
for i = 1:n_rows
    result_file = fullfile(pn.results, sprintf('sub-%03d', design.subject_id(i)), ...
        sprintf('sub-%03d_phantom%s.csv', design.subject_id(i), design.output_affix(i)));
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

fprintf('Density x thickness: %d/%d evaluations available (%d missing/incomplete).\n', ...
    n_rows - n_missing, n_rows, n_missing);

%% ── Plot: one figure per outcome, one line per density level ───────────────

for o = 1:numel(outcomes)
    plot_density_thickness(design, outcomes{o}, density_levels_name, level_colors, pn, ...
        sprintf('density_thickness_%s', outcomes{o}));
end

%% ── Local functions ──────────────────────────────────────────────────────

function plot_density_thickness(design, outcome, density_levels_name, level_colors, pn, fig_name)
    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.25 0.25 0.35 0.4]);
    hold on;

    legend_handles = gobjects(1, numel(density_levels_name));
    for d = 1:numel(density_levels_name)
        mask = design.skull_density_level == density_levels_name{d} & design.available;
        sub = sortrows(design(mask, :), 'thickness_mm');
        legend_handles(d) = plot(sub.thickness_mm, sub.(outcome), '-o', ...
            'Color', level_colors(d, :), 'MarkerFaceColor', level_colors(d, :), 'LineWidth', 2);
    end

    grid on;
    xlabel('Skull thickness (mm)');
    ylabel(strrep(outcome, '_', '\_'), 'Interpreter', 'tex');
    xlim([1 12]);
    legend(legend_handles, strcat('skull density: ', density_levels_name), 'Location', 'best');
    legend('boxoff');
    title(sprintf('%s vs. skull thickness, by skull density level', strrep(outcome, '_', '\_')), ...
        'FontSize', 10, 'Interpreter', 'tex');
    set(gca, 'FontSize', 11, 'Layer', 'top', 'Box', 'off');

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
