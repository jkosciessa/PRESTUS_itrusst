%% Tissue Property Ranges
%
% Visualizes the brain/skull tissue parameter ranges in data/tissue_parameters.csv
% (bottom/mid/top per tissue x property) as mini bar plots, in the same visual
% style as plot_metric_with_patches (used by PRESTUS_paramvar_02_benchmark2_plot.m):
%   1. One absolute-value mini bar plot per parameter (brain vs. skull bars)
%   2. One unified plot showing all parameters normalized to (a) % of default
%      (mid value) and (b) % of minimum (bottom value)

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.data    = fullfile(rootpath, 'data');
pn.figures = fullfile(rootpath, 'figures', 'tissue_ranges');
pn.code    = fullfile(rootpath, 'code'); addpath(pn.code)
if ~exist(pn.figures, 'dir'), mkdir(pn.figures); end

params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));

tissue_labels    = {'brain', 'skull'};
tissue_suffixes  = {'b', 's'};
property_labels  = {'sound_speed', 'density', 'alpha_coeff', 'specific_heat_capacity', ...
                     'thermal_conductivity', 'perfusion', 'absorption_fraction'};
property_symbols = {'c', '\rho', '\alpha', 'c_p', '\kappa', '\omega', '\beta'};
property_units   = {'m/s', 'kg/m^3', 'dB/(MHz^y cm)', 'J/(kg K)', ...
                     'W/(m K)', 'kg/(m^3 s)', 'fraction (0-1)'};

%% 1. Build a results table (bottom / mid / top per tissue x property)

results = table();
for p = 1:length(property_labels)
    prop = property_labels{p};
    for t = 1:length(tissue_labels)
        tissue = tissue_labels{t};
        tissue_mask = strcmp(params_table.tissue, tissue);

        mid_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'mid'));
        top_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'top'));
        bot_val = params_table.(prop)(tissue_mask & strcmp(params_table.level, 'bottom'));

        new_row = {tissue, tissue_suffixes{t}, prop, bot_val, mid_val, top_val};
        temp_table = cell2table(new_row, 'VariableNames', ...
            {'Tissue', 'TissueSuffix', 'Property', 'BOTTOM', 'MID', 'TOP'});
        results = [results; temp_table];
    end
end

% % of default (mid = 100 %) and % of minimum (bottom = 100 %)
results.PctDefault_bottom = 100 * results.BOTTOM ./ results.MID;
results.PctDefault_top    = 100 * results.TOP    ./ results.MID;
results.PctMin_mid        = 100 * results.MID    ./ results.BOTTOM;
results.PctMin_top        = 100 * results.TOP    ./ results.BOTTOM;

%% 2. Individual absolute-value mini bar plots (one per parameter, brain vs. skull)

for p = 1:length(property_labels)
    prop_mask = strcmp(results.Property, property_labels{p});
    row_brain = prop_mask & strcmp(results.Tissue, 'brain');
    row_skull = prop_mask & strcmp(results.Tissue, 'skull');

    % rows = [bottom; mid; top], columns = [brain, skull]
    data_matrix = [results.BOTTOM(row_brain), results.BOTTOM(row_skull); ...
                   results.MID(row_brain),    results.MID(row_skull); ...
                   results.TOP(row_brain),    results.TOP(row_skull)];

    xlabels = {[property_symbols{p} '_{b}'], [property_symbols{p} '_{s}']};

    plot_metric_with_patches(data_matrix, xlabels, property_units{p}, ...
        [property_labels{p} '_range'], [], [2 2], {[0 0 0], [0 0 0]}, pn);
end

% additional plot: same per-parameter absolute-value ranges, combined into a
% single integrated figure with one panel per property (own y-axis per panel,
% since absolute units differ across properties)
plot_all_properties_range(results, property_labels, property_symbols, property_units, ...
    pn, 'all_properties_range');

%% 3. Unified plot: % of default vs. % of minimum (all parameters x tissues)

xlabels    = cell(1, height(results));
tissue_ids = cell(1, height(results));
for r = 1:height(results)
    p = strcmp(property_labels, results.Property{r});
    xlabels{r}    = [property_symbols{p} '_{' results.TissueSuffix{r} '}'];
    tissue_ids{r} = results.TissueSuffix{r};
end

plot_unified_pct_ranges(results, xlabels, tissue_ids, pn, 'unified_pct_ranges');

% additional plot: same data, but with all brain (grey) bars grouped on the
% left and all skull (red) bars grouped on the right
tissue_order = [find(strcmp(results.TissueSuffix, 'b')); find(strcmp(results.TissueSuffix, 's'))];

plot_unified_pct_ranges(results(tissue_order, :), xlabels(tissue_order), tissue_ids(tissue_order), ...
    pn, 'unified_pct_ranges_by_tissue');

% additional plot: brain and skull split into separate panels (same colors),
% sharing a y-axis per row so the two tissues are directly comparable
idx_brain = find(strcmp(results.TissueSuffix, 'b'));
idx_skull = find(strcmp(results.TissueSuffix, 's'));

plot_unified_pct_ranges_split(results, xlabels, tissue_ids, idx_brain, idx_skull, ...
    pn, 'unified_pct_ranges_split');

%% Local functions

function plot_all_properties_range(results, property_labels, property_symbols, property_units, pn, fig_name)
    n_props = numel(property_labels);
    n_cols  = 4;
    n_rows  = ceil(n_props / n_cols);

    h = figure;
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.05 0.15 0.8 0.45 * n_rows]);

    for p = 1:n_props
        prop_mask = strcmp(results.Property, property_labels{p});
        row_brain = prop_mask & strcmp(results.Tissue, 'brain');
        row_skull = prop_mask & strcmp(results.Tissue, 'skull');

        bottom = [results.BOTTOM(row_brain), results.BOTTOM(row_skull)];
        mid    = [results.MID(row_brain),    results.MID(row_skull)];
        top    = [results.TOP(row_brain),    results.TOP(row_skull)];

        xlabels    = {[property_symbols{p} '_{b}'], [property_symbols{p} '_{s}']};
        tissue_ids = {'b', 's'};

        ax = subplot(n_rows, n_cols, p);
        draw_abs_panel(ax, bottom, top, mid, tissue_ids, xlabels, property_units{p});
        title(ax, property_labels{p}, 'Interpreter', 'none');
    end

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    save_vector_figure(h, pn.figures, fig_name);
end

function draw_abs_panel(ax, bottom, top, mid, tissue_ids, xlabels, y_label)
    hold(ax, 'on');

    barWidth  = 0.6;
    lightGrey = [0.5 0.5 0.5];
    red       = [0.7 0 0];

    for i = 1:numel(bottom)
        x = [i - barWidth/2, i + barWidth/2, i + barWidth/2, i - barWidth/2];
        y = [bottom(i), bottom(i), top(i), top(i)];
        if strcmp(tissue_ids{i}, 'b')
            patch(ax, x, y, lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
        else
            patch(ax, x, y, red, 'EdgeColor', 'none');
        end
    end

    scatter(ax, 1:numel(bottom), mid, 60, [0 0 0], 'filled', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1);

    set(ax, 'XTick', 1:numel(xlabels), 'XTickLabel', xlabels, ...
        'TickLabelInterpreter', 'tex', 'FontSize', 10);
    xtickangle(ax, 45);
    ylabel(ax, y_label);
    set(ax, 'Layer', 'top');
end

function plot_unified_pct_ranges(results, xlabels, tissue_ids, pn, fig_name)
    n = height(results);

    h = figure;
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.2 0.1 0.35 0.6]);

    ax1 = subplot(2, 1, 1);
    draw_pct_panel(ax1, results.PctDefault_bottom', results.PctDefault_top', ...
        100 * ones(1, n), tissue_ids, xlabels, '% of default (mid = 100%)', 100);
    title(ax1, 'Range relative to default (mid) value');

    ax2 = subplot(2, 1, 2);
    draw_pct_panel(ax2, 100 * ones(1, n), results.PctMin_top', ...
        results.PctMin_mid', tissue_ids, xlabels, '% of minimum (bottom = 100%)', 100);
    title(ax2, 'Range relative to minimum (bottom) value');

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    save_vector_figure(h, pn.figures, fig_name);
end

function plot_unified_pct_ranges_split(results, xlabels, tissue_ids, idx_brain, idx_skull, pn, fig_name)
    h = figure;
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.1 0.05 0.55 0.7]);

    ax1 = subplot(2, 2, 1);
    draw_pct_panel(ax1, results.PctDefault_bottom(idx_brain)', results.PctDefault_top(idx_brain)', ...
        100 * ones(1, numel(idx_brain)), tissue_ids(idx_brain), xlabels(idx_brain), ...
        '% of default (mid = 100%)', 100);
    title(ax1, 'Brain');

    ax2 = subplot(2, 2, 2);
    draw_pct_panel(ax2, results.PctDefault_bottom(idx_skull)', results.PctDefault_top(idx_skull)', ...
        100 * ones(1, numel(idx_skull)), tissue_ids(idx_skull), xlabels(idx_skull), ...
        '% of default (mid = 100%)', 100);
    title(ax2, 'Skull');
    linkaxes([ax1 ax2], 'y');

    ax3 = subplot(2, 2, 3);
    draw_pct_panel(ax3, 100 * ones(1, numel(idx_brain)), results.PctMin_top(idx_brain)', ...
        results.PctMin_mid(idx_brain)', tissue_ids(idx_brain), xlabels(idx_brain), ...
        '% of minimum (bottom = 100%)', 100);

    ax4 = subplot(2, 2, 4);
    draw_pct_panel(ax4, 100 * ones(1, numel(idx_skull)), results.PctMin_top(idx_skull)', ...
        results.PctMin_mid(idx_skull)', tissue_ids(idx_skull), xlabels(idx_skull), ...
        '% of minimum (bottom = 100%)', 100);
    linkaxes([ax3 ax4], 'y');

    h.Color = 'white';
    h.InvertHardcopy = 'off';
    save_vector_figure(h, pn.figures, fig_name);
end

function save_vector_figure(h, figdir, fig_name)
    % Force the painters (vector) renderer regardless of on-screen renderer,
    % so eps/pdf output is true vector graphics rather than a rasterized image.
    set(h, 'Renderer', 'painters');
    filepath = fullfile(figdir, fig_name);
    print(h, filepath, '-depsc', '-painters');
    print(h, filepath, '-dpdf', '-painters');
    print(h, filepath, '-dpng', '-r300');
end

function draw_pct_panel(ax, lower, upper, mid_vals, tissue_ids, xlabels, y_label, ref_value)
    hold(ax, 'on');

    barWidth  = 0.6;
    lightGrey = [0.5 0.5 0.5];
    red       = [0.7 0 0];

    for i = 1:numel(lower)
        x = [i - barWidth/2, i + barWidth/2, i + barWidth/2, i - barWidth/2];
        y = [lower(i), lower(i), upper(i), upper(i)];
        if strcmp(tissue_ids{i}, 'b')
            patch(ax, x, y, lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
        else
            patch(ax, x, y, red, 'EdgeColor', 'none');
        end
    end

    plot(ax, ax.XLim, [ref_value ref_value], '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    scatter(ax, 1:numel(lower), mid_vals, 60, [0 0 0], 'filled', ...
        'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1);

    set(ax, 'XTick', 1:numel(xlabels), 'XTickLabel', xlabels, ...
        'TickLabelInterpreter', 'tex', 'FontSize', 12);
    xtickangle(ax, 45);
    ylabel(ax, y_label);
    set(ax, 'Layer', 'top');
end
