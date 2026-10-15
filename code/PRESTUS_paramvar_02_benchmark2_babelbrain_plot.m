% PRESTUS_PARAMVAR_02_BENCHMARK2_BABELBRAIN_PLOT
%
% Differences between the BabelBrain runs of benchmark 2 (Sam Pichardo,
% data/pichardo, run at 6 and 9 points per wavelength) and PRESTUS (sub-002,
% itrusst_protocol1a) for the one-at-a-time parameter variation:
%
%   sub-002_babelbrain_ppw_difference   BabelBrain 9 vs 6 PPW, relative
%                                       difference per condition and metric
%   sub-002_babelbrain_reference        reference-medium outcomes per tool
%   sub-002_babelbrain_sensitivity      top - bottom level per condition, in
%                                       % of the tool's reference outcome
%   sub-002_babelbrain_ranges           bottom/top outcomes per condition
%                                       and tool, with reference lines
%
% BabelBrain reports the temperature rise at the target only, so PRESTUS
% focal heating is read at the fixed nominal target (riseT_target, see
% load_focus_heating_metrics) here.
%
% Output: figures/paramvar/babelbrain/

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data     = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.pichardo = itrusst_input_path(fullfile(rootpath, 'data', 'pichardo'), ...
    'ITRUSST_BABELBRAIN_DIR', 'BabelBrain benchmark 2 results (data/pichardo)');
pn.figures  = fullfile(rootpath, 'figures', 'paramvar', 'babelbrain');
if ~isfolder(pn.figures); mkdir(pn.figures); end

%% settings

subject_id        = 2;
tissue_labels     = {'brain', 'skull'};
tissue_suffixes   = {'_b', '_s'};
property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
sweep_values      = 1:2; % bottom, top
ppw               = [6 9];

% metrics shared by both tools
metrics = {'Isppa_brain',  'Brain I_{SPPA} (W/cm^2)'
           'riseT_brain',  'Brain temperature rise (°C)'
           'riseT_skull',  'Skull temperature rise (°C)'
           'riseT_target', 'Target temperature rise (°C)'};

S = paper_style();

%% load

bb = struct('res', {}, 'ref', {});
for k = 1:numel(ppw)
    [bb(k).res, bb(k).ref, xlabels] = load_babelbrain_benchmark2(pn.pichardo, ppw(k), ...
        tissue_labels, property_labels, property_symbols);
end

prestus = load_paramvar_results(pn.data, subject_id, tissue_labels, tissue_suffixes, ...
    property_labels, property_suffixes, sweep_values, property_symbols);
focus = load_focus_heating_metrics(pn.data, subject_id, tissue_labels, tissue_suffixes, ...
    property_labels, property_suffixes, sweep_values, property_symbols);
prestus.riseT_target = focus.riseT_target;
prestus_ref = paper_reference_metrics(pn.data, subject_id);

tools       = [{'PRESTUS'}, arrayfun(@(p) sprintf('BabelBrain %d PPW', p), ppw, 'UniformOutput', false)];
tool_res    = [{prestus}, {bb.res}];
tool_ref    = [{prestus_ref}, {bb.ref}];
tool_colors = [S.PRESTUS; 0.5 * S.BabelBrain + 0.5; S.BabelBrain];
n_conds     = numel(xlabels);

%% BabelBrain 9 vs 6 PPW: relative difference

bb_fields = fieldnames(bb(1).res);
bb_names  = {'I_{SPPA,brain}', 'MI', '\DeltaT_{brain}', '\DeltaT_{skull}', '\DeltaT_{target}', ...
             'CEM43_{brain}', 'CEM43_{skull}'};
rel = nan(2 * n_conds + 1, numel(bb_fields));
for f = 1:numel(bb_fields)
    v6 = [bb(1).res.(bb_fields{f})(:); bb(1).ref.(bb_fields{f})];
    v9 = [bb(2).res.(bb_fields{f})(:); bb(2).ref.(bb_fields{f})];
    rel(:, f) = 100 * (v9 - v6) ./ v6;
end
row_labels = [reshape([strcat(xlabels, ' bottom'); strcat(xlabels, ' top')], [], 1); {'reference'}];

h = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 22], 'InvertHardcopy', 'off');
ax = axes(h);
imagesc(ax, rel);
c_max = max(abs(rel(:)), [], 'omitnan');
colormap(ax, diverging_cmap(256));
clim(ax, [-c_max c_max]);
for r = 1:size(rel, 1)
    for f = 1:size(rel, 2)
        text(ax, f, r, sprintf('%+.1f', rel(r, f)), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end
set(ax, 'XTick', 1:numel(bb_names), 'XTickLabel', bb_names, 'YTick', 1:numel(row_labels), ...
    'YTickLabel', row_labels, 'TickLabelInterpreter', 'tex', 'TickLength', [0 0], 'FontSize', 8, ...
    'XAxisLocation', 'top');
yline(ax, n_conds + 0.5, 'k-', 'LineWidth', 1); % brain | skull
yline(ax, 2 * n_conds + 0.5, 'k-', 'LineWidth', 1); % reference
cb = colorbar(ax);
ylabel(cb, sprintf('BabelBrain %d vs %d PPW: (%d - %d) / %d PPW (%%)', ppw(2), ppw(1), ppw(2), ppw(1), ppw(1)));
save_fig(h, pn.figures, 'sub-002_babelbrain_ppw_difference');

%% reference medium per tool

h = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 20 7], 'InvertHardcopy', 'off');
tl = tiledlayout(h, 1, size(metrics, 1), 'TileSpacing', 'compact');
for m = 1:size(metrics, 1)
    ax = nexttile(tl);
    vals = cellfun(@(r) r.(metrics{m, 1}), tool_ref);
    b = bar(ax, vals, 0.6, 'FaceColor', 'flat', 'EdgeColor', 'none');
    b.CData = tool_colors;
    set(ax, 'XTickLabel', tools, 'XTickLabelRotation', 45, 'Box', 'off', 'FontSize', 9);
    ylabel(ax, metrics{m, 2}, 'Interpreter', 'tex');
end
title(tl, 'Reference medium (all properties at mid)', 'FontSize', 10);
save_fig(h, pn.figures, 'sub-002_babelbrain_reference');

%% sensitivity: top - bottom level relative to the reference outcome

h = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 24], 'InvertHardcopy', 'off');
tl = tiledlayout(h, size(metrics, 1), 1, 'TileSpacing', 'compact');
for m = 1:size(metrics, 1)
    ax = nexttile(tl);
    eff = nan(n_conds, numel(tools));
    for k = 1:numel(tools)
        r = tool_res{k}.(metrics{m, 1});
        eff(:, k) = 100 * (r(2, :) - r(1, :))' / tool_ref{k}.(metrics{m, 1});
    end
    b = bar(ax, eff, 'grouped', 'EdgeColor', 'none');
    for k = 1:numel(tools); b(k).FaceColor = tool_colors(k, :); end
    xline(ax, n_conds / 2 + 0.5, 'k:'); % brain | skull
    set(ax, 'XTick', 1:n_conds, 'XTickLabel', xlabels, 'TickLabelInterpreter', 'tex', ...
        'Box', 'off', 'FontSize', 9, 'XLim', [0.4 n_conds + 0.6]);
    ylim(ax, 'padded');
    ylabel(ax, {metrics{m, 2}, 'top - bottom (% of ref.)'}, 'Interpreter', 'tex');
    if m == 1; legend(ax, tools, 'Location', 'southwest', 'Box', 'off'); end
end
save_fig(h, pn.figures, 'sub-002_babelbrain_sensitivity');

%% bottom/top outcomes per condition and tool

offsets = linspace(-0.25, 0.25, numel(tools));
h = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 24], 'InvertHardcopy', 'off');
tl = tiledlayout(h, size(metrics, 1), 1, 'TileSpacing', 'compact');
for m = 1:size(metrics, 1)
    ax = nexttile(tl);
    hold(ax, 'on');
    hl = gobjects(1, numel(tools));
    for k = 1:numel(tools)
        r = tool_res{k}.(metrics{m, 1});
        yline(ax, tool_ref{k}.(metrics{m, 1}), '--', 'Color', tool_colors(k, :), 'LineWidth', 0.75);
        x = (1:n_conds) + offsets(k);
        hk = plot(ax, [x; x], r, '-', 'Color', tool_colors(k, :), 'LineWidth', 3);
        hl(k) = hk(1);
        scatter(ax, x, r(1, :), 14, S.low,  'filled', 'MarkerEdgeColor', tool_colors(k, :));
        scatter(ax, x, r(2, :), 14, S.high, 'filled', 'MarkerEdgeColor', tool_colors(k, :));
    end
    xline(ax, n_conds / 2 + 0.5, 'k:'); % brain | skull
    set(ax, 'XTick', 1:n_conds, 'XTickLabel', xlabels, 'TickLabelInterpreter', 'tex', ...
        'Box', 'off', 'FontSize', 9, 'XLim', [0.4 n_conds + 0.6]);
    ylim(ax, [0, ax.YLim(2)]);
    ylabel(ax, metrics{m, 2}, 'Interpreter', 'tex');
    if m == 1; legend(ax, hl, tools, 'Location', 'best', 'Box', 'off'); end
end
title(tl, 'bottom (white) / top (black) level; dashed: reference medium', 'FontSize', 10);
save_fig(h, pn.figures, 'sub-002_babelbrain_ranges');

%% summary

fprintf('\nReference medium, sub-002 / benchmark 2:\n');
fprintf('%-16s', ''); fprintf('%18s', tools{:}); fprintf('\n');
for m = 1:size(metrics, 1)
    fprintf('%-16s', metrics{m, 1});
    fprintf('%18.3f', cellfun(@(r) r.(metrics{m, 1}), tool_ref));
    fprintf('\n');
end
fprintf('\nBabelBrain %d vs %d PPW, max. |relative difference| over conditions:\n', ppw(2), ppw(1));
for f = 1:numel(bb_fields)
    fprintf('  %-14s %5.1f %%\n', bb_fields{f}, max(abs(rel(:, f)), [], 'omitnan'));
end

%% helpers

function save_fig(h, out_dir, name)
saveas(h, fullfile(out_dir, name), 'epsc');
saveas(h, fullfile(out_dir, name), 'png');
fprintf('  saved %s\n', fullfile(out_dir, name));
end

function cmap = diverging_cmap(n)
% blue - white - red
lo = [0.13 0.40 0.67]; mid = [1 1 1]; hi = [0.70 0.09 0.17];
t = linspace(0, 1, n / 2)';
cmap = [lo + t .* (mid - lo); mid + t .* (hi - mid)];
end
