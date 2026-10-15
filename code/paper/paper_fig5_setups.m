% PAPER_FIG5_SETUPS
%
% Paper Fig. 5.
%   left:  parallel coordinates of the skull and brain properties of the
%          joint setups (Setup A = setup1, Setup B = setup2; black) and the
%          emulated toolbox defaults (PRESTUS, k-Plan*, BabelBrain*). Axes
%          are flipped where Setup A uses the bottom bound, so that Setup A
%          runs along the top and Setup B along the bottom. Axis ends show
%          the min (blue) / max (red) value across setups. Data as in
%          PRESTUS_itrusst_setups_plot.m (attenuation at 500 kHz in dB/cm;
%          k-Plan perfusion is not reported and bridged dashed).
%   right: maximum brain intensity (Isppa_brain) and maximum temperature rise for
%          the 3D skull (sub-102) and no-skull (sub-101) benchmarks, as in
%          PRESTUS_pipeline_02_benchmark1_plot.m.
%
% Output: figures/paper/fig5/fig5_parallel_{skull,brain}.{pdf,png}
%         figures/paper/fig5/fig5_bars_{skull,noskull}.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath); % PRESTUS yaml reader
pn.tussim  = fullfile(rootpath, 'data', 'tussim');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'fig5');

S = paper_style();

%% setup properties

freq_mhz = 0.5;
property_fields = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_labels = {'Sound speed', 'Density', 'Attenuation', 'Thermal conduct.', ...
                   'Specific heat', 'Perfusion', 'Absorption'};

rows   = {'Setup A', 'Setup B', 'PRESTUS', 'k-Plan*', 'BabelBrain*'};
colors = {S.setup, S.setup, S.PRESTUS, S.kPlan, S.BabelBrain};
media  = {itrusst_tissue_properties('setup1', pn.data), itrusst_tissue_properties('setup2', pn.data)};
for tool = {'PRESTUS', 'kPlan', 'BabelBrain'}
    cfg = yaml.loadFile(fullfile(pn.configs, ['config_medium_' tool{1} '.yaml']));
    media{end+1} = cfg.medium_properties; %#ok<SAGROW>
end
setups_table = readtable(fullfile(pn.data, 'tissue_setups.csv'), 'TextType', 'char');

for tissue = {'skull', 'brain'}
    data = nan(numel(media), numel(property_fields));
    for s = 1:numel(media)
        m = media{s}.(tissue{1});
        for p = 1:numel(property_fields)
            data(s, p) = m.(property_fields{p});
        end
        data(s, 3) = m.alpha_coeff * freq_mhz ^ m.alpha_power;
    end
    data(strcmp(rows, 'k-Plan*'), strcmp(property_fields, 'perfusion')) = NaN; % not reported

    flip_cols = false(1, numel(property_fields));
    for p = 1:numel(property_fields)
        r = strcmp(setups_table.tissue, tissue{1}) & strcmp(setups_table.property, property_fields{p});
        flip_cols(p) = strcmp(setups_table.setup1{r}, 'bottom');
    end

    h = paper_figure([9.5 7.5]);
    ax = axes(h, 'Units', 'centimeters', 'Position', [1.6 1.9 6.0 4.6]);
    parallel_panel(ax, data, flip_cols, rows, colors, property_labels, S);
    paper_save(h, pn.figures, ['fig5_parallel_' tissue{1}]);
    close(h);
end

%% toolbox comparison (3D benchmarks)

setup_dirs = {'medium_setup1', 'medium_kPlan', 'medium_PRESTUS', 'medium_BabelBrain', 'medium_setup2'};
tick_labels = {'', 'k-Plan*', 'PRESTUS', 'BabelBrain*', ''};
benchmarks = struct('sub', {'sub-102', 'sub-101'}, 'name', {'skull', 'noskull'}, ...
    'title', {'Skull benchmark', 'No skull benchmark'});

vals = struct();
for b = 1:numel(benchmarks)
    vals(b).Isppa = nan(1, numel(setup_dirs));
    vals(b).riseT = nan(1, numel(setup_dirs));
    for s = 1:numel(setup_dirs)
        f = fullfile(pn.tussim, setup_dirs{s}, benchmarks(b).sub, [benchmarks(b).sub '_phantom.csv']);
        if ~isfile(f); warning('File %s not found.', f); continue; end
        T = readtable(f);
        vals(b).Isppa(s) = T.Isppa_brain(1);
        if ismember('maxT', T.Properties.VariableNames)
            vals(b).riseT(s) = T.maxT(1) - 37;
        else
            warning('File %s has no maxT (thermal simulation missing).', f);
        end
    end
end
riseT_lim = [0, max(0.5, ceil(2 * max([vals.riseT], [], 'omitnan')) / 2)];

for b = 1:numel(benchmarks)
    h = paper_figure([8.0 5.8]);
    panels = {vals(b).Isppa, 'Brain I_{SPPA} (W/cm^2)', S.bar_int, []; ...
              vals(b).riseT, 'Temperature rise (°C)', S.bar_temp, riseT_lim};
    for k = 1:2
        [v, ylab, clr, ylims] = panels{k, :};
        ax = axes(h, 'Units', 'centimeters', 'Position', [1.2 + (k-1) * 3.8, 1.6, 2.6, 3.6]);
        hold(ax, 'on');
        bh = bar(ax, 1:numel(v), v, 0.8, 'FaceColor', 'flat', 'EdgeColor', 'none');
        bh.CData = repmat(clr, numel(v), 1);
        bh.CData([1 end], :) = repmat(S.setup, 2, 1);
        if isempty(ylims); ylims = [0, max(5, 5 * ceil(max(v, [], 'omitnan') / 5))]; end
        text(ax, [1 numel(v)], [0 0] + 0.02 * diff(ylims), {'A', 'B'}, 'Color', 'w', ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        labels = cellfun(@(t) sprintf('\\color[rgb]{%.2f %.2f %.2f}\\it%s', clr, t), tick_labels, 'UniformOutput', false);
        set(ax, 'XTick', 1:numel(v), 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', 'XLim', [0.4 numel(v) + 0.6], ...
            'YLim', ylims, 'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, ...
            'TickLength', [0 0.02], 'Box', 'on');
        xtickangle(ax, 30);
        ylabel(ax, ylab, 'Interpreter', 'tex');
    end
    annotation(h, 'textbox', [0 0.92 1 0.08], 'String', benchmarks(b).title, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontName', S.font, 'FontSize', S.fontsize + 1, 'FontWeight', 'bold');
    paper_save(h, pn.figures, ['fig5_bars_' benchmarks(b).name]);
    close(h);
end

%% local functions

function parallel_panel(ax, data, flip_cols, rows, colors, property_labels, S)
% parallel coordinates in the layout of the paper figure: per-property
% min-max normalisation (cf. make_parallel_plot), setups labelled at the
% left, toolboxes at the right, min/max values at the axis ends
n_props = size(data, 2);
margin  = 0.08;
lo = min(data, [], 1, 'omitnan'); hi = max(data, [], 1, 'omitnan');
rng = hi - lo; rng(rng == 0) = 1;
y = (data - lo) ./ rng;
y(:, flip_cols) = 1 - y(:, flip_cols);
y = y * (1 - 2 * margin) + margin;

hold(ax, 'on');
for p = 1:n_props
    plot(ax, [p p], [0 1], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);
end
jitter = linspace(-0.12, 0.12, size(data, 1));
for s = 1:size(data, 1)
    x = (1:n_props) + jitter(s);
    plot(ax, x, y(s, :), '-', 'Color', colors{s}, 'LineWidth', 1.2);
    ok = ~isnan(y(s, :));
    for g = find(ok(1:end-1) & ~ok(2:end)) % bridge missing values
        e = g + find(ok(g+1:end), 1);
        plot(ax, x([g e]), y(s, [g e]), '--', 'Color', colors{s}, 'LineWidth', 1.2);
    end
    plot(ax, x(ok), y(s, ok), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colors{s}, 'MarkerEdgeColor', 'none');
    if startsWith(rows{s}, 'Setup')
        text(ax, 0.75, y(s, 1), rows{s}, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    end
end
% toolbox labels at the right end, spread apart where they would overlap
tools = find(~startsWith(rows, 'Setup'));
y_end = arrayfun(@(s) y(s, find(~isnan(y(s, :)), 1, 'last')), tools);
[y_sorted, order] = sort(y_end);
for k = 2:numel(y_sorted)
    y_sorted(k) = max(y_sorted(k), y_sorted(k-1) + 0.08);
end
for k = 1:numel(order)
    s = tools(order(k));
    text(ax, n_props + 0.35, y_sorted(k), rows{s}, 'Color', colors{s}, 'HorizontalAlignment', 'left');
end
min_clr = [0 0 0.85]; max_clr = [0.85 0 0];
for p = 1:n_props
    if flip_cols(p)
        b = {hi(p), max_clr}; t = {lo(p), min_clr};
    else
        b = {lo(p), min_clr}; t = {hi(p), max_clr};
    end
    text(ax, p, -0.05, fmtnum(b{1}), 'Color', b{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(ax, p, 1.05, fmtnum(t{1}), 'Color', t{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
set(ax, 'XLim', [0.6 n_props + 0.4], 'YLim', [-0.15 1.15], 'XTick', 1:n_props, ...
    'XTickLabel', property_labels, 'YColor', 'none', 'TickLength', [0 0]);
xtickangle(ax, 30);
end

function s = fmtnum(v)
% at most 3 significant digits after the point, no scientific notation
if v == round(v)
    s = sprintf('%d', round(v));
else
    s = regexprep(regexprep(sprintf('%.3f', v), '0+$', ''), '\.$', '');
end
end
