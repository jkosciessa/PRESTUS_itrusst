% PAPER_FIG5_SETUPS
%
% Paper Fig. 5.
%   left:  parallel coordinates of the skull and brain properties of the
%          heating setups Temp-H / Temp-L (heat_max / heat_min: per property the
%          bound that maximises / minimises brain heating; black) and the
%          emulated toolbox defaults (PRESTUS, k-Plan*, BabelBrain*), and the
%          acoustic-only transmission setup Trans-L (grey, dashed; no thermal
%          properties, so its line ends at attenuation). Trans-H / Trans-L
%          are listed in full in paper_supp_table_setups.m. Axes are flipped where Temp-H uses the bottom bound, so
%          that Temp-H runs along the top and Temp-L along the bottom. Axis ends show
%          the min (blue) / max (red) value across setups. Data as in
%          PRESTUS_itrusst_setups_plot.m (attenuation at 500 kHz in dB/cm;
%          k-Plan perfusion is not reported and bridged dashed). Absorption:
%          skull as coefficient at 500 kHz in dB/cm, brain as the absorbed
%          fraction of attenuation (beta/alpha) in %, the quantity the brain
%          bounds are set on (absorption_fraction in tissue_setups.csv). Skull panel above the brain panel; grey band and labels
%          for the acoustic, red for the thermal properties (the bar colours).
%   right: maximum brain intensity (Isppa_brain), bracketed by Trans-H /
%          Trans-L, and maximum temperature rise, bracketed by Temp-H / Temp-L,
%          for the 3D skull (sub-102) and no-skull (sub-101) benchmarks, as in
%          PRESTUS_pipeline_02_benchmark1_plot.m. Not saved while a setup
%          result is missing.
%
% Output: figures/paper/fig5/fig5_parallel.{pdf,png}
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
property_fields = {'sound_speed', 'density', 'alpha_coeff', 'absorption_coeff', ...
                   'thermal_conductivity', 'specific_heat_capacity', 'perfusion'};
property_labels = {'Sound speed', 'Density', 'Attenuation', 'Absorption', ...
                   'Thermal conduct.', 'Specific heat', 'Perfusion'};
i_att = strcmp(property_fields, 'alpha_coeff');
i_abs = strcmp(property_fields, 'absorption_coeff');

setups = {'heat_max', 'heat_min', 'trans_min'}; % tissue_setups.csv columns
rows   = {'Temp-H', 'Temp-L', 'Trans-L', 'PRESTUS', 'k-Plan*', 'BabelBrain*'};
trans_clr = [0.45 0.45 0.45];
colors = {S.setup, S.setup, trans_clr, S.PRESTUS, S.kPlan, S.BabelBrain};
styles = {'-', '-', '--', '-', '-', '-'};
is_setup = [true(1, numel(setups)), false(1, 3)];
media  = cellfun(@(s) itrusst_tissue_properties(s, pn.data), setups, 'UniformOutput', false);
for tool = {'PRESTUS', 'kPlan', 'BabelBrain'}
    cfg = yaml.loadFile(fullfile(pn.configs, ['config_medium_' tool{1} '.yaml']));
    media{end+1} = cfg.medium_properties; %#ok<SAGROW>
end
setups_table = readtable(fullfile(pn.data, 'tissue_setups.csv'), 'TextType', 'char');

h = paper_figure([9.5 13.0]);
tissues = {'skull', 'brain'};
for t = 1:numel(tissues)
    tissue = tissues(t);
    data = nan(numel(media), numel(property_fields));
    for s = 1:numel(media)
        m = media{s}.(tissue{1});
        for p = 1:numel(property_fields)
            data(s, p) = m.(property_fields{p});
        end
        data(s, i_att) = m.alpha_coeff * freq_mhz ^ m.alpha_power;
        if strcmp(tissue{1}, 'brain')
            data(s, i_abs) = 100 * m.absorption_coeff / m.alpha_coeff; % % of attenuation (same power law)
        else
            data(s, i_abs) = m.absorption_coeff * freq_mhz ^ m.alpha_power; % dB/cm at 500 kHz
        end
    end
    data(strcmp(rows, 'k-Plan*'), strcmp(property_fields, 'perfusion')) = NaN; % not reported

    flip_cols = false(1, numel(property_fields));
    for p = 1:numel(property_fields)
        r = strcmp(setups_table.tissue, tissue{1}) & strcmp(setups_table.property, property_fields{p});
        if ~any(r) && p == find(i_abs)
            % brain absorption is set as absorption_fraction (x alpha_coeff)
            r = strcmp(setups_table.tissue, tissue{1}) & strcmp(setups_table.property, 'absorption_fraction');
        end
        flip_cols(p) = strcmp(setups_table.heat_max{r}, 'bottom');
        for s = 1:numel(setups) % not specified in tissue_setups.csv (acoustic-only setups)
            if isempty(setups_table.(setups{s}){r}); data(s, p) = NaN; end
        end
    end

    ax = axes(h, 'Units', 'centimeters', 'Position', [1.6, 7.6 - (t-1) * 5.6, 6.0, 4.4]);
    opts = struct('thermal', ~ismember(property_fields, {'sound_speed', 'density', 'alpha_coeff'}), ...
                  'show_xlabels', strcmp(tissue{1}, 'brain'), 'italic_ends', true, 'connect', false, ...
                  'end_suffix', {repmat({''}, 1, numel(property_fields))});
    if strcmp(tissue{1}, 'brain'); opts.end_suffix{i_abs} = '%'; end
    paper_parallel_panel(ax, data, flip_cols, rows, colors, styles, is_setup, property_labels, opts);
end
annotation(h, 'line', [0.08 0.92], [0.555 0.555], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
paper_save(h, pn.figures, 'fig5_parallel');
close(h);

%% toolbox comparison (3D benchmarks)

% each panel is bracketed by the setup pair built for its outcome
tool_dirs = {'medium_kPlan', 'medium_PRESTUS', 'medium_BabelBrain'};
setup_dirs.Isppa = [{'medium_trans_max'}, tool_dirs, {'medium_trans_min'}];
setup_dirs.riseT = [{'medium_heat_max'},  tool_dirs, {'medium_heat_min'}];
pair_labels.Isppa = {'Trans-H', 'Trans-L'};
pair_labels.riseT = {'Temp-H', 'Temp-L'};
tick_labels = {'', 'k-Plan*', 'PRESTUS', 'BabelBrain*', ''};
benchmarks = struct('sub', {'sub-102', 'sub-101'}, 'name', {'skull', 'noskull'}, ...
    'title', {'Skull benchmark', 'No skull benchmark'});

vals = struct();
for b = 1:numel(benchmarks)
    vals(b).Isppa = nan(1, numel(setup_dirs.Isppa));
    vals(b).riseT = nan(1, numel(setup_dirs.riseT));
    for s = 1:numel(setup_dirs.Isppa)
        f = fullfile(pn.tussim, setup_dirs.Isppa{s}, benchmarks(b).sub, [benchmarks(b).sub '_phantom.csv']);
        if ~isfile(f); warning('File %s not found.', f); continue; end
        T = readtable(f);
        vals(b).Isppa(s) = T.Isppa_brain(1);
    end
    for s = 1:numel(setup_dirs.riseT)
        f = fullfile(pn.tussim, setup_dirs.riseT{s}, benchmarks(b).sub, [benchmarks(b).sub '_phantom.csv']);
        if ~isfile(f); warning('File %s not found.', f); continue; end
        T = readtable(f);
        if ismember('maxT', T.Properties.VariableNames)
            vals(b).riseT(s) = T.maxT(1) - 37;
        else
            warning('File %s has no maxT (thermal simulation missing).', f);
        end
    end
end
riseT_lim = [0, max(0.5, ceil(2 * max([vals.riseT], [], 'omitnan')) / 2)];

for b = 1:numel(benchmarks)
    if any(isnan([vals(b).Isppa, vals(b).riseT]))
        warning('Bars for the %s benchmark not saved: setup results missing.', benchmarks(b).name);
        continue
    end
    h = paper_figure([8.0 5.8]);
    panels = {vals(b).Isppa, 'Brain I_{SPPA} (W/cm^2)', S.bar_int, [], pair_labels.Isppa; ...
              vals(b).riseT, 'Temperature rise (°C)', S.bar_temp, riseT_lim, pair_labels.riseT};
    for k = 1:2
        [v, ylab, clr, ylims, pair] = panels{k, :};
        ax = axes(h, 'Units', 'centimeters', 'Position', [1.2 + (k-1) * 3.8, 1.6, 2.6, 3.6]);
        hold(ax, 'on');
        bh = bar(ax, 1:numel(v), v, 0.8, 'FaceColor', 'flat', 'EdgeColor', 'none');
        bh.CData = repmat(clr, numel(v), 1);
        bh.CData([1 end], :) = repmat(S.setup, 2, 1);
        if isempty(ylims); ylims = [0, max(5, 5 * ceil(max(v, [], 'omitnan') / 5))]; end
        labels = cellfun(@(t) sprintf('\\color[rgb]{%.2f %.2f %.2f}\\it%s', clr, t), tick_labels, 'UniformOutput', false);
        labels([1 end]) = cellfun(@(t) ['\bf ' t], pair, 'UniformOutput', false);
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
