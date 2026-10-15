% PAPER_FIG1_PROPERTY_RANGES
%
% Paper Fig. 1 (right): literature ranges of brain (grey) and skull (red)
% tissue properties from data/tissue_parameters.csv, one axis per property
% (bottom-top bar, reference 'mid' value as black dot). Attenuation is
% shown at the protocol frequency (500 kHz) in dB/cm, absorption fraction
% in %.
%
% The benchmark schematic (Fig. 1, left) is drawn in Illustrator.
%
% Output: figures/paper/fig1/fig1_property_ranges.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));
pn.data    = fullfile(rootpath, 'data');
pn.figures = fullfile(rootpath, 'figures', 'paper', 'fig1');

%% properties (order and axis limits of the paper panel)

freq_mhz = 0.5;

props = struct( ...
    'field', {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
              'specific_heat_capacity', 'perfusion', 'absorption_fraction'}, ...
    'label', {'sound speed', 'density', 'attenuation', {'thermal', 'conductivity'}, ...
              {'specific heat', 'capacity'}, 'perfusion', {'absorption', 'fraction'}}, ...
    'unit',  {'m/s', 'kg/m^3', 'dB/cm', 'W/(m K)', 'J/(kg K)', 'mL/min/kg', '%'}, ...
    'ylim',  {[1500 3000], [1000 2000], [0 12], [0.25 0.6], [500 4000], [0 900], [0 100]});

T = readtable(fullfile(pn.data, 'tissue_parameters.csv'), 'TextType', 'char');
tissues = {'brain', 'skull'};
levels  = {'bottom', 'mid', 'top'};

%% plot

S = paper_style();
ax_w = 1.1; ax_h = 4.6; ax_gap = 1.35; left0 = 1.0; bottom0 = 1.4;
h = paper_figure([left0 + numel(props) * (ax_w + ax_gap), bottom0 + ax_h + 1.0]);

for p = 1:numel(props)
    ax = axes(h, 'Units', 'centimeters', 'Position', [left0 + (p-1)*(ax_w + ax_gap), bottom0, ax_w, ax_h]);
    hold(ax, 'on');
    vals = nan(numel(tissues), 3); % bottom, mid, top
    for t = 1:numel(tissues)
        rows = strcmp(T.tissue, tissues{t});
        for l = 1:3
            row = rows & strcmp(T.level, levels{l});
            v = T.(props(p).field)(row);
            if strcmp(props(p).field, 'alpha_coeff')
                v = v * freq_mhz ^ T.alpha_power(row);
            elseif strcmp(props(p).field, 'absorption_fraction')
                v = 100 * v;
            end
            vals(t, l) = v;
        end
    end

    clrs = {S.brain, S.skull};
    for t = 1:numel(tissues)
        lo = vals(t, 1); hi = vals(t, 3);
        if hi - lo < 0.004 * diff(props(p).ylim) % keep zero-width ranges visible
            plot(ax, t + [-0.3 0.3], [lo lo], '-', 'Color', clrs{t}, 'LineWidth', 1.5);
        else
            patch(ax, t + [-0.3 0.3 0.3 -0.3], [lo lo hi hi], clrs{t}, 'EdgeColor', 'none');
        end
        plot(ax, t, vals(t, 2), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    end

    yl = props(p).ylim;
    yl = [min(yl(1), min(vals(:))), max(yl(2), max(vals(:)))];
    if ~isequal(yl, props(p).ylim)
        warning('%s: values exceed the panel limits, axis extended.', props(p).field);
    end
    yt = linspace(yl(1), yl(2), 6);
    ytl = repmat({''}, 1, numel(yt));
    ytl([1 end]) = {num2str(yl(1)), num2str(yl(2))};
    set(ax, 'XLim', [0.4 2.6], 'YLim', yl, 'YTick', yt, 'YTickLabel', ytl, ...
        'XTick', [], 'XColor', 'none', 'TickLength', [0.06 0.06]);
    title(ax, props(p).unit, 'Interpreter', 'tex', 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3]);
    text(ax, 1.5, yl(1) - 0.06 * diff(yl), props(p).label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', S.fontsize + 1);
    if p == 1
        text(ax, 2, vals(2, 1), 'skull', 'Color', S.skull, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        text(ax, 1.5, vals(1, 3), '  brain', 'Color', S.brain, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    end
end

paper_save(h, pn.figures, 'fig1_property_ranges');
close(h);
