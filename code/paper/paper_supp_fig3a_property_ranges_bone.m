% PAPER_SUPP_FIG3A_PROPERTY_RANGES_BONE
%
% Supplementary Fig. 3a: skull property ranges split by bone layer. Literature
% values come from the ITRUSST acoustic-properties table
% (itrusst_properties/literature_web/data/measurements.csv), using
% only rows with status 'included'. The table is not part of this
% repository; setenv('ITRUSST_MEASUREMENTS_CSV', '<file>') points to
% another copy. Per property and tissue class the bar spans the min-max of all values (incl. reported value_max), and each
% measurement is drawn as a dot (reported ranges as a vertical line).
%
% Skull is shown pooled over all divisions (as in Fig. 1) and split as
%   unclear    : whole, frontal, parietal, temporal, occipital (layer not given)
%   cortical   : cortical, cortical (outer table), cortical (inner table)
%   trabecular : cancellous, diploë
% Brain is not shown: the table gives brain absorption as a coefficient,
% whereas the simulations set it as a fraction of attenuation.
% Attenuation and absorption coefficient are tabulated at 0.5 MHz (dB/cm).
%
% Axis limits are fitted to the bone values; legend at the upper left corner.
%
% Output: figures/paper/supp_fig3/suppfig3a_property_ranges_bone_only.{pdf,png}

%% paths

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..', '..'));
rootpath = pwd;

addpath(fullfile(rootpath, 'code', 'functions'));

%% properties (order and axis limits as in paper_fig1_property_ranges)

props = struct( ...
    'field', {'sound_speed', 'density', 'attenuation', 'absorption_coefficient', ...
              'thermal_conductivity', 'specific_heat', 'perfusion'}, ...
    'label', {'sound speed', 'density', 'attenuation', {'absorption', 'coefficient'}, ...
              {'thermal', 'conductivity'}, {'specific heat', 'capacity'}, 'perfusion'}, ...
    'unit',  {'m/s', 'kg/m^3', 'dB/cm', 'dB/cm', 'W/(m K)', 'J/(kg K)', 'mL/min/kg'}, ...
    'scale', {1, 1, 1, 1, 1, 1, 1}, ...
    'ylim',  {[1500 3000], [1000 2000], [0 12], [0 4.03], [0.25 0.6], [500 4000], [0 900]}, ...
    'ylim_bone', {[2000 3000], [1000 2000], [0 12], [0 4.03], [0.25 0.4], [500 3000], [0 40]});

lit_file = itrusst_input_path(fullfile(rootpath, '..', 'itrusst_properties', 'literature_web', 'data', 'measurements.csv'), ...
    'ITRUSST_MEASUREMENTS_CSV', 'ITRUSST acoustic-properties table (measurements.csv)');
T = readtable(lit_file, 'TextType', 'char', 'Delimiter', ',');
T = T(strcmp(T.status, 'included'), :);

%% tissue classes

S = paper_style();
cls_all = struct( ...
    'name',  {'brain', 'skull', 'unclear', 'cortical', 'trabecular'}, ...
    'tissue',{'brain', 'skull', 'skull', 'skull', 'skull'}, ...
    'div',   {{}, {}, {'whole', 'frontal', 'parietal', 'temporal', 'occipital'}, ...
              {'cortical', 'cortical (outer table)', 'cortical (inner table)'}, ...
              {'cancellous', 'diploë'}}, ...
    'color', {S.brain, S.skull, [0.85 0.62 0.62], [0.42 0.00 0.00], [0.93 0.40 0.15]}, ...
    'label', {'brain', 'skull (all)', 'skull (layer unclear)', 'skull cortical', 'skull trabecular'});

V = struct( ...
    'name',     {'suppfig3a_property_ranges_bone_only'}, ...
    'dir',      {fullfile(rootpath, 'figures', 'paper', 'supp_fig3')}, ...
    'brain',    {false}, ...
    'ylim',     {'ylim_bone'}, ...
    'legend_x', {0.15});   % legend start (cm) at the figure corner

for iv = 1:numel(V)
cls = cls_all;
if ~V(iv).brain
    cls = cls_all(~strcmp({cls_all.name}, 'brain'));
end

%% plot

n_cls = numel(cls);
ax_w = 2.3; ax_h = 4.6; ax_gap = 0.9; left0 = 1.0; bottom0 = 1.4;
h = paper_figure([left0 + numel(props) * (ax_w + ax_gap), bottom0 + ax_h + 1.4]);
bw = 0.34; % half bar width

for p = 1:numel(props)
    ax = axes(h, 'Units', 'centimeters', 'Position', [left0 + (p-1)*(ax_w + ax_gap), bottom0, ax_w, ax_h]);
    hold(ax, 'on');
    yl = props(p).(V(iv).ylim);
    for c = 1:n_cls
        rows = strcmp(T.tissue, cls(c).tissue) & strcmp(T.property, props(p).field);
        if ~isempty(cls(c).div)
            rows = rows & ismember(T.division, cls(c).div);
        end
        v  = props(p).scale * T.value(rows);
        vx = props(p).scale * T.value_max(rows);
        if isempty(v)
            text(ax, c, yl(1) + 0.02 * diff(yl), 'n/a', 'Rotation', 90, 'FontSize', S.fontsize - 1, ...
                'Color', [0.6 0.6 0.6], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            continue
        end
        lo = min(v); hi = max([v; vx(~isnan(vx))]);
        if hi - lo < 0.004 * diff(yl) % keep zero-width ranges visible
            plot(ax, c + [-bw bw], [lo lo], '-', 'Color', cls(c).color, 'LineWidth', 1.5);
        else
            patch(ax, c + [-bw bw bw -bw], [lo lo hi hi], cls(c).color, 'EdgeColor', 'none');
        end
        % individual measurements, spread horizontally when values coincide
        [v, ord] = sort(v); vx = vx(ord);
        xo = zeros(size(v));
        for i = 2:numel(v)
            near = abs(v(1:i-1) - v(i)) < 0.012 * diff(yl);
            k = mod(sum(near), 7); % 0, +d, -d, +2d, ... kept within the bar
            xo(i) = 0.09 * ceil(k / 2) * (-1)^k;
        end
        for i = 1:numel(v)
            if ~isnan(vx(i))
                plot(ax, c + xo(i) + [0 0], [v(i) vx(i)], '-', 'Color', 'k', 'LineWidth', 0.75);
            end
        end
        plot(ax, c + xo, v, 'o', 'MarkerSize', 2.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        if lo < yl(1) || hi > yl(2)
            warning('%s (%s): values exceed the panel limits, axis extended.', props(p).field, cls(c).name);
            yl = [min(yl(1), lo), max(yl(2), hi)];
        end
    end

    % ticks follow the nominal limits; the axis only grows to show outliers
    yt = linspace(props(p).(V(iv).ylim)(1), props(p).(V(iv).ylim)(2), 6);
    ytl = repmat({''}, 1, numel(yt));
    ytl([1 end]) = {num2str(yt(1)), num2str(yt(end))};
    set(ax, 'XLim', [0.4 n_cls + 0.6], 'YLim', yl, 'YTick', yt, 'YTickLabel', ytl, ...
        'XTick', [], 'XColor', 'none', 'TickLength', [0.04 0.04]);
    title(ax, props(p).unit, 'Interpreter', 'tex', 'FontAngle', 'italic', 'Color', [0.3 0.3 0.3]);
    text(ax, (n_cls + 1) / 2, yl(1) - 0.06 * diff(yl), props(p).label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontSize', S.fontsize + 1);
end

% legend (top row)
lx = V(iv).legend_x; ly = bottom0 + ax_h + 0.9;
ann = axes(h, 'Units', 'centimeters', 'Position', [0 0 h.Position(3:4)], 'Visible', 'off', ...
    'XLim', [0 h.Position(3)], 'YLim', [0 h.Position(4)]);
hold(ann, 'on');
for c = 1:n_cls
    patch(ann, lx + [0 0.35 0.35 0], ly + [-0.12 -0.12 0.12 0.12], cls(c).color, 'EdgeColor', 'none');
    th = text(ann, lx + 0.5, ly, cls(c).label, 'Color', cls(c).color, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle', 'Units', 'data');
    lx = th.Extent(1) + th.Extent(3) + 0.5;
end
plot(ann, lx + 0.2, ly, 'o', 'MarkerSize', 2.2, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
text(ann, lx + 0.4, ly, 'literature value', 'VerticalAlignment', 'middle');

paper_save(h, V(iv).dir, V(iv).name);
close(h);
end
