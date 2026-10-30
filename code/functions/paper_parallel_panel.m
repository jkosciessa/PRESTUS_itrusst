function paper_parallel_panel(ax, data, flip_cols, rows, colors, styles, is_setup, property_labels, opts)
% PAPER_PARALLEL_PANEL
%
% Parallel coordinates in the layout of the paper figures (Fig. 5, Supp.
% Fig. 4): per-property min-max normalisation (cf. make_parallel_plot),
% setups (is_setup) labelled at the left, toolboxes at the right, min (blue)
% / max (red) values at the axis ends. flip_cols reverses an axis.
%
% Inputs:
%   ax              - target axes
%   data            - rows x properties, NaN = missing (bridged dashed)
%   flip_cols       - logical, per property
%   rows, colors, styles - per row: label, RGB, line style
%   is_setup        - logical, per row
%   property_labels - per property
%   opts            - (optional) struct with fields
%                     .thermal      logical per property: draws a band below
%                                   the axes (acoustic grey / thermal red, the
%                                   Fig. 5 bar colours) and colours the axis
%                                   labels the same way
%                     .show_xlabels false to hide the axis labels (upper
%                                   panel of a stacked figure); default true
%                     .italic_ends  true for italic axis-end values
%                     .end_suffix   cell per property, appended to the
%                                   axis-end values (e.g. '%'); default ''
%                     .connect      false for markers only, without lines
%                                   between the properties; default true
if nargin < 9; opts = struct(); end
if ~isfield(opts, 'show_xlabels'); opts.show_xlabels = true; end
if ~isfield(opts, 'italic_ends'); opts.italic_ends = false; end
if ~isfield(opts, 'connect'); opts.connect = true; end
if ~isfield(opts, 'end_suffix'); opts.end_suffix = repmat({''}, 1, size(data, 2)); end
S = paper_style();
n_props = size(data, 2);
margin  = 0.08;
lo = min(data, [], 1, 'omitnan'); hi = max(data, [], 1, 'omitnan');
rng = hi - lo; rng(rng == 0) = 1;
y = (data - lo) ./ rng;
y(:, flip_cols) = 1 - y(:, flip_cols);
y = y * (1 - 2 * margin) + margin;

hold(ax, 'on');
if isfield(opts, 'thermal')
    % bands below the axes, one per run of acoustic / thermal properties
    band_clr = {1 - 0.35 * (1 - S.bar_int), 1 - 0.45 * (1 - S.bar_temp)}; % light tints
    edges = [0, find(diff(opts.thermal)), n_props];
    for g = 1:numel(edges) - 1
        x0 = edges(g) + 1; x1 = edges(g + 1);
        patch(ax, [x0 - 0.2, x1 + 0.2, x1 + 0.2, x0 - 0.2], [0 0 margin margin] - 0.01, ...
            band_clr{opts.thermal(x0) + 1}, 'EdgeColor', 'none');
    end
end
for p = 1:n_props
    plot(ax, [p p], [0 1], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);
end
jitter = linspace(-0.12, 0.12, size(data, 1));
for s = 1:size(data, 1)
    x = (1:n_props) + jitter(s);
    ok = ~isnan(y(s, :));
    if opts.connect
        plot(ax, x, y(s, :), styles{s}, 'Color', colors{s}, 'LineWidth', 1.2);
        for g = find(ok(1:end-1) & ~ok(2:end)) % bridge missing values
            e = g + find(ok(g+1:end), 1);
            plot(ax, x([g e]), y(s, [g e]), '--', 'Color', colors{s}, 'LineWidth', 1.2);
        end
        plot(ax, x(ok), y(s, ok), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colors{s}, 'MarkerEdgeColor', 'none');
    else
        plot(ax, x(ok), y(s, ok), 'o', 'MarkerSize', 4.5, 'MarkerFaceColor', colors{s}, 'MarkerEdgeColor', 'none');
    end
end
% setup labels at the left end, toolbox labels at the right end, spread
% apart where they would overlap
setups = find(is_setup);
y_start = spread_labels(y(setups, 1));
for k = 1:numel(setups)
    text(ax, 0.75, y_start(k), rows{setups(k)}, 'Color', colors{setups(k)}, ...
        'HorizontalAlignment', 'right', 'FontWeight', 'bold');
end
tools = find(~is_setup);
y_end = spread_labels(arrayfun(@(s) y(s, find(~isnan(y(s, :)), 1, 'last')), tools));
for k = 1:numel(tools)
    s = tools(k);
    text(ax, n_props + 0.35, y_end(k), rows{s}, 'Color', colors{s}, 'HorizontalAlignment', 'left');
end
min_clr = [0 0 0.85]; max_clr = [0.85 0 0];
for p = 1:n_props
    if flip_cols(p)
        b = {hi(p), max_clr}; t = {lo(p), min_clr};
    else
        b = {lo(p), min_clr}; t = {hi(p), max_clr};
    end
    angle = 'normal'; if opts.italic_ends; angle = 'italic'; end
    text(ax, p, -0.05, [fmtnum(b{1}) opts.end_suffix{p}], 'Color', b{2}, 'FontAngle', angle, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(ax, p, 1.05, [fmtnum(t{1}) opts.end_suffix{p}], 'Color', t{2}, 'FontAngle', angle, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
labels = property_labels;
if isfield(opts, 'thermal')
    label_clr = {S.bar_int * 0.8, S.bar_temp};
    labels = arrayfun(@(p) sprintf('\\color[rgb]{%.2f %.2f %.2f}%s', label_clr{opts.thermal(p) + 1}, ...
        property_labels{p}), 1:n_props, 'UniformOutput', false);
end
if ~opts.show_xlabels; labels = repmat({''}, 1, n_props); end
set(ax, 'XLim', [0.6 n_props + 0.4], 'YLim', [-0.15 1.15], 'XTick', 1:n_props, ...
    'XTickLabel', labels, 'TickLabelInterpreter', 'tex', 'YColor', 'none', 'TickLength', [0 0]);
xtickangle(ax, 30);
end

function y_out = spread_labels(y_in)
% label positions in the order of y_in, moved apart to at least 0.08
[y_sorted, order] = sort(y_in(:)');
for k = 2:numel(y_sorted)
    y_sorted(k) = max(y_sorted(k), y_sorted(k-1) + 0.08);
end
y_out(order) = y_sorted;
end

function s = fmtnum(v)
% at most 3 significant digits after the point, no scientific notation
if v == round(v)
    s = sprintf('%d', round(v));
else
    s = regexprep(regexprep(sprintf('%.3f', v), '0+$', ''), '\.$', '');
end
end
