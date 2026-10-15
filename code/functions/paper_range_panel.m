function ax = paper_range_panel(ax, data, symbols, tissues, y_label, y_limits, refs)
% Range panel of the one-at-a-time parameter variation (paper Figs. 2-4):
% per condition a bar from bottom- to top-level outcome (brain grey, skull
% red) with the bottom (white) and top (black) level as dots, and dashed
% horizontal reference lines behind the bars.
%
%   data     [2 x n]  row 1 = bottom level, row 2 = top level
%   symbols  {1 x n}  TeX property symbols, e.g. '\rho'
%   tissues  {1 x n}  'b' or 's' (tick labels are coloured by tissue)
%   y_limits [lo hi]  or [] for [0, rounded max]
%   refs     struct array with fields value, color, (label): reference
%            lines; NaN values are skipped
%
% Same encoding as plot_metric_with_patches, at print size.

S = paper_style();
n = size(data, 2);
bar_w = 0.6;
hold(ax, 'on');

if isempty(y_limits)
    y_limits = [0, nice_ceil(max(data(:), [], 'omitnan') * 1.1)];
end

% reference lines first, so they sit behind the bars
if nargin < 7; refs = struct('value', {}, 'color', {}); end
for r = 1:numel(refs)
    if isnan(refs(r).value); continue; end
    plot(ax, [0.4 n+0.6], refs(r).value * [1 1], '--', 'Color', refs(r).color, 'LineWidth', 0.75);
    if isfield(refs, 'label') && ~isempty(refs(r).label)
        text(ax, 0.5, refs(r).value, refs(r).label, 'Color', refs(r).color, ...
            'VerticalAlignment', 'bottom', 'FontSize', S.fontsize - 1);
    end
end

labels = cell(1, n);
for i = 1:n
    clr = S.brain; if strcmp(tissues{i}, 's'); clr = S.skull; end
    lo = min(data(:, i)); hi = max(data(:, i));
    if ~any(isnan([lo hi]))
        patch(ax, i + bar_w/2 * [-1 1 1 -1], [lo lo hi hi], clr, 'EdgeColor', 'none');
    end
    sym = symbols{i}; sub = tissues{i};
    if contains(sym, '_') % e.g. c_p -> c_{p,s}
        parts = strsplit(sym, '_'); sym = parts{1}; sub = [parts{2} ',' sub];
    end
    labels{i} = sprintf('\\color[rgb]{%.2f %.2f %.2f}%s_{%s}', clr, sym, sub);
end
for i = 1:n
    scatter(ax, i, data(1, i), 22, S.low,  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    scatter(ax, i, data(2, i), 22, S.high, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

set(ax, 'XTick', 1:n, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', ...
    'XLim', [0.4 n+0.6], 'YLim', y_limits, 'Layer', 'top', 'TickLength', [0 0.02], ...
    'XTickLabelRotation', 0);
ylabel(ax, y_label, 'Interpreter', 'tex');
end

function v = nice_ceil(x)
% round up to 1, 2 or 5 x 10^k steps (axis limit)
if ~(x > 0); v = 1; return; end
p = 10^floor(log10(x));
steps = [1 1.5 2 2.5 3 4 5 6 8 10];
v = p * steps(find(steps * p >= x, 1));
end
