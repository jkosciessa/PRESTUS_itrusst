function make_parallel_plot(data, property_labels, setup_labels, setup_colors, ...
    flip_cols, pn, fig_name)
% MAKE_PARALLEL_PLOT
%
% Build one parallel-coordinates figure comparing setups (rows of `data`)
% across properties (columns of `data`). Each property gets its own
% vertical axis, scaled independently to its own min/max across setups
% (classic parallel-coordinates convention), so that a fast-varying
% property is not swamped by a slow-varying one on a shared axis.
% Connecting lines trace each setup across all properties.
%
% Columns listed in flip_cols have their axis direction reversed (max at
% bottom, min at top); all others run min (bottom) to max (top).
%
% Missing values (NaN) are skipped: the line breaks across that property
% and is bridged with a dashed segment instead of being interpolated.
%
% Used by PRESTUS_itrusst_setups_plot.m.

    n_setups = size(data, 1);
    n_props  = size(data, 2);

    % per-property min-max normalisation, with a small margin so extreme
    % points aren't drawn on the axis frame
    margin = 0.08;
    prop_min = min(data, [], 1, 'omitnan');
    prop_max = max(data, [], 1, 'omitnan');
    prop_range = prop_max - prop_min;
    prop_range(prop_range == 0) = 1; % guard against constant columns

    norm_data = (data - prop_min) ./ prop_range;
    norm_data(:, flip_cols) = 1 - norm_data(:, flip_cols);
    norm_data = norm_data * (1 - 2*margin) + margin;

    h = figure;
    set(gcf, 'renderer', 'Painters');
    set(h, 'Units', 'normalized');
    set(h, 'Position', [0.2 0.2 0.35 0.45]);
    hold on;

    x = 1:n_props;

    % vertical axis lines, one per property
    for p = 1:n_props
        plot([x(p) x(p)], [0 1], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % connecting lines per setup (NaNs break the line automatically).
    % Setups with (near-)identical values would otherwise draw exactly on
    % top of each other, so each setup's line is nudged by a small,
    % consistent horizontal jitter, and drawn with slight transparency so
    % any residual overlap still shows the lines underneath.
    jitter_span = 0.45;
    if n_setups > 1
        jitter = linspace(-jitter_span/2, jitter_span/2, n_setups);
    else
        jitter = 0;
    end

    % Setups are colour-coded, but when two or more share the same colour
    % (e.g. because a shared colour is used to group them visually) they
    % are told apart by marker shape instead: the first occurrence of a
    % colour gets a circle, the next occurrence of that same colour gets
    % a square, and so on through the marker list below.
    marker_cycle = {'o', 's', '^', 'd', 'v', 'p'};
    setup_markers = cell(1, n_setups);
    for s = 1:n_setups
        occurrence = sum(cellfun(@(c) isequal(c, setup_colors{s}), setup_colors(1:s)));
        setup_markers{s} = marker_cycle{mod(occurrence - 1, numel(marker_cycle)) + 1};
    end

    lines = gobjects(n_setups, 1);
    for s = 1:n_setups
        y = norm_data(s, :);
        lines(s) = plot(x + jitter(s), y, ['-', setup_markers{s}], ...
            'Color', [setup_colors{s}, 0.85], 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'MarkerFaceColor', setup_colors{s}, 'MarkerEdgeColor', 'none');

        % A missing value leaves a gap that the solid line above simply
        % skips. Bridge each such gap with a dashed segment connecting the
        % last valid point before it to the first valid point after it, to
        % show the missing value explicitly rather than silently dropping
        % that stretch of line.
        is_valid = ~isnan(y);
        gap_starts = find(is_valid(1:end-1) & ~is_valid(2:end));
        for gs = gap_starts
            ge = gs + find(is_valid(gs+1:end), 1, 'first');
            if ~isempty(ge)
                plot(x([gs ge]) + jitter(s), y([gs ge]), '--', ...
                    'Color', [setup_colors{s}, 0.85], 'LineWidth', 1.8, ...
                    'HandleVisibility', 'off');
            end
        end
    end

    % per-property min/max value labels at the bottom/top of each axis
    % (flipped columns show max at the bottom and min at the top).
    % The minimum value is always coloured blue, the maximum always red,
    % regardless of which end of the axis it ends up on.
    min_color = [0.00 0.00 0.85];
    max_color = [0.85 0.00 0.00];
    for p = 1:n_props
        if ismember(p, flip_cols)
            bottom_val = prop_max(p); bottom_color = max_color;
            top_val    = prop_min(p); top_color    = min_color;
        else
            bottom_val = prop_min(p); bottom_color = min_color;
            top_val    = prop_max(p); top_color    = max_color;
        end
        text(x(p), -0.06, fmtnum(bottom_val), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', 'FontSize', 9, 'Color', bottom_color, 'FontWeight', 'bold');
        text(x(p), 1.06, fmtnum(top_val), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 9, 'Color', top_color, 'FontWeight', 'bold');
    end

    set(gca, 'XTick', x, 'XTickLabel', property_labels, 'TickLabelInterpreter', 'tex', 'FontSize', 12);
    xtickangle(30);
    set(gca, 'YTick', [], 'YColor', 'none');
    xlim([0.5, n_props + 0.5]);
    ylim([-0.12, 1.12]);
    set(gca, 'Layer', 'top', 'Box', 'off');

    legend(lines, setup_labels, 'Location', 'eastoutside', 'Box', 'off');

    h.Color = 'white';
    h.InvertHardcopy = 'off';

    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end

function s = fmtnum(v)
% Plain fixed-point formatting (never scientific notation): integers
% print without a decimal point, non-integers keep up to 3 decimals
% with trailing zeros trimmed.
    if v == round(v)
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.3f', v);
        s = regexprep(s, '0+$', '');
        s = regexprep(s, '\.$', '');
    end
end
