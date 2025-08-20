function plot_metric_with_patches(data_matrix, xlabels, y_label, fig_name, y_limits, ...
    scatter_levels, scatter_colors, pn)

    parameters = 1:size(data_matrix,2);
    result_mean = nanmean(data_matrix, 1);
    result_std = nanstd(data_matrix, 0, 1); % sample std
    result_max = nanmax(data_matrix, [], 1);
    result_min = nanmin(data_matrix, [], 1);

    barWidth = 0.6;

    h = figure;
    set(h, 'Units', 'normalized');
    if numel(xlabels) <10
        set(h, 'Position', [0.25 0.25 0.15 0.3]);
    else
        set(h, 'Position', [0.25 0.25 0.25 0.3]);
    end
    hold on;
    
    lower_bounds = result_min;
    upper_bounds = result_max;

    lightGrey = [0.5 0.5 0.5];  % Light grey color for background

    for i = 1:numel(parameters)
        x = [parameters(i)-barWidth/2, parameters(i)+barWidth/2, parameters(i)+barWidth/2, parameters(i)-barWidth/2];
        y = [lower_bounds(i), lower_bounds(i), upper_bounds(i), upper_bounds(i)];
        % If xlabel starts with 'b', draw a light grey rectangle behind the bar
        if startsWith(xlabels{i}, 'b', 'IgnoreCase', true)
            patch(x, y, lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
        elseif startsWith(xlabels{i}, 's', 'IgnoreCase', true)
            % Fill bar with red
            patch(x, y, [0.7 0 0], 'EdgeColor', 'none');
        end
    end

    % % Redraw bars on top of background patches to keep bar visible
    % for i = 1:numel(parameters)
    %     x = [parameters(i)-barWidth/2, parameters(i)+barWidth/2, parameters(i)+barWidth/2, parameters(i)-barWidth/2];
    %     y = [lower_bounds(i), lower_bounds(i), upper_bounds(i), upper_bounds(i)];
    %     if startsWith(xlabels{i}, 'b', 'IgnoreCase', true)
    %         patch(x, y, lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
    %     elseif startsWith(xlabels{i}, 's', 'IgnoreCase', true)
    %         % Fill bar with red
    %         patch(x, y, [0.7 0 0], 'EdgeColor', 'none');
    %     end
    % end

    for p = 1:numel(parameters)
        scatter(parameters(p), data_matrix(scatter_levels(1),p), ...
            60, scatter_colors{1}, 'filled', ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1);
        scatter(parameters(p), data_matrix(scatter_levels(2),p), ...
            60, scatter_colors{2}, 'filled', ...
            'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1);
    end

    set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
    xtickangle(45);
    ylabel(y_label);
    if ~isempty(y_limits)
        ylim(y_limits);
    end
    set(gca, 'Layer', 'top');
    h.Color = 'white';
    h.InvertHardcopy = 'off';
    saveas(h, fullfile(pn.figures, fig_name), 'epsc');
    saveas(h, fullfile(pn.figures, fig_name), 'png');
end
