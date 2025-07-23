
% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1', 'sub-002');

% Define tissues and suffixes
tissue_labels = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};

% Define properties and suffixes
property_labels = {'sound_speed', 'density', 'alpha_0_true', ...
                   'alpha_power_true', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', ...
                   'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_y_', '_k_', '_heatc_', '_perf_', '_abs_'};

% Define sweep values (10% to 200%)
sweep_values = 0.1:0.2:2.0;

% Preallocate results structure
results = struct();
xlabels = {};

i_cond = 0;
for t = 1:length(tissue_labels)
    tissue = tissue_labels{t};
    tissue_suffix = tissue_suffixes{t};

    for p = 1:length(property_labels)
        i_cond = i_cond+1;
        property = property_labels{p};
        prop_suffix = property_suffixes{p};

        % Prepare to store results for this tissue/property
        isppa_vec = nan(length(sweep_values),1);
        maxT_vec = nan(length(sweep_values),1);

        for v = 1:length(sweep_values)
            sweep_percent = 100 * sweep_values(v);
            % Construct the filename affix
            affix = [tissue_suffix prop_suffix num2str(sweep_percent)];
            xlabels{i_cond} = [tissue_suffix(2:end) '_' property];
            % Construct the filename
            filename = fullfile(pn.data, ['sub-002_phantom_output_table' affix '.csv']);

            if isfile(filename)
                T = readtable(filename);
                if all(ismember({'isppa_at_target','maxT'}, T.Properties.VariableNames))
                    isppa_vec(v) = T.isppa_at_target(1);
                    maxT_vec(v) = T.maxT(1);
                else
                    warning('File %s does not contain required columns.', filename);
                end
            else
                warning('File %s not found.', filename);
            end
        end
        % Store results in a structured way
        results.sweep(:,i_cond) = sweep_values;
        results.isppa_at_target(:,i_cond) = isppa_vec;
        results.maxT(:,i_cond) = maxT_vec;
    end
end

%% restrict the ranges based on empirical ranges
% This affects the individual data points.

rangeToPlot{1} = 10;
rangeToPlot{2} = 10;
rangeToPlot{3} = 1:20;
rangeToPlot{4} = 2:20;
rangeToPlot{5} = 3:11;
rangeToPlot{6} = 10;
rangeToPlot{7} = 4:14;
rangeToPlot{8} = 10:12;

rangeToPlot{8+1} = 7:11;
rangeToPlot{8+2} = 6:10;
rangeToPlot{8+3} = 1:20;
rangeToPlot{8+4} = 3:20;
rangeToPlot{8+5} = 10:16;
rangeToPlot{8+6} = 10:17;
rangeToPlot{8+7} = 5:15;
rangeToPlot{8+8} = 5:20;

a = 1; b = 20; c = 1; d = 10;
for x = 1:16
    rangeToPlot{x} = unique(round(c + (d-c) * (rangeToPlot{x} - a) / (b - a)));
    rangeLowerUpper(x,:) = [min(rangeToPlot{x}), max(rangeToPlot{x})];
end

%% Plot maximum temperature

% Example plot
% figure;
% imagesc(results.maxT)

parameters = 1:size(results.maxT,2);
sweep_values = results.sweep(:,1);

result_lower = []; result_upper = [];
result_mean = median(results.maxT, 1, 'omitnan');
% result_mean = results.maxT(5,:);
for i_param = 1:size(results.maxT,2)
    result_lower(i_param) = results.maxT(rangeLowerUpper(i_param,1),i_param);%-result_mean(i_param);%std(results.maxT, 0, 1, 'omitnan'); % 0 for sample standard deviation
    result_upper(i_param) = results.maxT(rangeLowerUpper(i_param,2),i_param);%-result_mean(i_param);%std(results.maxT, 0, 1, 'omitnan'); % 0 for sample standard deviation
end

h = figure;
set(gcf, 'Renderer', 'painters');
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.35 0.3]); % 50% width and height, starting at 25% from left and bottom
    hold on;
    
    % % 1. Plot error bars (mean ± std)
    % bar(parameters, result_mean, 'FaceColor', [.9 .9 .9], 'EdgeColor',[1 1 1])
    % errorbar(parameters, result_mean, result_lower, result_upper, 'Color', [.7 .5 .5], 'LineStyle', 'none', ...
    %     'CapSize', 3, 'LineWidth', 1);

    for i_param = 1:16
        try
        rectangle('Position', ...
            [i_param-0.4, result_lower(i_param), 0.8, result_upper(i_param)-result_lower(i_param)],  ...
            'FaceColor', [.9 .9 .9]);
        catch
        end
    end
    
    % 2. Scatter plot of individual sweeps with color coding
    cmap = gray(length(sweep_values)); % Color map for sweeps
    for p = 1:numel(parameters)
        scatter(...
            repmat(parameters(p), [1 numel(rangeToPlot{p})]),... % X-values
            results.maxT(rangeToPlot{p},p),...                                % Y-values
            40,...                                               % Marker size
            cmap(rangeToPlot{p},:), ...
            'filled',...
            'MarkerEdgeColor', [.5 .5 .5],...
            'MarkerFaceAlpha', 0.6...
        );
    end
    
    % Add colorbar to show sweep values
    c = colorbar;
    caxis([10 200]);
    c.Label.String = 'Parameter scaling';
    colormap(gray);
    
    % Set x-axis labels and rotate
    set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
    xtickangle(45);
    
    ylabel('max. Temperature (deg C)');
    ylim([37 40])
    
    set(gca, 'Layer', 'top');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);

figureName = ['sub-002_temperature'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');

%% Plot target intensity

parameters = 1:size(results.isppa_at_target,2);
sweep_values = results.sweep(:,1);

result_lower = []; result_upper = [];
result_mean = median(results.isppa_at_target, 1, 'omitnan');
% result_mean = results.maxT(5,:);
for i_param = 1:size(results.isppa_at_target,2)
    result_lower(i_param) = results.isppa_at_target(rangeLowerUpper(i_param,1),i_param);%-result_mean(i_param);%std(results.maxT, 0, 1, 'omitnan'); % 0 for sample standard deviation
    result_upper(i_param) = results.isppa_at_target(rangeLowerUpper(i_param,2),i_param);%-result_mean(i_param);%std(results.maxT, 0, 1, 'omitnan'); % 0 for sample standard deviation
end

h = figure;
set(gcf, 'Renderer', 'painters');
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.35 0.3]); % 50% width and height, starting at 25% from left and bottom
    hold on;
    
    % % 1. Plot error bars (mean ± std)
    % bar(parameters, result_mean, 'FaceColor', [.9 .9 .8], 'EdgeColor',[1 1 1])
    % errorbar(parameters, result_mean, result_std, 'Color', [.7 .5 .5], 'LineStyle', 'none', ...
    %     'CapSize', 3, 'LineWidth', 1);
    
    for i_param = 1:16
        try
        rectangle('Position', ...
            [i_param-0.4, result_lower(i_param), 0.8, result_upper(i_param)-result_lower(i_param)],  ...
            'FaceColor', [.9 .9 .9]);
        catch
        end
    end

    % 2. Scatter plot of individual sweeps with color coding
    cmap = gray(length(sweep_values)); % Color map for sweeps
    for p = 1:numel(parameters)
        scatter(...
            repmat(parameters(p), [1 numel(rangeToPlot{p})]),... % X-values
            results.isppa_at_target(rangeToPlot{p},p),...                                % Y-values
            40,...                                               % Marker size
            cmap(rangeToPlot{p},:), ...
            'filled',...
            'MarkerEdgeColor', [.5 .5 .5],...
            'MarkerFaceAlpha', 0.6...
        );
    end
    
    % Add colorbar to show sweep values
    c = colorbar;
    caxis([10 200]);
    c.Label.String = 'Parameter scaling';
    colormap(gray);
    
    % Set x-axis labels and rotate
    set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
    xtickangle(45);
    
    ylabel('target intensity');
    yl = yline(0, '-k', 'LineWidth', 1); % dashed black line at y=0
    
    set(gca, 'Layer', 'top');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);

figureName = ['sub-002_intensity'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');