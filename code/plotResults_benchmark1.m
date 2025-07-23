
% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1', 'sub-001');

% Define tissues and suffixes
tissue_labels = {'brain'};
tissue_suffixes = {'_b'};

% Define properties and suffixes
property_labels = {'sound_speed', 'density', 'alpha_0_true', ...
                   'alpha_power_true', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', ...
                   'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_y_', '_k_', '_heatc_', '_perf_', '_abs_'};

% Define sweep values (10% to 200%)
sweep_values = 0.1:0.1:2.0;

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
            filename = fullfile(pn.data, ['sub-001_phantom_output_table' affix '.csv']);

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

%% Plot maximum temperature

% Example plot
% figure;
% imagesc(results.maxT)

parameters = 1:size(results.maxT,2);
sweep_values = results.sweep(:,1);

result_mean = nanmean(results.maxT, 1);
result_std = nanstd(results.maxT, 0, 1); % 0 for sample standard deviation

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.5 0.3]); % 50% width and height, starting at 25% from left and bottom
    hold on;
    
    % 1. Plot error bars (mean ± std)
    bar(parameters, result_mean, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1])
    errorbar(parameters, result_mean, result_std, 'Color', [.7 .5 .5], 'LineStyle', 'none', ...
        'CapSize', 3, 'LineWidth', 1);
    
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
    ylim([37 43])
    
    set(gca, 'Layer', 'top');

figureName = ['sub-001_temperature'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');

%% Plot target intensity

parameters = 1:size(results.isppa_at_target,2);
sweep_values = results.sweep(:,1);

result_mean = nanmean(results.isppa_at_target, 1);
result_std = nanstd(results.isppa_at_target, 0, 1); % 0 for sample standard deviation

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.5 0.3]); % 50% width and height, starting at 25% from left and bottom
    hold on;
    
    % 1. Plot error bars (mean ± std)
    bar(parameters, result_mean, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1])
    errorbar(parameters, result_mean, result_std, 'Color', [.7 .5 .5], 'LineStyle', 'none', ...
        'CapSize', 3, 'LineWidth', 1);
    
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
    
    yl = yline(0, '-k', 'LineWidth', 1); % dashed black line at y=0

    ylabel('target intensity');
    
    set(gca, 'Layer', 'top');

figureName = ['sub-001_intensity'];
h.Color = 'white';
h.InvertHardcopy = 'off';
saveas(h, fullfile(pn.figures, figureName), 'epsc');
saveas(h, fullfile(pn.figures, figureName), 'png');