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
                   'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', ...
                   'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

% Define sweep values (2 levels)
sweep_values = 1:2;

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
            % Construct the filename affix
            affix = [tissue_suffix prop_suffix num2str(sweep_values(v))];
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

%% Restrict to only plot levels 1 and 2
scatter_levels = [1, 2]; % Only first and second sweep levels
scatter_colors = {[1 1 1], [0 0 0]}; % white and black

%% Plot maximum temperature as difference between levels
parameters = 1:size(results.maxT,2);
result_diff = results.maxT(2,:) - results.maxT(1,:); % difference between high and low levels

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.5 0.3]); % 50% width and height
hold on;

bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);

set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
xtickangle(45);
ylabel('Difference in max. Temperature (deg C)');
set(gca, 'Layer', 'top');
h.Color = 'white';
h.InvertHardcopy = 'off';

%% Plot target intensity as difference between levels
parameters = 1:size(results.isppa_at_target,2);
result_diff = results.isppa_at_target(2,:) - results.isppa_at_target(1,:); % difference between high and low levels

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.5 0.3]); % 50% width and height
hold on;

bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);

set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
xtickangle(45);
yline(0, '-k', 'LineWidth', 1); % line at y=0
ylabel('Difference in target intensity');
set(gca, 'Layer', 'top');
h.Color = 'white';
h.InvertHardcopy = 'off';
