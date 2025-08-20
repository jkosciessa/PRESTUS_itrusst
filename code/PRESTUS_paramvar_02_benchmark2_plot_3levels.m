% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1', 'sub-002');
pn.code = fullfile(rootpath, 'code'); addpath(pn.code)

% Define tissues and suffixes
tissue_labels = {'brain', 'skull'};
tissue_suffixes = {'_b', '_s'};

% Define properties and suffixes
property_labels = {'soundSpeed', 'density', 'attenuation', ...
                   'thermConductivity', ...
                   'heatCapacity', 'perfusion', ...
                   'absorption'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};

% Define sweep values (2 levels)
sweep_values = 1:2;

% Preallocate results structure
results = struct();
xlabels = {};

% Additional metrics storage
max_Isppa_brain = nan(length(sweep_values), length(property_labels));
riseT_brain = nan(length(sweep_values), length(property_labels));

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
        maxIsppa_vec = nan(length(sweep_values),1);
        riseT_vec = nan(length(sweep_values),1);

        for v = 1:length(sweep_values)
            affix = [tissue_suffix prop_suffix num2str(sweep_values(v))];
            xlabels{i_cond} = [tissue_suffix(2:end) '_' property];
            filename = fullfile(pn.data, ['sub-002_phantom_output_table' affix '.csv']);

            if isfile(filename)
                T = readtable(filename);
                if all(ismember({'isppa_at_target'}, T.Properties.VariableNames))
                    isppa_vec(v) = T.isppa_at_target(1);
                    maxIsppa_vec(v) = T.max_Isppa_brain(1);
                else
                    warning('File %s does not contain required columns.', filename);
                end
                if all(ismember({'maxT'}, T.Properties.VariableNames))
                    maxT_vec(v) = T.maxT(1);
                    riseT_vec(v) = T.riseT_brain(1);
                else
                    warning('File %s does not contain required columns.', filename);
                end
            else
                warning('File %s not found.', filename);
            end
        end
        results.sweep(:,i_cond) = sweep_values;
        results.isppa_at_target(:,i_cond) = isppa_vec;
        results.maxT(:,i_cond) = maxT_vec - 37;
        results.max_Isppa_brain(:, i_cond) = maxIsppa_vec;
        results.riseT_brain(:, i_cond) = riseT_vec;
    end
end

%% Restrict to only plot levels 1 and 2
scatter_levels = [1, 2]; % Only first and second sweep levels
scatter_colors = {[1 1 1], [0 0 0]}; % white and black

%% Plot maximum temperature
plot_metric_with_patches(results.maxT, xlabels, 'max. Temp rise (deg C)', 'sub-002_temperature', [0 4], scatter_levels, scatter_colors, pn);

%% Plot target intensity
plot_metric_with_patches(results.isppa_at_target, xlabels, 'target intensity', 'sub-002_intensity', [], scatter_levels, scatter_colors, pn);

%% Plot max_Isppa_brain
plot_metric_with_patches(results.max_Isppa_brain, xlabels, 'max Isppa brain', 'sub-002_maxIsppa_brain', [], scatter_levels, scatter_colors, pn);

%% Plot riseT_brain
plot_metric_with_patches(results.riseT_brain, xlabels, 'rise T brain (deg C)', 'sub-002_riseT_brain', [], scatter_levels, scatter_colors, pn);
