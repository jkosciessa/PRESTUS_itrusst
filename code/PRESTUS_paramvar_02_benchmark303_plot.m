% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1', 'sub-303');
pn.code    = fullfile(rootpath, 'code'); addpath(pn.code)

subject_id        = 303;
tissue_labels     = {'brain', 'skull'};
tissue_suffixes   = {'_b', '_s'};
property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
sweep_values      = 1:2;

[results, xlabels] = load_paramvar_results(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values);

%% Restrict to only plot levels 1 and 2
scatter_levels = [1, 2]; % Only first and second sweep levels
scatter_colors = {[1 1 1], [0 0 0]}; % white and black

plot_metric_with_patches(results.maxT,                   xlabels, 'max. Temp rise (deg C)',  'sub-303_temperature',        [0 4],  scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.Ipa_target,             xlabels, 'target intensity',          'sub-303_intensity',          [0 20], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.Isppa_brain,            xlabels, 'max Isppa brain',           'sub-303_maxIsppa_brain',     [0 20], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.riseT_brain,            xlabels, 'rise T brain (deg C)',      'sub-303_riseT_brain',        [0 2],  scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.real_focal_distance_mm, xlabels, 'real focal distance (mm)',  'sub-303_real_focal_distance', [],    scatter_levels, scatter_colors, pn);

%% Plot differences

% Plot maximum temperature as difference between levels
parameters = 1:size(results.maxT,2);
result_diff = results.maxT(2,:) - results.maxT(1,:); % difference between high and low levels

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.25 0.3]); % 50% width and height
hold on;

bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);

set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
xtickangle(45);
ylabel('Difference in max. Temperature (deg C)');
set(gca, 'Layer', 'top');
h.Color = 'white';
h.InvertHardcopy = 'off';

% Plot target intensity as difference between levels
parameters = 1:size(results.Ipa_target,2);
result_diff = results.Ipa_target(2,:) - results.Ipa_target(1,:); % difference between high and low levels

h = figure;
set(h, 'Units', 'normalized');
set(h, 'Position', [0.25 0.25 0.25 0.3]); % 50% width and height
hold on;

bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);

set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', strrep(xlabels, '_', ' '));
xtickangle(45);
yline(0, '-k', 'LineWidth', 1); % line at y=0
ylabel('Difference in target intensity');
set(gca, 'Layer', 'top');
h.Color = 'white';
h.InvertHardcopy = 'off';
