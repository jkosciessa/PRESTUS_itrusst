% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures', 'paramvar'); mkdir(pn.figures);
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-002');
pn.code    = fullfile(rootpath, 'code'); addpath(fullfile(pn.code, 'functions'))

subject_id        = 002;
tissue_labels     = {'brain', 'skull'};
tissue_suffixes   = {'_b', '_s'};
property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
sweep_values      = 1:2;

[results, xlabels] = load_paramvar_results(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols);

%% Restrict to only plot levels 1 and 2
scatter_levels = [1, 2]; % Only first and second sweep levels
scatter_colors = {[1 1 1], [0 0 0]}; % white and black

% plot_metric_with_patches(results.maxT,                   xlabels, 'max. Temp rise (deg C)',  'sub-002_temperature',        [0 2.5], scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.Ipa_target,             xlabels, 'target intensity',          'sub-002_intensity',          [0 15],  scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.Isppa_brain,            xlabels, 'max Isppa brain',           'sub-002_maxIsppa_brain',     [0 17],  scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.riseT_brain,            xlabels, 'rise T brain (deg C)',      'sub-002_riseT_brain',        [0 1.4], scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.Isppa_skull,            xlabels, 'max Isppa skull',           'sub-002_maxIsppa_skull',     [0 11],  scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.riseT_skull,            xlabels, 'rise T skull (deg C)',      'sub-002_riseT_skull',        [0 2.5], scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results.real_focal_distance_mm, xlabels, 'real focal distance (mm)',  'sub-002_real_focal_distance', [],     scatter_levels, scatter_colors, pn);

%% Alternative: temperature — brain tissue only; intensity — acoustic properties only

[results_temp, xlabels_temp] = load_paramvar_results(pn.data, subject_id, ...
    {'skull'}, {'_s'}, ...
    property_labels, property_suffixes, sweep_values, property_symbols);

% plot_metric_with_patches(results_temp.maxT,        xlabels_temp, 'max. Temp rise (deg C)', 'sub-002_temperature_skullonly',    [0 2.5], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results_temp.riseT_brain,  xlabels_temp, 'rise T brain (deg C)',   'sub-002_riseT_brain_skullonly',    [0 0.8], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results_temp.riseT_skull,  xlabels_temp, 'rise T skull (deg C)',   'sub-002_riseT_skull_skullonly',    [0 2.6], scatter_levels, scatter_colors, pn);

[results_int, xlabels_int] = load_paramvar_results(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, ...
    {'sound_speed', 'density', 'alpha_coeff'}, {'_c_', '_rho_', '_a0_'}, sweep_values, {'c', '\rho', '\alpha'});

plot_metric_with_patches(results_int.Ipa_target,             xlabels_int, 'target intensity',        'sub-002_intensity_acousticonly',          [0 15], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results_int.Isppa_brain,            xlabels_int, 'max Isppa brain',          'sub-002_maxIsppa_brain_acousticonly',     [0 17], scatter_levels, scatter_colors, pn, 0:2:16);
% plot_metric_with_patches(results_int.Isppa_skull,            xlabels_int, 'max Isppa skull',          'sub-002_maxIsppa_skull_acousticonly',     [0 11], scatter_levels, scatter_colors, pn);
% plot_metric_with_patches(results_int.real_focal_distance_mm, xlabels_int, 'real focal distance (mm)', 'sub-002_real_focal_distance_acousticonly', [],    scatter_levels, scatter_colors, pn);

%% Focal geometry and heating-location comparison — all tissues x all properties

[focus_results, ~] = load_focus_heating_metrics(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols);

plot_metric_with_patches(results.real_focal_distance_mm, xlabels, 'real focal distance (mm)',       'sub-002_real_focal_distance_mm', [],  scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.focal_vol_3dB_mm3,      xlabels, 'focal volume, -3dB/FWHM (mm^3)', 'sub-002_focal_vol_3dB_mm3',      [],  scatter_levels, scatter_colors, pn);

% Heating-location comparison: whole-brain max (results.riseT_brain, plotted
% above for the skull-only subset) vs at the empirical acoustic focus vs at
% the fixed nominal target — see load_focus_heating_metrics.m
plot_metric_with_patches(focus_results.riseT_focus,  xlabels, 'rise T at empirical focus (deg C)',      'sub-002_riseT_focus',  [0 0.8], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(focus_results.riseT_target, xlabels, 'rise T at fixed nominal target (deg C)', 'sub-002_riseT_target', [0 0.8], scatter_levels, scatter_colors, pn);

%% Plot differences

% % Plot maximum temperature as difference between levels
% parameters = 1:size(results.maxT,2);
% result_diff = results.maxT(2,:) - results.maxT(1,:); % difference between high and low levels
% 
% h = figure;
% set(h, 'Units', 'normalized');
% set(h, 'Position', [0.25 0.25 0.25 0.3]); % 50% width and height
% hold on;
% 
% bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);
% 
% set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', xlabels, 'TickLabelInterpreter', 'tex', 'FontSize', 15);
% xtickangle(45);
% ylabel('Difference in max. Temperature (deg C)');
% set(gca, 'Layer', 'top');
% h.Color = 'white';
% h.InvertHardcopy = 'off';
% 
% % Plot target intensity as difference between levels
% parameters = 1:size(results.Ipa_target,2);
% result_diff = results.Ipa_target(2,:) - results.Ipa_target(1,:); % difference between high and low levels
% 
% h = figure;
% set(h, 'Units', 'normalized');
% set(h, 'Position', [0.25 0.25 0.25 0.3]); % 50% width and height
% hold on;
% 
% bar(parameters, result_diff, 'FaceColor', [.7 .0 .0], 'EdgeColor',[1 1 1]);
% 
% set(gca, 'XTick', 1:numel(parameters), 'XTickLabel', xlabels, 'TickLabelInterpreter', 'tex', 'FontSize', 15);
% xtickangle(45);
% yline(0, '-k', 'LineWidth', 1); % line at y=0
% ylabel('Difference in target intensity');
% set(gca, 'Layer', 'top');
% h.Color = 'white';
% h.InvertHardcopy = 'off';

