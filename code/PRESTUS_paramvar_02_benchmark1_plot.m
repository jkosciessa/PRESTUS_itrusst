% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.figures = fullfile(rootpath, 'figures');
pn.data    = fullfile(rootpath, 'data', 'tussim', 'itrusst_protocol1a', 'sub-001');
pn.code    = fullfile(rootpath, 'code'); addpath(fullfile(pn.code, 'functions'))

subject_id       = 1;
tissue_labels    = {'brain'};
tissue_suffixes  = {'_b'};
property_labels   = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                     'specific_heat_capacity', 'perfusion', 'absorption_fraction'};
property_suffixes = {'_c_', '_rho_', '_a0_', '_k_', '_heatc_', '_perf_', '_abs_'};
property_symbols  = {'c', '\rho', '\alpha', '\kappa', 'c_p', '\omega', '\beta'};
sweep_values      = 1:2;

[results, xlabels] = load_paramvar_results(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols);
[focus_results, ~] = load_focus_heating_metrics(pn.data, subject_id, ...
    tissue_labels, tissue_suffixes, property_labels, property_suffixes, sweep_values, property_symbols);

%% Restrict to only plot levels 1 and 2
scatter_levels = [1, 2]; % Only first and second sweep levels
scatter_colors = {[1 1 1], [0 0 0]}; % white and black

plot_metric_with_patches(results.maxT,                   xlabels, 'max. Temp rise (deg C)',   'sub-001_temperature',           [0 4],  scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.Ipa_target,             xlabels, 'target intensity',           'sub-001_intensity',             [0 60], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.Isppa_brain,            xlabels, 'max Isppa brain',            'sub-001_maxIsppa_brain',        [0 60], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.riseT_brain,            xlabels, 'rise T brain, whole-brain max (deg C)', 'sub-001_riseT_brain', [0 1.4], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.real_focal_distance_mm, xlabels, 'real focal distance (mm)',   'sub-001_real_focal_distance_mm', [],    scatter_levels, scatter_colors, pn);
plot_metric_with_patches(results.focal_vol_3dB_mm3,      xlabels, 'focal volume, -3dB/FWHM (mm^3)', 'sub-001_focal_vol_3dB_mm3', [],    scatter_levels, scatter_colors, pn);

% Heating-location comparison: whole-brain max (above) vs at the empirical
% acoustic focus vs at the fixed nominal target — see load_focus_heating_metrics.m
plot_metric_with_patches(focus_results.riseT_focus,  xlabels, 'rise T at empirical focus (deg C)', 'sub-001_riseT_focus',  [0 1.4], scatter_levels, scatter_colors, pn);
plot_metric_with_patches(focus_results.riseT_target, xlabels, 'rise T at fixed nominal target (deg C)', 'sub-001_riseT_target', [0 1.4], scatter_levels, scatter_colors, pn);
