clear all; close all; clc;

user = 'julkos';

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% settings

subject_id        = 002;
setup_name        = 'itrusst_protocol1';
variant_name      = 'itrusst_protocol1a';

%% load parameter file and adjust paths

cd(fullfile(rootpath, 'data', 'configs'));
parameters = load_parameters(['config_', setup_name, '.yaml']);

if ismac
    parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
else
    parameters.startup.simnibs_bin_path = pn.simnibs;
end

pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
if ~exist(pn.path_sim, 'dir'); mkdir(pn.path_sim); end

parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";
parameters.path.anat            = pn.path_seg;
parameters.path.seg             = pn.path_seg;
parameters.path.t1_pattern      = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T1.nii.gz');
parameters.path.t2_pattern      = fullfile(sprintf('m2m_sub-%03d', subject_id), 'T2_reg.nii.gz');
parameters.path.sim             = pn.path_sim;
parameters.startup.paths_to_add = {pn.kwave, pn.minimize};

pn.outputs_folder = fullfile(pn.path_sim, sprintf('sub-%03d', subject_id));
if ~exist(pn.outputs_folder, 'dir'); mkdir(pn.outputs_folder); end

%% apply mid-level (default) tissue properties
% tissue_parameters.csv level='mid' — the baseline of the paramvar sweep;
% the config hardcodes slightly different values.

tissue_labels   = {'brain', 'skull'};
property_labels = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                   'specific_heat_capacity', 'perfusion', 'absorption_fraction'};

params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));
for t = 1:length(tissue_labels)
    for p = 1:length(property_labels)
        parameters.medium_properties.(tissue_labels{t}).(property_labels{p}) = ...
            params_table.(property_labels{p})(strcmp(params_table.tissue, tissue_labels{t}) & ...
                                               strcmp(params_table.level, 'mid'));
    end
end

%% simulation settings

parameters.simulation.medium      = 'phantom';
parameters.simulation.interactive = 0;
parameters.simulation.code_type   = 'cpp_gpu';
parameters.simulation.precision   = 'single';

parameters.io.overwrite_files    = 'never';
parameters.io.overwrite_simnibs  = 0;
parameters.io.save_heatingvideo  = 0;
parameters.io.save_matrices      = 0;

parameters.grid.resolution_mm  = 0.5;
parameters.grid.default_dims   = [70, 120] / parameters.grid.resolution_mm;
parameters.grid.pml_size       = [20]  / parameters.grid.resolution_mm; % 20 mm PML, matches primary pipeline configs
parameters.grid.axisymmetric   = 1;

parameters.pct.enabled = 0;

parameters.platform = 'matlab';

% all-mid baseline (counterpart to the bottom/top sweep of paramvar_01)
parameters.io.output_affix = '_baseline';

parameters.subject_id = subject_id;

%% run

prestus_pipeline_start(parameters);
