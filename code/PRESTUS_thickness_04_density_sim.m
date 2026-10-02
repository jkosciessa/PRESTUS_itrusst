clear all; close all; clc;

% PRESTUS_THICKNESS_04_DENSITY_SIM
%
% Skull density (bottom/mid/top) x skull thickness (1-12 mm, 0.5 mm steps;
% phantoms from PRESTUS_thickness_01_createPhantoms.m), all other tissue
% properties at their mid level: 3 x 23 = 69 independent evaluations with
% heating. Follow-up to the non-zero Morris sigma for skull density
% (explore/PRESTUS_paramvar_09/10). The 'mid' level reproduces the
% acoustic-only thickness_02 sweep, now with heating.
%
% Thickness is baked into each phantom, so it is varied via subject_id.
% Shortened thermal protocol (20 s stimulation, no post-stimulation period).
%
% Output: data/tussim/<variant_name>/sub-XXX/sub-XXX_phantom<affix>.csv,
% plus the design table read by PRESTUS_thickness_05_density_plot.m.

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

setup_name    = 'itrusst_protocol';
variant_name  = 'itrusst_protocol1a_density_thickness';
[thickness_subjects, thicknesses_mm] = itrusst_thickness_subjects(1:0.5:12);
density_levels_name = {'bottom', 'mid', 'top'};   % skull density levels to sweep

%% ── Detect job-submission platform (auto mode) ──────────────────────────────

try
    hpc_type = hpc_detect_system();
catch
    hpc_type = 'matlab';
end
fprintf('Job submission platform: %s\n', hpc_type);

%% ── Load literature bounds + base config ────────────────────────────────────

params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));

cd(fullfile(rootpath, 'data', 'configs'))
base_parameters = load_parameters(['config_', setup_name, '.yaml']);

if ismac
    base_parameters.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
else
    base_parameters.startup.simnibs_bin_path = pn.simnibs;
end
pn.path_sim = fullfile(rootpath, 'data', 'tussim', variant_name);
if ~exist(pn.path_sim, 'dir'); mkdir(pn.path_sim); end
base_parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";
base_parameters.path.anat = pn.path_seg;
base_parameters.path.seg  = pn.path_seg;
base_parameters.path.sim  = pn.path_sim;
base_parameters.startup.paths_to_add = {pn.kwave, pn.minimize};
base_parameters.pct.enabled = 0;

base_parameters.platform = hpc_type;
base_parameters.simulation.medium = 'phantom';
base_parameters.simulation.interactive = 0;
base_parameters.simulation.code_type = 'cpp_cpu';
base_parameters.simulation.precision = 'single';
base_parameters.io.overwrite_files = 'never';
base_parameters.io.overwrite_simnibs = 0;
base_parameters.io.save_heatingvideo = 0;
base_parameters.io.save_matrices = 0;
base_parameters.io.save_thermal_matrices = 0;
base_parameters.io.save_acoustic_matrices = 0;
base_parameters.modules.run_nifti_creation = 0;
base_parameters.modules.generate_report = 0;
base_parameters.modules.run_source_setup  = 1;
base_parameters.modules.run_acoustic_sims = 1;
base_parameters.modules.run_heating_sims  = 1;
base_parameters.grid.resolution_mm = 0.5;
base_parameters.grid.default_dims = [70, 120] / base_parameters.grid.resolution_mm;
base_parameters.grid.pml_size = [20] / base_parameters.grid.resolution_mm;
base_parameters.grid.axisymmetric = 1;

% shortened thermal protocol: 20 s stimulation, no post-stimulation period
base_parameters.timing.ptrd = 20;
base_parameters.timing.post_ptri_dur = 0;

base_parameters.hpc.timelimit = '00:30:00';  % 00:15:00 was too short for the thermal stage

tissue_labels = {'brain', 'skull'};
property_labels = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                    'specific_heat_capacity', 'perfusion', 'absorption_fraction'};

%% ── Build the design table + job specs ──────────────────────────────────────

n_rows = numel(thickness_subjects) * numel(density_levels_name);
design_rows = table('Size', [n_rows 4], ...
    'VariableTypes', {'double', 'double', 'string', 'string'}, ...
    'VariableNames', {'subject_id', 'thickness_mm', 'skull_density_level', 'output_affix'});
row_i = 0;

for s = 1:numel(thickness_subjects)
    for d = 1:numel(density_levels_name)
        row_i = row_i + 1;
        affix = sprintf('_s_rho_%s', density_levels_name{d});
        design_rows(row_i, :) = {thickness_subjects(s), thicknesses_mm(s), ...
            string(density_levels_name{d}), string(affix)};
    end
end

writetable(design_rows, fullfile(pn.path_sim, 'density_thickness_design.csv'));
fprintf('Density x thickness design: %d evaluations (%d thicknesses x %d density levels).\n', ...
    n_rows, numel(thickness_subjects), numel(density_levels_name));

%% ── Submit every evaluation (all independent) ───────────────────────────────

for row_i = 1:n_rows
    subject_id  = design_rows.subject_id(row_i);
    density_lvl = char(design_rows.skull_density_level(row_i));
    affix       = char(design_rows.output_affix(row_i));

    parameters = base_parameters;
    parameters.subject_id = subject_id;
    parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
    parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");

    for t = 1:numel(tissue_labels)
        tissue = tissue_labels{t};
        for p = 1:numel(property_labels)
            property = property_labels{p};
            if strcmp(tissue, 'skull') && strcmp(property, 'density')
                level_name = density_lvl;   % the swept factor
            else
                level_name = 'mid';         % everything else held at midpoint
            end
            parameters.medium_properties.(tissue).(property) = ...
                params_table.(property)(strcmp(params_table.tissue, tissue) & ...
                                         strcmp(params_table.level, level_name));
        end
    end

    parameters.io.output_affix = affix;
    prestus_pipeline_start(parameters);
end
