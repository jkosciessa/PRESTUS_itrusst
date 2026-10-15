clear all; close all; clc;

% PRESTUS_THICKNESS_06_SOUNDSPEED_SIM
%
% Skull sound speed x skull thickness: c_skull at 9 evenly spaced values
% between its literature bounds plus the mid value (10 values,
% data/tissue_parameters.csv) x 3 skull thicknesses (3, 6.5, 9 mm; phantoms
% from PRESTUS_thickness_01_createPhantoms.m), all other tissue properties
% at their mid level: 30 independent evaluations with heating.
%
% Follow-up to the non-monotonic skull sound-speed effect (sub-002: brain
% I_SPPA at mid below both bounds, in PRESTUS and BabelBrain). A skull layer
% transmits least when its thickness is an odd multiple of lambda/4 and
% most at multiples of lambda/2, so the c_skull curve should change shape
% with thickness: monotonic at 3 mm, one dip near mid at 6.5 mm, an interior
% transmission peak at 9 mm (plane-wave layer estimate at 500 kHz).
%
% Thickness is baked into each phantom, so it is varied via subject_id
% (6.5 mm = sub-356, the 2D counterpart of sub-002). Shortened thermal
% protocol (20 s stimulation, no post-stimulation period), as in
% PRESTUS_thickness_04_density_sim.m, whose 'mid' runs this reproduces.
%
% Output: data/tussim/<variant_name>/sub-XXX/sub-XXX_phantom_s_c_<c>.csv,
% plus the design table read by PRESTUS_thickness_07_soundspeed_plot.m.

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code', 'functions'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

setup_name    = 'itrusst_protocol';
variant_name  = 'itrusst_protocol1a_soundspeed_thickness';
[thickness_subjects, thicknesses_mm] = itrusst_thickness_subjects([3 6.5 9]);
n_steps       = 9;   % evenly spaced c_skull values from bottom to top (+ mid)

%% ── Detect job-submission platform (auto mode) ──────────────────────────────

try
    hpc_type = hpc_detect_system();
catch
    hpc_type = 'matlab';
end
fprintf('Job submission platform: %s\n', hpc_type);

%% ── Load literature bounds + base config ────────────────────────────────────

params_table = readtable(fullfile(pn.data, 'tissue_parameters.csv'));

skull_c = @(level) params_table.sound_speed(strcmp(params_table.tissue, 'skull') & ...
                                            strcmp(params_table.level, level));
c_values = unique(round([linspace(skull_c('bottom'), skull_c('top'), n_steps), skull_c('mid')]));

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

base_parameters.hpc.timelimit = '00:30:00';  % as PRESTUS_thickness_04

tissue_labels = {'brain', 'skull'};
property_labels = {'sound_speed', 'density', 'alpha_coeff', 'thermal_conductivity', ...
                    'specific_heat_capacity', 'perfusion', 'absorption_fraction'};

%% ── Build the design table ──────────────────────────────────────────────────

n_rows = numel(thickness_subjects) * numel(c_values);
design_rows = table('Size', [n_rows 4], ...
    'VariableTypes', {'double', 'double', 'double', 'string'}, ...
    'VariableNames', {'subject_id', 'thickness_mm', 'skull_sound_speed', 'output_affix'});
row_i = 0;

for s = 1:numel(thickness_subjects)
    for c = 1:numel(c_values)
        row_i = row_i + 1;
        affix = sprintf('_s_c_%d', c_values(c));
        design_rows(row_i, :) = {thickness_subjects(s), thicknesses_mm(s), c_values(c), string(affix)};
    end
end

writetable(design_rows, fullfile(pn.path_sim, 'soundspeed_thickness_design.csv'));
fprintf('Sound speed x thickness design: %d evaluations (%d thicknesses x %d sound speeds: %s m/s).\n', ...
    n_rows, numel(thickness_subjects), numel(c_values), mat2str(c_values));

%% ── Submit every evaluation (all independent) ───────────────────────────────

for row_i = 1:n_rows
    subject_id = design_rows.subject_id(row_i);

    parameters = base_parameters;
    parameters.subject_id = subject_id;
    parameters.path.t1_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
    parameters.path.t2_pattern = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");

    % all properties at midpoint, then the swept factor
    for t = 1:numel(tissue_labels)
        tissue = tissue_labels{t};
        for p = 1:numel(property_labels)
            property = property_labels{p};
            parameters.medium_properties.(tissue).(property) = ...
                params_table.(property)(strcmp(params_table.tissue, tissue) & ...
                                         strcmp(params_table.level, 'mid'));
        end
    end
    parameters.medium_properties.skull.sound_speed = design_rows.skull_sound_speed(row_i);

    parameters.io.output_affix = char(design_rows.output_affix(row_i));
    prestus_pipeline_start(parameters);
end
