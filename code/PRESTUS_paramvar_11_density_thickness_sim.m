clear all; close all; clc;

% PRESTUS_PARAMVAR_11_DENSITY_THICKNESS_SIM
%
% Targeted follow-up on a signal from PRESTUS_paramvar_09_morris_
% interactions_sim.m's pilot run: skull density was the only one of the
% 14 tissue-property factors with a non-zero sigma (interaction/context-
% dependence) in that pilot, based on just n=2 trajectories - too thin to
% trust, but worth a cheap, targeted, controlled check rather than waiting
% for a much larger Morris r to pin it down indirectly.
%
% This script asks one specific version of "does skull density interact
% with something else": does its effect on Isppa_brain/Ipa_target/riseT_
% brain/riseT_skull change across the range of skull thickness? All other
% tissue properties (brain and skull) are held at their literature
% midpoint throughout, isolating skull density x skull thickness as a
% clean 2-way design (3 density levels x 12 thicknesses = 36 evaluations),
% not a statistical screen.
%
% GEOMETRY: reuses the existing thickness-sweep phantoms (subjects 301-
% 312, skull thickness 1-12 mm, built by PRESTUS_thickness_01_
% createPhantoms.m) rather than generating new ones. PRESTUS_thickness_02_
% sim.m already ran an acoustic-only sweep across these phantoms at
% all-mid tissue properties (modules.run_heating_sims=0) - this script
% reruns that same geometry x mid-properties baseline (the 'mid' density
% level here IS that reference point) but WITH heating enabled and with
% skull density additionally varied to bottom/top, so this is a superset
% of that earlier sweep's acoustic-only baseline, not a duplicate of it.
%
% Skull thickness itself is NOT a `parameters` field - it is baked into
% each subject's final_tissues.nii.gz at phantom-creation time (see
% PRESTUS_thickness_01_createPhantoms.m), so varying it here means
% looping over subject_id, not a config field.
%
% CACHING: each of the 36 (thickness, density-level) combinations needs
% its own fresh acoustic AND thermal solve - unlike PRESTUS_paramvar_08/
% 09's acoustic/thermal decoupling trick, there is no shared acoustic
% combo to reuse here, since both thickness (geometry) and density (an
% acoustic-relevant property) change together. All 36 evaluations are
% therefore submitted independently with no phase ordering / wait step.
%
% JOB SUBMISSION: "auto mode", same mechanism as PRESTUS_paramvar_08/09 -
% hpc_detect_system() picks 'slurm'/'qsub' if a scheduler is on PATH, else
% falls back to running locally via 'matlab'; every evaluation is
% submitted individually via prestus_pipeline_start(parameters).
%
% Total: 36 evaluations (3 density levels x 12 thicknesses), well under
% budget. Shortened thermal protocol (20s stim, 0s post-stim) and the same
% storage-minimizing flags as PRESTUS_paramvar_08/09.
%
% Output: data/tussim/<variant_name>/sub-<subject_id>/sub-<subject_id>_
% phantom<affix>.csv per evaluation; figures/paramvar/density_thickness_
% design.csv (design table for PRESTUS_paramvar_12_density_thickness_
% plot.m).

user = 'julkos';

currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
rootpath = fileparts(pathstr);

addpath(fullfile(rootpath, 'code'));
pn = setup_itrusst_paths(rootpath, user);

%% ── Configuration ───────────────────────────────────────────────────────────

setup_name    = 'itrusst_protocol1';
variant_name  = 'itrusst_protocol1a_density_thickness';
thickness_subjects = 301:312;              % thickness = subject_id - 300, in mm
thicknesses_mm      = thickness_subjects - 300;
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
pn.figures = fullfile(rootpath, 'figures', 'paramvar');
if ~exist(pn.figures, 'dir'); mkdir(pn.figures); end
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

% Shortened thermal protocol (see PRESTUS_paramvar_08/09 / doc/PRESTUS_
% itrusst_methods.md, Analysis 4): 20s stimulation, 0s post-stim cool-off.
base_parameters.timing.ptrd = 20;
base_parameters.timing.post_ptri_dur = 0;

base_parameters.hpc.timelimit = '00:30:00';  % see PRESTUS_paramvar_09's pilot run: 00:15:00 was NOT
                                              % a safe margin for the thermal stage at this protocol

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

writetable(design_rows, fullfile(pn.figures, 'density_thickness_design.csv'));
fprintf('Density x thickness design: %d evaluations (%d thicknesses x %d density levels).\n', ...
    n_rows, numel(thickness_subjects), numel(density_levels_name));

%% ── Submit every evaluation (no phase ordering / wait needed - see header) ──

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
