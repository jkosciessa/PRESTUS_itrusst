% CALIBRATION_BENCHMARK_01_SIM
% Compare axisymmetric (2D) vs. full 3D free-water calibration.
%
% For each variant the script runs calibration_setup with the shared
% config_calibration.yaml, overriding only the axisymmetric2D flag and
% the output path.  Wall time is recorded for each variant and saved
% together with the calibrated parameters to a .mat summary file.
%
% Results are consumed by calibration_benchmark_02_plot.m.
%
% Run from the code/ directory or via F5 in the IDE.

clear all; close all; clc;

user = 'julkos';

% get root path (script must be run from file location)
currentFile = mfilename('fullpath');
[pathstr, ~, ~] = fileparts(currentFile);
cd(fullfile(pathstr, '..'))
rootpath = pwd;

addpath(fullfile(rootpath, 'code'));
pn = setup_itrusst_paths(rootpath, user);

%% Variant definitions

variant_names  = {'3D',           'axisymmetric'};
axisym_flags   = [0,               1              ];

%% Equipment library (needed by calibration_setup)

equip_param = load_equipment_config();

%% Base calibration parameters

cd(pn.configs);
parameters_base.calibration = load_parameters('config_calibration.yaml');

if ismac
    parameters_base.startup.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
else
    parameters_base.startup.simnibs_bin_path = pn.simnibs;
end

parameters_base.startup.paths_to_add   = {pn.kwave, pn.minimize};
parameters_base.hpc.ld_library_path    = "/opt/gcc/7.2.0/lib64";
parameters_base.simulation.code_type   = 'cpp_cpu';
parameters_base.simulation.interactive = 0;
parameters_base.io.overwrite_files     = 'always';
parameters_base.platform               = 'slurm';
parameters_base.hpc.wait_for_job       = true;

%% Results container

n_variants = numel(variant_names);
results(n_variants) = struct( ...
    'variant',           '', ...
    'axisymmetric2D',    NaN, ...
    'wall_time_s',       NaN);

%% Run calibration for each variant

for i = 1:n_variants
    parameters = parameters_base;

    parameters.calibration.axisymmetric2D = axisym_flags(i);
    parameters.calibration.path_output = fullfile( ...
        rootpath, 'data', 'calibration_benchmark', variant_names{i});
    parameters.calibration.path_output_profiles = fullfile( ...
        parameters.calibration.path_output, 'PRESTUS_virtual_parameters');

    fprintf('\n=== Variant: %s (axisymmetric2D=%d) ===\n', ...
        variant_names{i}, axisym_flags(i));

    t_start = tic;
    calibration_setup(parameters, equip_param);
    elapsed = toc(t_start);

    results(i).variant        = variant_names{i};
    results(i).axisymmetric2D = axisym_flags(i);
    results(i).wall_time_s    = elapsed;

    fprintf('Finished %s in %.1f s\n', variant_names{i}, elapsed);
end

%% Save summary

out_dir = fullfile(rootpath, 'data', 'calibration_benchmark');
if ~isfolder(out_dir); mkdir(out_dir); end
save(fullfile(out_dir, 'benchmark_results.mat'), 'results', 'variant_names');
fprintf('\nSummary saved to %s\n', fullfile(out_dir, 'benchmark_results.mat'));
